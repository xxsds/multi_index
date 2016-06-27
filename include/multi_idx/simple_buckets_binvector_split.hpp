#pragma once

#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <bitset>
#include "multi_idx/perm.hpp"
#include "sdsl/io.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/sd_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "multi_idx/multi_idx_helper.hpp"
#include "simd_utils.hpp"

#define LIKELY(x)   (__builtin_expect((x), 1))
#define UNLIKELY(x) (__builtin_expect((x), 0))

namespace multi_index {

  template<uint8_t t_b,
           uint8_t t_k,
           uint8_t t_id, // id of the permutation managed by this instance
           typename perm_b_k,
           typename t_bv,
           typename t_sel,
           bool use_simd> 
  class _simple_buckets_binvector_split_common {
    public:
        typedef uint64_t size_type;
        typedef uint64_t entry_type;
        typedef perm_b_k perm;
        enum {id = t_id};

    protected:  
        static constexpr uint8_t init_splitter_bits(size_t i=0){
            return i < perm_b_k::match_len ? perm_b_k::mi_permute_block_widths[t_id][t_b-1-i] + init_splitter_bits(i+1) : 0;
        }

        /* Low_* stuff control how many of the less signifigant bits form the lower part */
        static constexpr uint8_t    low_bits    = 32; // PLS, keep this a power of 2, better if word aligned
        static constexpr uint64_t   low_mask    = (1ULL<<low_bits)-1;
        static constexpr uint8_t    splitter_bits = init_splitter_bits();
        static constexpr uint8_t    mid_bits = 64 - (low_bits + splitter_bits);
        static constexpr uint8_t    mid_shift   = low_bits; 
        static constexpr uint64_t   mid_mask = (1ULL<<mid_bits)-1;
        static constexpr uint8_t    high_shift = (64-splitter_bits);
        using  mid_entries_type = typename mid_entries_trait<mid_bits>::type; 

        uint64_t                    m_n;      // number of items
        sdsl::int_vector<low_bits>  m_low_entries;
        mid_entries_type            m_mid_entries; 
        //sdsl::int_vector<64>  fake_entries; for DEBUG USE ONLY
        t_bv                        m_C;     // bit vector for prefix sums of meta-symbols
        t_sel                       m_C_sel; // select1 structure for m_C

    public:
        _simple_buckets_binvector_split_common() = default;

        _simple_buckets_binvector_split_common(const std::vector<entry_type> &input_entries){
            std::cout << "Splitter bits " << (uint16_t) splitter_bits << std::endl; 
            m_n = input_entries.size();
            m_low_entries = sdsl::int_vector<low_bits>(input_entries.size(), 0); 
            m_mid_entries = mid_entries_trait<mid_bits>::get_instance(input_entries.size(), 0);
            build_small_universe(input_entries);
        }

        inline std::pair<std::vector<uint64_t>, uint64_t> match(const entry_type q, uint8_t errors=t_k, const bool find_only_candidates=false) const {
          
            const uint64_t bucket = get_bucket_id(q);
            
            auto l = bucket == 0 ? 0 : m_C_sel(bucket) - bucket +1; 
            const auto r = m_C_sel(bucket+1) - (bucket+1) + 1;  
    
           // std::cout << "q " << q << " b " << bucket << " l " << l << " r " <<  r << std::endl;
    
            const uint64_t candidates = r-l;
            std::vector<entry_type> res;
            
            if(find_only_candidates) return {res, candidates};
            if(errors >= 6) res.reserve(128);
  
            const uint64_t q_permuted   = perm_b_k::mi_permute[t_id](q); 
            const uint64_t q_high       = (q_permuted>>(high_shift))<<high_shift;
            const uint64_t q_low        = q_permuted & low_mask;
                                
            const auto begin  = m_low_entries.begin() + l;
            const auto end    = m_low_entries.begin() + r;

            if ( use_simd ) {
                auto it = begin;
                for (; ((size_t)it)%16 != 0 and it != end; ++it, ++l) {
                  const uint32_t item_low = ((uint64_t) *it);   
                  if UNLIKELY(_mm_popcnt_u32(q_low^item_low) <= errors) {
                    const uint64_t curr_el = q_high | (((uint64_t) m_mid_entries[l]) << mid_shift) | item_low;
                    if (_mm_popcnt_u64(q_permuted^curr_el) <= errors)
                      res.push_back(perm_b_k::mi_rev_permute[t_id](curr_el));
                  }
                }
                
                const __m128i query = _mm_set1_epi32(q_low);
                const __m128i tk    = _mm_set1_epi32(errors+1);

                for(; it < end-4; it+=4, l+=4) {
                  _mm_prefetch((const char *) (it+4), _MM_HINT_T0);
                  const __m128i vec     = _mm_xor_si128(_mm_load_si128(reinterpret_cast<const __m128i*>(it)), query); // need to be aligned!
                  const __m128i popcounts = popcount_epi32(vec); // not an intrinsics, see simd_utils.hpp
     
                  uint16_t mask = (_mm_movemask_epi8(_mm_cmpgt_epi32(tk, popcounts))) & 0x1111;

                  while(UNLIKELY(mask)) {
                    size_t i =  (__builtin_ffs(mask))-1;
                    mask = mask ^ (1 << i);
                    i = i/4;
                    const uint64_t curr_el = q_high | (((uint64_t) m_mid_entries[l+i]) << mid_shift) | it[i];
                    if (_mm_popcnt_u64(q_permuted^curr_el) <= errors)
                      res.push_back(perm_b_k::mi_rev_permute[t_id](curr_el));
                  }          
                }
                
                for (; it != end; ++it, ++l) {
                  const uint32_t item_low = ((uint64_t) *it);   
                  if UNLIKELY(_mm_popcnt_u32(q_low^item_low) <= errors) {
                    const uint64_t curr_el = q_high | (((uint64_t) m_mid_entries[l]) << mid_shift) | item_low;
                    if (_mm_popcnt_u64(q_permuted^curr_el) <= errors)
                      res.push_back(perm_b_k::mi_rev_permute[t_id](curr_el));
                  }
                }
            } else {
                for (auto it = begin; it != end; ++it, ++l) {
                   const uint64_t item_low = ((uint64_t) *it);
                   // uint64_t curr_el = q_high | (((uint64_t) m_mid_entries[l]) << mid_shift) | item_low;
                   // if(curr_el != fake_entries[i]){
                      // std::cout <<std::endl;
                      // std::cout << "q  " << std::bitset<64>(q_high).to_string() << std::endl;
                      // std::cout << "mi " << std::bitset<64>(m_mid_entries[l]).to_string() << std::endl;
                      // std::cout << "lo " << std::bitset<64>(item_low).to_string() << std::endl;
                      // std::cout << "cu " << std::bitset<64>(curr_el).to_string() << std::endl;
                      //std::cout << "fe " << std::bitset<64>(fake_entries[i]).to_string() << std::endl;
                   // }
                
                   if (sdsl::bits::cnt(q_low^item_low) <= errors) {
                     const uint64_t curr_el = q_high | (((uint64_t) m_mid_entries[l]) << mid_shift) | item_low;
                     if (sdsl::bits::cnt(q_permuted^curr_el) <= errors)
                       res.push_back(perm_b_k::mi_rev_permute[t_id](curr_el));
                   }
                }
            }
            return {res, candidates};
        }
  
        _simple_buckets_binvector_split_common& operator=(const _simple_buckets_binvector_split_common& idx) {
            if ( this != &idx ) {
                m_n       = std::move(idx.m_n);
                m_low_entries   = std::move(idx.m_low_entries);
                m_mid_entries   = std::move(idx.m_mid_entries);
                //fake_entries   = std::move(idx.fake_entries);           
                m_C           = std::move(idx.m_C);
                m_C_sel       = std::move(idx.m_C_sel);
                m_C_sel.set_vector(&m_C);
            }
            return *this;
        }

        _simple_buckets_binvector_split_common& operator=(_simple_buckets_binvector_split_common&& idx) {
            if ( this != &idx ) {
                m_n       = std::move(idx.m_n);
                m_low_entries   = std::move(idx.m_low_entries);
                m_mid_entries   = std::move(idx.m_mid_entries);
                //fake_entries   = std::move(idx.fake_entries);                               
                m_C           = std::move(idx.m_C);
                m_C_sel       = std::move(idx.m_C_sel);
                m_C_sel.set_vector(&m_C);
            }
            return *this;
        }

        _simple_buckets_binvector_split_common(const _simple_buckets_binvector_split_common& idx) {
            *this = idx;
        }

        _simple_buckets_binvector_split_common(_simple_buckets_binvector_split_common&& idx){
            *this = std::move(idx);
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const {
            using namespace sdsl;
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            uint64_t written_bytes = 0;
            written_bytes += write_member(m_n, out, child, "n");
            written_bytes += m_low_entries.serialize(out, child, "low_entries"); 
            written_bytes += m_mid_entries.serialize(out, child, "mid_entries"); 
            written_bytes += m_C.serialize(out, child, "C");
            written_bytes += m_C_sel.serialize(out, child, "C_sel");  
            //written_bytes += fake_entries.serialize(out, child, "fake");  
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            using namespace sdsl;
            read_member(m_n, in);
            m_low_entries.load(in);
            m_mid_entries.load(in);
            m_C.load(in);
            m_C_sel.load(in, &m_C);
            //fake_entries.load(in);
        }

        size_type size() const{
            return m_n;
        }

protected:

    inline uint64_t get_bucket_id(const uint64_t x) const {
        return perm_b_k::mi_permute[t_id](x) >> (64-splitter_bits);
    }

    void build_small_universe(const std::vector<entry_type> &input_entries) {
        // Implement a countingSort-like strategy to order entries accordingly to
        // their splitter_bits MOST significant bits
        // Ranges of keys having the same MSB are not sorted. 
        uint64_t splitter_universe = ((uint64_t) 1) << (splitter_bits);

        std::vector<uint64_t> prefix_sums(splitter_universe + 1, 0); // includes a sentinel
        for (auto x: input_entries) {
            prefix_sums[get_bucket_id(x)]++;
        }
        
        m_C = t_bv(splitter_universe+input_entries.size(), 0);
        size_t idx = 0;
        for(auto x : prefix_sums) {         
          for(size_t i = 0; i < x; ++i, ++idx)
            m_C[idx] = 0; 
          m_C[idx++] = 1;
        }
        m_C_sel = t_sel(&m_C);
        
        uint64_t sum = prefix_sums[0];
        prefix_sums[0] = 0;
        for(uint64_t i = 1; i < prefix_sums.size(); ++i) {
            uint64_t curr = prefix_sums[i];
            prefix_sums[i] = sum + i; // +i is to obtain a striclty monotone sequence as we would have with binary vectors
            sum += curr;
        }

        // Partition elements into buckets accordingly to their less significant bits
        for (auto x : input_entries) {
            uint64_t bucket = get_bucket_id(x);
            uint64_t permuted_item = perm_b_k::mi_permute[t_id](x);
            //fake_entries[prefix_sums[bucket]-bucket]   = permuted_item;
            //m_mid_entries[prefix_sums[bucket]-bucket]  = (permuted_item>>mid_shift) & mid_mask; // -bucket is because we have a striclty monotone sequence 
            mid_entries_trait<mid_bits>::assign(m_mid_entries, prefix_sums[bucket]-bucket, (permuted_item>>mid_shift) & mid_mask); 
            m_low_entries[prefix_sums[bucket]-bucket] = (permuted_item & low_mask);
            prefix_sums[bucket]++;
        }
    }
};


template<uint8_t t_b=4,
       uint8_t t_k=3,
       uint8_t t_id=0, // id of the permutation managed by this instance
       typename perm_b_k=perm<t_b,t_b-t_k>,
       typename t_bv=sdsl::bit_vector,
       typename t_sel=typename t_bv::select_1_type,
       bool use_simd=false> 
class _simple_buckets_binvector_split : public _simple_buckets_binvector_split_common<t_b, t_k, t_id, perm_b_k, t_bv, t_sel, use_simd>{
    typedef _simple_buckets_binvector_split_common<t_b,t_k,t_id,perm_b_k,t_bv,t_sel,use_simd> base;
    using base::base;
};


template<typename t_bv=sdsl::bit_vector,
        bool use_simd=false,
       typename t_sel=typename t_bv::select_1_type> 
struct simple_buckets_binvector_split {
template<uint8_t t_b, uint8_t t_k, uint8_t t_id, typename t_perm>
using type = _simple_buckets_binvector_split<t_b, t_k, t_id, t_perm, t_bv, t_sel, use_simd>;
};
}
