#pragma once

#include <iostream>
#include <algorithm>
#include <vector>
#include "multi_idx/perm.hpp"
#include "sdsl/io.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/sd_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "multi_idx/multi_idx_helper.hpp"

namespace multi_index {
 
// we can match everything with less than t_k errors
  template<uint8_t t_b=4,
           uint8_t t_k=3,
           size_t  t_id=0, // id of the permutation managed by this instance
           typename perm_b_k=perm<t_b,t_b-t_k>,
           typename t_bv=sdsl::bit_vector,
           typename t_sel=typename t_bv::select_1_type> 
  class _simple_buckets_binvector {
    public:
        typedef uint64_t size_type;
        typedef uint64_t entry_type;
        typedef perm_b_k perm;
        enum {id = t_id};

    private:        
        static constexpr uint8_t init_splitter_bits(size_t i=0){
            return i < perm_b_k::match_len ? perm_b_k::mi_permute_block_widths[t_id][t_b-1-i] + init_splitter_bits(i+1) : 0;
        }

        uint64_t              m_n;      // number of items
        sdsl::int_vector<64>  m_entries;
        t_bv                  m_C;     // bit vector for prefix sums of meta-symbols
        t_sel                 m_C_sel; // select1 structure for m_C 

    public:
        static constexpr uint8_t    splitter_bits = init_splitter_bits();

        _simple_buckets_binvector() = default;

        _simple_buckets_binvector(const std::vector<entry_type> &input_entries) {
            std::cout << "Splitter bits " << (uint16_t) splitter_bits << std::endl; 
    
            m_n = input_entries.size();
            m_entries = sdsl::int_vector<64>(input_entries.size(), 0);
            
            build_small_universe(input_entries);
        }

        // assert(errors <= t_k)
        inline std::pair<std::vector<uint64_t>, uint64_t> match(const entry_type q, uint8_t errors=t_k, const bool find_only_candidates=false) {
            uint64_t bucket = get_bucket_id(q);
    
            const auto l = bucket == 0 ? 0 : m_C_sel(bucket) - bucket +1; 
            const auto r = m_C_sel(bucket+1) - (bucket+1) +1;  
            
            auto begin = m_entries.begin() + l;
            auto end = m_entries.begin() + r;
    
            uint64_t candidates = std::distance(begin,end);
            std::vector<entry_type> res;
            
            if(find_only_candidates) return {res, candidates};
            for (auto it = begin; it != end; ++it) {
               if (sdsl::bits::cnt(q^*it) <= errors)
                 res.push_back(*it);
            }
            return {res, candidates};
        }
  
        _simple_buckets_binvector& operator=(const _simple_buckets_binvector& idx) {
            if ( this != &idx ) {
                m_n       = std::move(idx.m_n);
                m_entries   = std::move(idx.m_entries);
                m_C           = std::move(idx.m_C);
                m_C_sel       = std::move(idx.m_C_sel);
                m_C_sel.set_vector(&m_C);
            }
            return *this;
        }

        _simple_buckets_binvector& operator=(_simple_buckets_binvector&& idx) {
            if ( this != &idx ) {
                m_n       = std::move(idx.m_n);
                m_entries   = std::move(idx.m_entries);
                m_C           = std::move(idx.m_C);
                m_C_sel       = std::move(idx.m_C_sel);
                m_C_sel.set_vector(&m_C);
            }
            return *this;
        }

        _simple_buckets_binvector(const _simple_buckets_binvector& idx) {
            *this = idx;
        }

        _simple_buckets_binvector(_simple_buckets_binvector&& idx){
            *this = std::move(idx);
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const {
            using namespace sdsl;
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            uint64_t written_bytes = 0;
            written_bytes += write_member(m_n, out, child, "n");
            written_bytes += m_entries.serialize(out, child, "entries"); 
            written_bytes += m_C.serialize(out, child, "C");
            written_bytes += m_C_sel.serialize(out, child, "C_sel");     
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            using namespace sdsl;
            read_member(m_n, in);
            m_entries.load(in);
            m_C.load(in);
            m_C_sel.load(in, &m_C);
        }

        size_type size() const{
            return m_n;
        }

private:

    inline uint64_t get_bucket_id(const uint64_t x) const {
        return perm_b_k::mi_permute[t_id](x) >> (64-splitter_bits);
    }

    void build_small_universe(const std::vector<entry_type> &input_entries) {
        // countingSort-like strategy to order entries accordingly to bucket_id
        uint64_t splitter_universe = ((uint64_t) 1) << (splitter_bits);

        std::vector<uint64_t> prefix_sums(splitter_universe + 1, 0); 
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
            prefix_sums[i] = sum;
            sum += curr;
        }
        for (auto x : input_entries) {
            uint64_t bucket = get_bucket_id(x);
            m_entries[prefix_sums[bucket]] = x; 
            prefix_sums[bucket]++;
        }
    }
};

  template<typename t_bv=sdsl::bit_vector,
           typename t_sel=typename t_bv::select_1_type> 
  struct simple_buckets_binvector {
      template<uint8_t t_b, uint8_t t_k, size_t t_id, typename t_perm>
      using type = _simple_buckets_binvector<t_b, t_k, t_id, t_perm, t_bv, t_sel>;
  };

}
