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
  template<uint8_t t_b=7,
           uint8_t t_k=3,
           size_t  t_id=0,
           typename perm_b_k=perm<t_b,t_b-t_k>
           >
  class _simple_buckets_vector {
    public:
        typedef uint64_t size_type;
        typedef uint64_t entry_type;
        typedef perm_b_k perm;
        enum {id = t_id};

    private:        
        static constexpr uint8_t init_splitter_bits(size_t i=0){
            return i < perm_b_k::match_len ? perm_b_k::mi_permute_block_widths[t_id][t_b-1-i] + init_splitter_bits(i+1) : 0;
        }
        static constexpr uint8_t    splitter_bits = init_splitter_bits();

        uint64_t              m_n;      // number of items
        sdsl::int_vector<64>  m_entries;
        sdsl::int_vector<64>  m_prefix_sums;

    public:

        _simple_buckets_vector() = default;

        _simple_buckets_vector(const std::vector<entry_type> &input_entries) {
            std::cout << "Splitter bits " << (uint16_t) splitter_bits << std::endl; 
    
            m_n = input_entries.size();
//            check_permutation<_simple_buckets_vector, t_id>(input_entries);
            m_entries = sdsl::int_vector<64>(input_entries.size(), 0);
            
            build_small_universe(input_entries);
        }

        // k with passed to match function
        // assert(k<=t_k)
        inline std::pair<std::vector<uint64_t>, uint64_t> match(const entry_type q, uint8_t errors=t_k, const bool find_only_candidates=false) {
            uint64_t bucket = get_bucket_id(q);
    
            /* DEBUG
            std::cout << "---------------------\n";
            for(auto x: m_perm.mi_perms[t_id])
              std::cout << (uint64_t) x << " ";
            std::cout << std::endl;
            std::cout << "Splitter_bits " << (uint64_t) splitter_bits << std::endl;
            std::bitset<64> a(q);
            uint64_t l = m_perm.mi_permute_block_sizes[0], j = 0; 
            for(size_t i = 0; i < 64; i++) {
              if (i == l) {std::cout << "  |  "; j++; l += m_perm.mi_permute_block_sizes[j];}
              std::cout << (uint64_t) a[i];
            }
            std::cout << std::endl;
                
            std::cout << " q " <<  std::bitset<64>(q) <<std::endl;
            std::cout << "rq " << std::bitset<64>(m_perm.mi_permute[t_id](q)) <<std::endl;
            std::cout << "bq " << std::bitset<64>(bucket) <<std::endl;
            std::cout << "---------------------\n";
            */
            
            const auto l = m_prefix_sums[(bucket)] - bucket; 
            const auto r = m_prefix_sums[bucket+1] - (bucket+1) + 1;  
             
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
  
        _simple_buckets_vector& operator=(const _simple_buckets_vector& idx) {
            if ( this != &idx ) {
                m_n       = std::move(idx.m_n);
                m_entries   = std::move(idx.m_entries);
                m_prefix_sums = std::move(idx.m_prefix_sums);
            }
            return *this;
        }

        _simple_buckets_vector& operator=(_simple_buckets_vector&& idx) {
            if ( this != &idx ) {
                m_n       = std::move(idx.m_n);
                m_entries   = std::move(idx.m_entries);
                m_prefix_sums  = std::move(idx.m_prefix_sums);
            }
            return *this;
        }

        _simple_buckets_vector(const _simple_buckets_vector& idx) {
            *this = idx;
        }

        _simple_buckets_vector(_simple_buckets_vector&& idx){
            *this = std::move(idx);
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const {
            using namespace sdsl;
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            uint64_t written_bytes = 0;
            written_bytes += write_member(m_n, out, child, "n");
            written_bytes += m_entries.serialize(out, child, "entries"); 
            written_bytes += m_prefix_sums.serialize(out, child, "prefix_sums");   
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            using namespace sdsl;
            read_member(m_n, in);
            m_entries.load(in);
            m_prefix_sums.load(in);
        }

        size_type size() const{
            return m_n;
        }

private:

    inline uint64_t get_bucket_id(const uint64_t x) const {
        return perm_b_k::mi_permute[t_id](x) >> (64-splitter_bits); // take the most significant bits
    }

    void build_small_universe(const std::vector<entry_type> &input_entries) {
        // Implement a countingSort-like strategy to order entries accordingly to
        // their splitter_bits MOST significant bits
        // Ranges of keys having the same MSB are not sorted. 
        uint64_t splitter_universe = ((uint64_t) 1) << (splitter_bits);

        std::vector<uint64_t> prefix_sums(splitter_universe + 1, 0); // includes a sentinel
        m_prefix_sums = sdsl::int_vector<64>(prefix_sums.size(), 0);
        
        for (auto x: input_entries) {
            prefix_sums[get_bucket_id(x)]++;
        }

        uint64_t sum = prefix_sums[0];
        prefix_sums[0] = 0;
        m_prefix_sums[0] = prefix_sums[0];
        for(uint64_t i = 1; i < prefix_sums.size(); ++i) {
            uint64_t curr = prefix_sums[i];
            prefix_sums[i] = sum + i; // +i is to obtain a striclty monotone sequence as we would have with binary vectors
            m_prefix_sums[i] = prefix_sums[i];
            sum += curr;
        }

        // Partition elements into buckets accordingly to their less significant bits
        for (auto x : input_entries) {
            uint64_t bucket = get_bucket_id(x);
            m_entries[prefix_sums[bucket]-bucket] = x; // -bucket is because we have a striclty monotone sequence 
            prefix_sums[bucket]++;
        }
    }
};

    struct simple_buckets_vector {
        template<uint8_t t_b, uint8_t t_k, uint8_t t_id, typename t_perm>
        using type = _simple_buckets_vector<t_b, t_k, t_id, t_perm>;
    };

}
