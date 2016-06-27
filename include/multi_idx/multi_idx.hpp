#pragma once

#include <algorithm>
#include <vector>
#include <future>
#include "sdsl/io.hpp"
#include "multi_idx/perm.hpp"
#include "multi_idx/tuple_foreach.hpp"
#include "multi_idx/multi_idx_helper.hpp"
#include "multi_idx/simple_buckets_binsearch.hpp"
#include "multi_idx/simple_buckets_binvector.hpp"
#include "multi_idx/simple_buckets_vector.hpp"
#include "multi_idx/simple_buckets_binvector_unaligned.hpp"
#include "multi_idx/simple_buckets_binvector_split.hpp"
#include "multi_idx/simple_buckets_binvector_split_xor.hpp"
#include "multi_idx/triangle_buckets_binvector_split_simd.hpp"
#include "multi_idx/triangle_clusters_binvector_split.hpp"
#include "multi_idx/triangle_clusters_binvector_split_threshold.hpp"
#include "multi_idx/xor_buckets_binvector_split.hpp"

namespace multi_index {

/*! Basic multi-index implementation.
 *  \tparam t_idx_strategy Index class for matching under a fixed permutation.
 *  \tparam t_k            Number of allowed errors.
 *  \tparam t_b            Keys are split into t_b blocks.
 *
 *  \par The default value of blocks t_b is t_k+1. This way, a key K in the
 *       index has at least one block in common with the query key Q if
 *       H(K,Q) <= t_k. It is also possible to specify t_b > t_k+1. In general
 *       a key K with H(K,Q) has at least t_b-t_k blocks in common with the
 *       query key Q.
 *  \sa multi_index_red
 */
template<typename t_idx_strategy,
           uint8_t t_k = 3,
           uint8_t t_b = t_k+1
           >
class multi_idx {
    public:    
        typedef uint64_t          size_type;
        typedef perm<t_b,t_b-t_k> perm_b_k;
    private:
        static constexpr size_t m_num_perms = std::tuple_size<decltype(perm_b_k::mi_perms)>::value;
        typename perm_type_gen<m_num_perms, t_b, t_k, t_idx_strategy>::type m_idx;

    public:
        multi_idx() = default;
        multi_idx(const multi_idx &) = default;
        multi_idx(multi_idx &&) = default;
        multi_idx &operator=(const multi_idx &) = default;
        multi_idx &operator=(multi_idx &&) = default;

        /*!
        *  \param keys  Vector of hash values
        *  \pre Items are all different (no duplicates)
        */
        multi_idx(const std::vector<uint64_t>& keys, bool async=false) {
            constructor c{keys, async};
            tuple_foreach(m_idx, c);
        }

        std::pair<std::vector<uint64_t>,uint64_t> match(const uint64_t query, const bool find_only_candidates=false) {
            std::vector<uint64_t> matches;
            uint64_t candidates = 0;
            matcher m{matches, candidates, query, find_only_candidates};
            tuple_foreach(m_idx, m);
            return {matches, candidates};
        }

        //! Serializes the data structure into the given ostream
        uint64_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                          std::string name = "") const {
            using namespace sdsl;
            structure_tree_node *child =
                structure_tree::add_child(v, name, util::class_name(*this));
            uint64_t written_bytes = 0;
            serializer s{written_bytes, child, out};
            tuple_foreach(m_idx, s);
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream &in) {
            loader l{in};
            tuple_foreach(m_idx, l);
        }

        uint64_t size() const {
            if (m_num_perms > 0) {
                return std::get<0>(m_idx).size();
            }
            return 0;
        }

    private:
        // Functors which do the actual work on the tuple of indexes

        struct constructor {
            const std::vector<uint64_t>& items;
            bool is_async;
            std::vector<std::future<int>> futures;

            constructor(const std::vector<uint64_t>& f_items, bool asy=false):items(f_items),is_async(asy){};
            template <typename T>
            void operator()(T&& t, std::size_t i)  {
                if ( is_async ) {
                    std::cout<<"start construction thread"<<std::endl;
                    futures.push_back( std::async(std::launch::async,[&](){ 
                        t = (typename std::remove_reference<T>::type){items}; 
                        return 1;
                    }));
                } else {
                    t = (typename std::remove_reference<T>::type){items}; 
                }
            }
        };

        struct serializer {
            uint64_t& wb;
            sdsl::structure_tree_node *c;
            std::ostream& o;
            serializer(uint64_t& written_bytes, sdsl::structure_tree_node *child, std::ostream &out) : 
                wb(written_bytes), c(child), o(out) {};
            template <typename T>
            void operator()(T&& t, std::size_t i) const {
                wb += sdsl::serialize(t, o, c, "multi_idx");
            }
        };

        struct loader {
            std::istream& in;
            loader(std::istream &f_in) : in(f_in) {};
            template <typename T>
            void operator()(T&& t, std::size_t i) const {
                sdsl::load(t, in);
            }
        };

        struct matcher {
            std::vector<uint64_t>& matches;
            uint64_t& candidates;
            uint64_t query;
            bool only_cands;
            matcher(std::vector<uint64_t>& mats, uint64_t& cands, uint64_t qry, bool only_cand) : 
                matches(mats), candidates(cands), query(qry), only_cands(only_cand) 
            { };

            template<typename T>
            void operator()(T&& t, std::size_t i) const {
                auto res = t.match(query, t_k, only_cands);
                matches.insert(matches.end(), std::get<0>(res).begin(), std::get<0>(res).end());    
                candidates += std::get<1>(res);
            }
        };

};

}
