#pragma once

#include <algorithm>
#include <vector>
#include <future>
#include "sdsl/io.hpp"
#include "multi_idx/perm.hpp"
#include "multi_idx/multi_idx.hpp"

namespace multi_index {

/*!  multi index-index implementation
 *  \tparam t_idx_strategy Index class for matching in a fixed permutation.
 *  \tparam t_k            Number of allowed errors.
 *  \tparam t_b            The hash is split into t_b blocks.
 *  \sa multi_idx
 */
template<typename t_idx_strategy,
         uint8_t t_k=3,
         uint8_t t_block_errors=1>
class multi_idx_red {
    static_assert(t_k >= t_block_errors,"It should hold that t_k >= t_block_errors");
    public:    
        static constexpr uint8_t t_b = (t_k/(t_block_errors+1)) + 1;
        typedef uint64_t  size_type;
        typedef perm<t_b,1> perm_b_k;

    private:
        static constexpr size_t m_num_perms = std::tuple_size<decltype(perm_b_k::mi_perms)>::value;
    public:
        typename perm_type_gen<m_num_perms, t_b, 1, t_idx_strategy, perm<t_b,1>>::type m_idx;

    public:
        multi_idx_red() = default;
        multi_idx_red(const multi_idx_red &) = default;
        multi_idx_red(multi_idx_red &&) = default;
        multi_idx_red &operator=(const multi_idx_red &) = default;
        multi_idx_red &operator=(multi_idx_red &&) = default;

        /*!
        *  \param keys  Vector of hash values
        *  \pre Items are all different (no duplicates)
        */
        multi_idx_red(const std::vector<uint64_t>& keys, bool async=false) {
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

            constructor(const std::vector<uint64_t>& f_items, bool asy=false):items(f_items),is_async(asy) {};

            template <typename T>
            void operator()(T&& t, std::size_t i) {
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
                wb += sdsl::serialize(t, o, c, "multi_idx_red");
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
                using TT = typename std::remove_reference<T>::type;
                  // For all block_errors <= t_block_errors match
                  // with flipping block_errors bits for t_k errors
                  for (auto block_mask : splitter_mask<TT::splitter_bits, t_block_errors>::precomp.data) {
                        uint64_t query_flipped = TT::perm::mi_permute[TT::id](query);
                        query_flipped = query_flipped ^ block_mask;
                        query_flipped = TT::perm::mi_rev_permute[TT::id](query_flipped);
                        uint32_t block_errors = sdsl::bits::cnt(block_mask);
                        auto res = t.match(query_flipped, t_k-block_errors, only_cands);
                        matches.insert(matches.end(), std::get<0>(res).begin(), std::get<0>(res).end());    
                        candidates += std::get<1>(res);
                  }
            }
        };

};

}
