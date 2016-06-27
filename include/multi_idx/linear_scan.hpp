#pragma once

#include <algorithm>
#include <vector>
#include "sdsl/io.hpp"

namespace multi_index {

template<uint8_t t_k=3>
class linear_scan {
    public:    
        typedef uint64_t size_type;
    private:
        std::vector<uint64_t> m_keys;

    public:
        linear_scan() = default;
        linear_scan(const linear_scan &) = default;
        linear_scan(linear_scan &&) = default;
        linear_scan &operator=(const linear_scan &) = default;
        linear_scan &operator=(linear_scan &&) = default;

        /*!
        *  \param keys  Vector of hash values
        *  \pre Items are all different (no duplicates)
        */
        linear_scan(const std::vector<uint64_t>& keys, bool async=false) {
            m_keys = keys;  
        }

        std::pair<std::vector<uint64_t>,uint64_t> match(const uint64_t query, const bool find_only_candidates=false) {
            std::vector<uint64_t> matches;
            uint64_t candidates = 0;
            for(auto it = m_keys.begin(); it!=m_keys.end(); ++it){
                if ( sdsl::bits::cnt((*it)^query) < t_k ) {
                    matches.push_back(*it);
                }
            }
            candidates = m_keys.size();
            return {matches, candidates};
        }

        //! Serializes the data structure into the given ostream
        uint64_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                          std::string name = "") const {
            using namespace sdsl;
            structure_tree_node *child =
                structure_tree::add_child(v, name, util::class_name(*this));
            structure_tree::add_size(child, 0);
            return 0;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream &in) {
            uint64_t x;
            m_keys.clear();
            while(in.read((char*) &x,sizeof(x))){
                m_keys.push_back(x);
            }
            std::cout<<"loaded "<<m_keys.size()<<" keys"<<std::endl;
        }

        uint64_t size() const {
           return m_keys.size();
        }
};

}
