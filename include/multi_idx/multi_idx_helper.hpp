#pragma once

#include <unordered_map>
#include "multi_idx/tuple_foreach.hpp"

namespace multi_index {

static std::unordered_map<std::string, uint64_t> stat_counter;
//static std::map<uint64_t, uint64_t> bucket_distr;
//static std::map<uint64_t, uint64_t> cluster_distr;

template<typename t_index>
std::map<uint64_t, uint64_t>
get_bucket_dist(const t_index& index){
    std::map<uint64_t, uint64_t> res;
    uint64_t buckets = (typename decltype(index.m_C)::rank_1_type){&(index.m_C)}(index.m_C.size());
    for(uint64_t bucket=0; bucket<buckets; ++bucket){
        const auto l = bucket == 0 ? 0 : index.m_C_sel(bucket) - bucket +1; 
        const auto r = index.m_C_sel(bucket+1) - (bucket+1) + 1;  
        ++res[r-l];            
    }
    return res;
}

struct bucket_accumulator {
    std::map<uint64_t, uint64_t> cnt;
    template <typename T>
    void operator()(T&& t, std::size_t) {
        auto res = get_bucket_dist(t);
        for(auto x : res){
            cnt[x.first] += x.second;
        }
    }
};

template<typename t_index>
std::map<uint64_t, uint64_t>
get_cluster_dist(const t_index& index){
    std::map<uint64_t, uint64_t> res;
    for(auto it = index.m_first_level.begin(); it+3 < index.m_first_level.end(); it+=3){
        ++res[(*(it+3)) - (*it)];
    }
    return res;
}

struct cluster_accumulator {
    std::map<uint64_t, uint64_t> cnt;
    template <typename T>
    void operator()(T&& t, std::size_t) {
        auto res = get_cluster_dist(t);
        for(auto x : res){
            cnt[x.first] += x.second;
        }
    }
};


template<typename t_index>
std::map<uint64_t, uint64_t>
get_error_dist(const t_index& index){
    std::map<uint64_t, uint64_t> res;
    for(auto it = index.m_first_level.begin()+2; it < index.m_first_level.end(); it+=3){
        ++res[*it];
    }
    return res;
}

struct error_accumulator {
    std::map<uint64_t, uint64_t> cnt;
    template <typename T>
    void operator()(T&& t, std::size_t) {
        auto res = get_error_dist(t);
        for(auto x : res){
            cnt[x.first] += x.second;
        }
    }
};


template<typename t_accum, typename t_index>
std::map<uint64_t, uint64_t>
get_dist(const t_index& index){
    t_accum acc;
    tuple_foreach(index.m_idx, acc);
    return acc.cnt;
}


// Helper structs to generate a tuple of strategy classes
template<size_t t_idx_id, uint8_t t_b, uint8_t t_k, class t_idx_strategy, typename t_perm=perm<t_b,t_b-t_k>>
struct perm_type_gen{
    using perm_type = std::tuple<typename t_idx_strategy::template type<t_b, t_k, t_idx_id-1, t_perm>>;
    using type = decltype(std::tuple_cat(perm_type(), std::declval<typename perm_type_gen<t_idx_id-1,t_b,t_k,t_idx_strategy,t_perm>::type>())); 
};

template<uint8_t t_b, uint8_t t_k, class t_idx_strategy, typename t_perm>
struct perm_type_gen<0, t_b, t_k, t_idx_strategy, t_perm>{
    using type = std::tuple<>;
};

// Trait for the type of mid_entries for the split strategy classes
template<uint8_t t_w>
    struct mid_entries_trait{
    using type = sdsl::int_vector<>;
    static type get_instance(typename type::size_type n, typename type::value_type x) {
        return type(n, x, t_w);
    } 
    static void assign(type& rac, uint64_t i, uint64_t x) { rac[i] = x; }
};

template<>
    struct mid_entries_trait<16>{
    using type = sdsl::int_vector<16>;
    static type get_instance(typename type::size_type n, typename type::value_type x) {
        return type(n, x);
    } 
    static void assign(type& rac, uint64_t i, uint64_t x) { rac[i] = x; }
};

template<>
struct mid_entries_trait<8>{
    using type = sdsl::int_vector<8>;
    static type get_instance(typename type::size_type n, typename type::value_type x) {
        return type(n, x);
    } 
    static void assign(type& rac, uint64_t i, uint64_t x) { rac[i] = x; }
};

template<>
struct mid_entries_trait<0>{
    struct type : public sdsl::int_vector<>{
        uint64_t operator[](uint64_t) const{
            return 0;
        }
    };

    static type get_instance(typename type::size_type n, typename type::value_type x) {
        return type{};
    }
    static void assign(type&, uint64_t, uint64_t) { }
};

template<typename t_strat, size_t t_id>
void check_permutation(const std::vector<typename t_strat::entry_type> &input_entries) {
    std::cout << "Check permuting functions\n";
    for(auto x : input_entries) {
        if (x != t_strat::perm_b_k::mi_rev_permute[t_id](t_strat::perm_b_k::mi_permute[t_id](x)))
            std::cout << "ERROR permuting " << x << std::endl; 
    }
}

}
