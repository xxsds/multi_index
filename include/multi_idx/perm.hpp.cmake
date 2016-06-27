#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>
#include "sdsl/bit_vectors.hpp"

template<uint8_t t_splitter_bits, uint8_t t_block_errors>
struct splitter_mask{
    static constexpr size_t binomial(size_t n, size_t k){
        if ( n < k ) 
            return 0;
        if ( k == 0 or n == k )
            return 1;
        if ( k == 1 )
            return n;
        return (n*binomial(n-1, k-1))/k;
    }

    static constexpr size_t all_binomial(size_t n, size_t k) {
        return k==0 ? 1 : binomial(n, k) + all_binomial(n, k-1);
    }

    static struct impl {
        std::array<uint64_t, all_binomial(t_splitter_bits, t_block_errors)> data;

        impl(){
            size_t idx=0;
            data[idx++] = 0ULL;
            for (size_t block_errors=1; block_errors <= t_block_errors; ++block_errors) {
                sdsl::bit_vector mask(t_splitter_bits, 0);
                for(size_t i=0; i<block_errors; ++i) mask[mask.size()-i-1] = 1;
                do {
                    data[idx++] = (mask.get_int(0, mask.size())) << (64-mask.size());
                } while ( next_permutation( mask.begin(), mask.end() ) );
            }
        }
    } precomp;
};

template<uint8_t t_splitter_bits, uint8_t t_block_errors>
typename splitter_mask<t_splitter_bits, t_block_errors>::impl splitter_mask<t_splitter_bits,t_block_errors>::precomp;

template<uint8_t t_b, uint8_t t_k>
struct perm{
 	typedef uint64_t (*perm_fun_t)(uint64_t);
	typedef std::vector<perm_fun_t> t_fun_vec;

	static const uint8_t max_dist = 0;
	static const uint8_t match_len = 0;

	static const t_fun_vec mi_permute;
	static const t_fun_vec mi_rev_permute;
	static const std::vector<std::vector<uint8_t>> mi_perms;
	static const std::vector<uint8_t> mi_permute_block_sizes;

	static const t_fun_vec chi_permute;
	static const t_fun_vec chi_rev_permute;
	static const std::vector<std::vector<uint8_t>> chi_perms;
	static const std::vector<uint8_t> chi_permute_block_sizes;
};

