#pragma once

#include <sdsl/int_vector.hpp>

constexpr uint64_t rol(uint64_t x, int rot) {
// Rotate left a complete word.
// x: value to be rotated
// rot: rotate count, negative values will rotate right
// 1 cycle (should be if inlined)

  // A shift by >= bits is undefined by the C/C++ standard.
  // We take care for this case with the "& (bits - 1)" below.
  // For many CPUs this stunt is not neccessary.
  return (x << (rot & (64 - 1))) | (x >> ((-rot) & (64 - 1)));
  }


// Generalized Bit Reversal

const uint64_t a_bfly_mask[]={
  // 0..ld_bits
  // For butterfly ops
  // = all_bits / ((1 << (1 << i)) + 1)
  0x5555555555555555,   // 0
  0x3333333333333333,   // 1
  0x0f0f0f0f0f0f0f0f,   // 2
  0x00ff00ff00ff00ff,   // 3
  0x0000ffff0000ffff,   // 4
  0x00000000ffffffff,   // 5
  0xffffffffffffffff};  // 6

inline uint64_t bit_permute_step_simple(uint64_t x, uint64_t m, uint64_t shift) {
  // INLINE
  // Simplified replacement of bit_permute_step
  // Can always be replaced by bit_permute_step (not vice-versa).
  // x86: >= 5/4 (5/3) cycles
  // ARM: >= 3/2 cycles

    // assert(((m << shift) & m) == 0);
    // assert(((m << shift) >> shift) == m);
    // assert(((m << shift) | m) == all_bits);  // for permutations
    return ((x & m) << shift) | ((x >> shift) & m);
    }


inline uint64_t general_reverse_bits(uint64_t x, uint64_t k) {
// Swap all subwords of given levels.
// See Hacker's Delight, 7.1 "Generalized Bit Reversal"
// k: set of t_subword, i.e. one bit per subword size.

  int i,j;
  uint64_t m;
  for (i = 0; i <= 6-1; ++i) {  // UNROLL
    j = 1 << i;
    if ((k & j) != 0) {
      // x = bit_index_complement(x,j);
      m = a_bfly_mask[i];
      x = bit_permute_step_simple(x,m,j);
      }
    }
  return x;
}

inline uint64_t bswap(uint64_t x) {
  // INLINE
  // Exchange byte order.
  // This can be expressed in assembler:
  // bits = 8: n/a
  // bits = 16: "xchg al,ah" or "rol ax,16"
  // bits = 32: "bswap eax"
  // bits = 64: "bswap rax"
  // bits = 128: "xchg rax,rdx; bswap rax; bswap rdx"

  return general_reverse_bits(x, ~7);
}
