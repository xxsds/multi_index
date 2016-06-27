
//////
// Intro

// Some bit hacks and permutations

// (c) 2011..2014 by Jasper L. Neumann
// www.sirrida.de / programming.sirrida.de
// E-Mail: info@sirrida.de

// Granted to the public domain
// First version: 2011-02
// Last change: 2014-09-25

// This is a collection of some bit fiddling procedures.
// All routines act on a word size depending on the
// steering file included before (perm_b*.c).
// These routines naturally adapt to several word (working) sizes
// of a bit length of a power of 2 (as long as the compiler allows).

// Auxiliary routines
//   - odd
//   - gcd
//   - mul_inv
//   - rol
//   - rol_lo
//   - rolc_lo
//   - gray_code
//   - inv_gray_code
//   - nr_1bits
//   - nr_leading_0bits
//   - nr_trailing_0bits
//   - is_contiguous_1bits
//   - tbm
//   - blend
//   - simd_odd
//   - bit_permute_step
//   - bit_permute_step_simple
//   - bit_permute_step2
//   - identity_perm
//   - invert_perm
//   - random_perm
//   - used_source_bits
//   - used_target_bits
// Bit index operations
//   - bit_index_complement
//   - bit_index_swap
//   - bit_index_swap_complement
//   - bit_index_ror
//   - transpose
//   - shuffle_power
//   - unshuffle_power
//   - permute_bpc
//   - invert_bpc
// Generalized Bit Reversal
//   - general_reverse_bits
//   - bswap
// Swap by primitives
//   - prim_swap
// Shuffle and unshuffle
//   - shuffle
//   - unshuffle
// Compress and expand
//   * Compress bit masks
//     - compress_mask_right
//     - compress_mask_left
//     - compress_mask
//   * Generate configuration
//     - gen_ce_right
//     - gen_ce_left
//   * Usage
//     - apply_compress_right
//     - apply_compress_left
//     - apply_expand_right
//     - apply_expand_left
//   * Compound
//     - compress_right
//     - compress_left
//     - compress
//     - expand_right
//     - expand_left
//     - expand
// Butterfly network
//   - butterfly
//   - bfly
//   - ibfly
//   - bfly_parity
// Rotate via butterfly
//   * Generate configuration
//     - gen_frot
//     - gen_vrot
//   * Compound
//     - fror_bfly
//     - frol_bfly
//     - frot_bfly
//     - vror_bfly
//     - vrol_bfly
//     - vrot_bfly
// SWAR rotate
//   * frol
//   * fror
//   * frot
//   * frolc
//   * frorc
//   * frotc
//   * vrol
//   * vror
//   * vrot
// Compress/expand-flip via butterfly
//   * Generate configuration
//     - gen_cef_right
//     - gen_cef_left
//   * Compound
//     - compress_flip_right
//     - compress_flip_left
//     - compress_flip
//     - expand_flip_right
//     - expand_flip_left
//     - expand_flip
// Omega/flip
//   - omega
//   - flip
// Permutations via Benes network
//   - gen_benes_ex
//   - gen_benes
//   - benes_fwd
//   - benes_bwd
//   - benes_fwd_ex
//   - benes_bwd_ex
//   - benes_parity

#include "perm_bas.h"


//////
// Random replacement

mycall t_bits random_bits() {
// rand() only gives 15 bits for gcc.
// Replace with your own generator if you want.

  // t_bits x;
  // x = rand();
  // x = (x << 12) + rand();
  // x = (x << 12) + rand();
  // return x;

  my_randseed = my_randseed * 0x08088405 + 1;  // random as in Delphi (32 bit)
    // this alternates between even and odd numbers
  return my_randseed + ((t_bits)(rand() & 32767));  // we correct this
  }


//////
// Auxiliary stuff

mycall t_bool odd(t_int x) {
// Is x an odd number?
// INLINE

  return (x & 1) != 0;
  }

mycall t_longint gcd(t_longint a, t_longint b) {
// Greatest common divisor.

  t_longint t;

  // a = abs(a);
  // b = abs(b);
  while (b != 0) {
    t = a % b;
    a = b;
    b = t;
    }
  return a;
  }

mycall t_bits mul_inv(t_bits x) {
// Calculate multiplicative inverse, i.e. mul_inv(x)*x==1.
// The result is defined for all odd values of x.
// For even x we simply return 0.
// See Hacker's Delight, 10.16 "Exact division by constants"
// Multiplicative inverse modulo 2**bits by Newton's method.

  t_bits xn,t;

  if (!odd(x))
    return 0;  // only defined for odd numbers

  xn = x;
  while (true) {
    t = x * xn;
    if (t == 1)
      return xn;
    xn = xn * (2 - t);
    }
  }

mycall t_bits rol(t_bits x, t_int rot) {
// INLINE
// Rotate left a complete word.
// x: value to be rotated
// rot: rotate count, negative values will rotate right
// 1 cycle (should be if inlined)

  // A shift by >= bits is undefined by the C/C++ standard.
  // We take care for this case with the "& (bits - 1)" below.
  // For many CPUs this stunt is not neccessary.
  return (x << (rot & (bits - 1))) | (x >> ((-rot) & (bits - 1)));
  }

mycall t_bits rol_lo(t_bits x, t_int rot, t_subword sw) {
// Rotate left. This is not a SWAR operation.
// x: value to be rotated
// rot: rotate count, negative values will rotate right
// sw: log_2 of the number of affected bits, must be <= ld_bits
// Only the low (1 << sw) bits of x are processed and returned.

  t_int b;  // # affected bits
  t_int r;  // rot % b
  t_bits m;  // mask for affected bits

  if (sw >= ld_bits) {
    // Prevent shifting by >= bits at [*].
    return rol(x, rot);
    }
  else {
    b = 1 << sw;

    // m = a_sw_base[sw];
    m = (lo_bit << b) - 1;  // [*] m must become -1 for sw=bits

    r = rot & (b - 1);
    x = x & m;
    x = (x << r) | (x >> (b - r));  // [*]
    return x & m;
    }
  }

mycall t_bits rolc_lo(t_bits x, t_int rot, t_subword sw) {
// Rotate left and complement. This is not a SWAR operation.
// x: value to be rotated
// rot: rotate count, negative values will rotate/complement right
// sw: log_2 of the number of affected bits, must be <= ld_bits
// Only the low (1 << sw) bits of x are processed and returned.

  t_int b;  // # affected bits
  t_int r;  // rot % b
  t_bits m;  // mask for affected bits

  b = 1 << sw;

  m = a_sw_base[sw];  // m must become all_bits for sw=bits

  r = rot & (b-1);
  x = x & m;

  if (b-r >= bits) {
    // Prevent shifting by b-r >= bits.
    x = (x << r) ^ ((lo_bit << r)-1);
    }
  else {
    // x = (x << r) | ((x >> (b-r)) ^ ((lo_bit << r) - 1));  // same result
    x = ((x << r) | (x >> (b-r))) ^ ((lo_bit << r) - 1);
    // x = rol_lo(x,r) ^ ((lo_bit << r) - 1);
    }

  if ((rot & b) != 0) {
    x = ~x;
    }
  return x & m;
  }

mycall t_bits gray_code(t_bits x) {
// See Hacker's Delight, 13.1 "Gray Code"

  return x ^ (x >> 1);
  }

mycall t_bits inv_gray_code(t_bits x) {
// See Hacker's Delight, 13.1 "Gray Code"

  t_int i;

  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    x = x ^ (x >> (1 << i));
    }
  return x;
  }

mycall t_int nr_1bits(t_bits x) {
// Count the set bits, aka population count.
// See Hacker's Delight, 5.1 "Counting 1-Bits"

  t_int i;

  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    x = (x & a_bfly_mask[i]) + ((x >> (1 << i)) & a_bfly_mask[i]);
    }
  return (t_int)(x);
  }

mycall t_int nr_leading_0bits(t_bits x) {
// See Hacker's Delight, 5.3 "Counting Leading 0's"
// floor(log_2(x))
// 0 => bits

  t_int res, i;

  if (x == 0) {
    res = bits;
    }
  else {
    res = bits-1;
    for (i = ld_bits-1; i >= 0; --i) {  // UNROLL
      if ((x & ~a_bfly_mask[i]) != 0) {
        x = x >> (1 << i);
        res = res - (1 << i);
        }
      }
    }
  return res;
  }

mycall t_int nr_trailing_0bits(t_bits x) {
// See Hacker's Delight, 5.4 "Counting Trailing 0's"
// 0 => bits

  t_int res, i;

  if (x == 0) {
    res = bits;
    }
  else {
    res = 0;
    for (i = ld_bits-1; i >= 0; --i) {  // UNROLL
      if ((x & a_sw_base[i]) == 0) {
        x = x >> (1 << i);
        res = res + (1 << i);
        }
      }
    }
  return res;
  }

mycall t_bool is_contiguous_1bits(t_bits x) {
// Is x a contiguous bit string?

  return ((((x - 1) | x) + 1) & x) == 0;
  }

mycall t_bits tbm(t_bits x, t_int mode) {
// General trailing bit modifications.
// 2014-02-11 by Jasper L. Neumann

  switch (mode & 0x1f) {
    case 0x00: return 0;
    case 0x01: return x & ~(x+1);
    case 0x02: return ~x & (x+1);
    case 0x03: return x ^ (x+1);
    case 0x04: return ~(x ^ (x+1));
    case 0x05: return x | ~(x+1);
    case 0x06: return ~x | (x+1);
    case 0x07: return all_bits;
    case 0x08: return x & (x+1);
    case 0x09: return x;
    case 0x0a: return x+1;
    case 0x0b: return x | (x+1);
    case 0x0c: return ~(x | (x+1));
    case 0x0d: return ~x-1;
    case 0x0e: return ~x;
    case 0x0f: return ~(x & (x+1));
    case 0x10: return 0;
    case 0x11: return ~x & (x-1);
    case 0x12: return x & -x;
    case 0x13: return x ^ (x-1);
    case 0x14: return x ^ -x;
    case 0x15: return ~x | (x-1);
    case 0x16: return x | -x;
    case 0x17: return all_bits;
    case 0x18: return x & (x-1);
    case 0x19: return x-1;
    case 0x1a: return x;
    case 0x1b: return x | (x-1);
    case 0x1c: return ~(x | (x-1));
    case 0x1d: return ~x;
    case 0x1e: return -x;
    case 0x1f: return ~(x & (x-1));
    default:   return 0;  // this can not happen
    }
  }

mycall t_bits blend(t_bits m, t_bits x, t_bits y) {
// The bit equivalent to m?x:y.
// O(1)

  return (m & x) | (~m & y);  // | can be replaced by ^ or +
  // return ((x | y) & m) ^ y;
  }

mycall t_bits simd_odd(t_bits x, t_subword sw) {
// Set mask denoting odd values, i.e. propagate low bit to the left.
// 2014-07-12 by Jasper L. Neumann
// O(1)

  t_bits m;
  t_bits lsb;
  t_int shift;

  shift = 1 << sw;
  m = a_bfly_lo[sw];
  lsb = x & m;
  lsb = ((lsb << (shift-1)) << 1) - lsb;  // transform lsb: 1 ==> -1
  // A shift by >= bits is undefined by the C/C++ standard.
  // Here we do a double shift to avoid shift by bits.
  return lsb;
  }

mycall t_bits bit_permute_step(t_bits x, t_bits m, t_uint shift) {
// INLINE
// Can be replaced by bit_permute_step_simple,
// if for the relevant bits n the following holds:
// nr_1bits(bit_permute_step_simple(n,m,shift)) == nr_1bits(n)
// x86: >= 6/5 cycles
// ARM: >= 4/4 cycles

  t_bits t;

  // assert(((m << shift) & m) == 0);
  // assert(((m << shift) >> shift) == m);
  t = ((x >> shift) ^ x) & m;
  x = x ^ t;  t = t << shift;  x = x ^ t;  // x = (x ^ t) ^ (t << shift);
  return x;
  }

mycall t_bits bit_permute_step_simple(t_bits x, t_bits m, t_uint shift) {
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

mycall void bit_permute_step2(t_bits* x1, t_bits* x2, t_bits m, t_uint shift) {
// INLINE
// Extended variant of bit_permute_step.
// Will be slow if not inlined.

  t_bits t;

  t = ((*x2 >> shift) ^ *x1) & m;
  *x1 = *x1 ^ t;
  *x2 = *x2 ^ (t << shift);
  }

mycall void identity_perm(ta_index tgt) {

  t_int i;

  for (i = 0; i <= bits-1; ++i) {
    tgt[i] = i;
    }
  }

mycall void invert_perm(const ta_index src, ta_index tgt) {

  t_int i;

  for (i = 0; i <= bits-1; ++i) {
    tgt[i] = no_index;
    }
  for (i = 0; i <= bits-1; ++i) {
    if (src[i] != no_index) {
      tgt[src[i]] = i;
      }
    }
  }

mycall void random_perm(ta_index tgt) {

  t_int j,q;
  t_int x;

  for (q = 0; q <= bits-1; ++q) {
    tgt[q] = q;
    }

  for (q = 1; q <= bits-1; ++q) {
    //j = random_int(q+1);  // 0..q
    x = tgt[q];
    tgt[q] = tgt[j];
    tgt[j] = x;
    }
  }

mycall t_bits used_source_bits(const ta_index perm) {
// 2012-09-14 by Jasper L. Neumann

  t_int i;
  t_bits n_src;

  n_src = 0;  // set of needed source bits
  for (i = 0; i <= bits-1; ++i) {
    if (perm[i] != no_index) {
      n_src = n_src | (lo_bit << perm[i]);
      }
    }
  return n_src;
  }

mycall t_bits used_target_bits(const ta_index perm) {
// 2012-09-14 by Jasper L. Neumann

  t_int i;
  t_bits n_tgt;
  t_bits mask;

  n_tgt = 0;  // set of needed target bits
  mask = 1;
  for (i = 0; i <= bits-1; ++i) {
    if (perm[i] != no_index) {
      n_tgt = n_tgt | mask;
      }
    mask = mask << 1;
    }
  return n_tgt;
  }


//////
// Bit index operations

mycall t_bits bit_index_complement(t_bits x, t_subword k) {
// INLINE
// See ARM System Developer's Guide, 7.6.2 "Bit Permutations"
// as used in the loop of general_reverse_bits.

  t_int shift;
  t_bits m;

  shift = 1 << k;
  m = a_bfly_mask[k];
  return bit_permute_step_simple(x,m,shift);
  }

mycall t_bits bit_index_swap(t_bits x, t_subword j, t_subword k) {
// INLINE
// See ARM System Developer's Guide, 7.6.2 "Bit Permutations"
// => shuffle, unshuffle

  t_int shift;
  t_subword q;
  t_bits m;

  if (j != k) {
    if (j < k) {
      q = j;
      j = k;
      k = q;
      }
    shift = (1 << (t_int)(j)) - (1 << (t_int)(k));
    m = a_bfly_mask[(t_int)(j)] & ~a_bfly_mask[(t_int)(k)];  // b_j==0 & b_k==1
    x = bit_permute_step(x,m,shift);
    }
  return x;
  }

mycall t_bits bit_index_swap_complement(t_bits x, t_subword j, t_subword k) {
// INLINE
// See ARM System Developer's Guide, 7.6.2 "Bit Permutations"

  t_int shift;
  t_subword q;
  t_bits m;

  if (j != k) {
    if (j < k) {
      q = j;
      j = k;
      k = q;
      }
    shift = (1 << j) + (1 << k);
    m = a_bfly_mask[j] & a_bfly_mask[k];  // b_j==0 & b_k==0
    x = bit_permute_step(x,m,shift);
    }
  return x;
  }

mycall t_bits bit_index_ror(t_bits x, t_subword ofs, t_subword field, t_int rot) {
// Rotate an bit index field to the right by rot.
// q+field+ofs=ld_bits
// q: upper bit indexes (unchanged)
// field: size of the affected bit string
// ofs: lower bit indexes (unchanged)
// Bit-parallel implementation: 2011-10-04 by Jasper L. Neumann

  t_int i,j,d,g,k,n;

  if (field > 0) {
    rot = rot % (t_int)(field);  // rot might be negative
    if (rot != 0) {
      if (rot < 0) {
        // we need a real modulo operation yielding 0..field-1
        rot = rot + field;
        }
      g = gcd(field,rot);
      d = field / g;
      for (i = 0; i <= g-1; ++i) {
        k = i;
        for (j = 0; j <= d-2; ++j) {
          n = k + rot;
          if (n >= (t_int)(field)) {
            // avoid mod
            n = n - field;
            }
          x = bit_index_swap(x,n+ofs,k+ofs);
          k = n;
          }
        }
      }
    }
  return x;
  }

mycall t_bits transpose(t_bits x, t_subword ld_fields, t_subword ld_col, t_subword ld_row) {
// Transpose bit matrixes.
// ld_fields: width of the bit fields
// ld_col: ld(bit columns)
// ld_row: ld(bit rows)
// 2011-10-04 by Jasper L. Neumann

  return bit_index_ror(x,ld_fields,ld_col+ld_row,ld_col);
  }

mycall t_bits shuffle_power(t_bits x, t_subword sw1, t_subword sw2, t_int pwr) {
// pwr times shuffle(x,sw1,sw2)
// See Hacker's Delight, 7.6/7.8 "Rearrangements and Index Transformations"
// 2011-10-04 by Jasper L. Neumann

  return bit_index_ror(x,sw1,sw2-sw1,-pwr);
  }

mycall t_bits unshuffle_power(t_bits x, t_subword sw1, t_subword sw2, t_int pwr) {
// pwr times unshuffle(x,sw1,sw2)
// See Hacker's Delight, 7.6/7.8 "Rearrangements and Index Transformations"
// 2011-10-04 by Jasper L. Neumann

  return bit_index_ror(x,sw1,sw2-sw1,pwr);
  }

mycall t_bits permute_bpc(t_bits x, const ta_subword tgt, t_subword_set k) {
// Do a BPC permutation via bit index operation primitives.
// tgt: permutation vector, must hold 0..ld_bits-1, all different
// k: complement set, see general_reverse_bits
// 2012-02-07 by Jasper L. Neumann

  ta_subword src,inv_src;
  t_subword i,j,ii;
  t_subword_set n,m,kk;

  for (i = 0; i <= ld_bits-1; ++i) {
    src[i] = i;
    inv_src[i] = i;
    }

  kk = 0;  // k for generated

  for (i = 0; i <= ld_bits-1; ++i) {  // any order
    n = 1 << i;
    ii = src[i];
    if (tgt[i] == ii) {
      if (((k ^ kk) & n) != 0) {
        x = bit_index_complement(x,i);
        }
      }
    else {
      j = inv_src[tgt[i]];
      m = 1 << j;
      if (((k & n) != 0) ^ ((kk & m) != 0)) {  // boolean ^
        x = bit_index_swap_complement(x,i,j);
        if ((kk & n) == 0) {
          kk = kk | m;
          }
        else {
          kk = kk & ~m;
          }
        }
      else {
        x = bit_index_swap(x,i,j);
        if ((kk & n) != 0) {
          kk = kk | m;
          }
        else {
          kk = kk & ~m;
          }
        }

      inv_src[ii] = j;
      src[j] = ii;
      }
    }
  return x;
  }

mycall void invert_bpc(const ta_subword src, t_subword_set src_k, ta_subword tgt, t_subword_set* tgt_k) {

  t_subword i;

  for (i = 0; i <= ld_bits-1; ++i) {
    tgt[src[i]] = i;
    }
  *tgt_k = 0;
  for (i = 0; i <= ld_bits-1; ++i) {
    if ((src_k & (1 << tgt[i])) != 0) {
      *tgt_k = *tgt_k | (1 << i);
      }
    }
  }


//////
// Generalized Bit Reversal

mycall t_bits general_reverse_bits(t_bits x, t_int k) {
// Swap all subwords of given levels.
// See Hacker's Delight, 7.1 "Generalized Bit Reversal"
// k: set of t_subword, i.e. one bit per subword size.

  t_int i,j;
  t_bits m;

  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    j = 1 << i;
    if ((k & j) != 0) {
      // x = bit_index_complement(x,j);
      m = a_bfly_mask[i];
      x = bit_permute_step_simple(x,m,j);
      }
    }
  return x;
  }

mycall t_bits bswap(t_bits x) {
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


//////
// Swap by primitives

mycall t_bits prim_swap(t_bits x, t_bits m) {
// Swap by primitives.
// Generalization of general_reverse_bits.
// See "Matters Computational" by Joerg Arndt, "A restricted method"
// See Hacker's Delight, 7.1 "Generalized Bit Reversal"
// Bit-parallel implementation: 2011-02 by Jasper L. Neumann

  t_int i,s;
  t_bits q;

  if ((m & hi_bit) != 0) {  // highest bit set?
    // normal operation
    for (i = ld_bits-1; i >= 0; --i) {  // UNROLL
      q = m & a_prim_swap[i];
      s = 1 << i;
      q = (q << 1) - (q >> (s-1));  // broadcast bits
      x = bit_permute_step(x,q,s);
      }
    }
  else {
    // inverse operation
    // same as above but with reversed loop
    for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
      q = m & a_prim_swap[i];
      s = 1 << i;
      q = (q << 1) - (q >> (s-1));  // broadcast bits
      x = bit_permute_step(x,q,s);
      }
    }
  return x;
  }


//////
// Shuffle and unshuffle

mycall t_bits shuffle(t_bits x, t_subword sw1, t_subword sw2) {
// Shuffle/zip/interlace entities.
// See Hacker's Delight, 7.2 "Shuffling Bits"
// sw1: log_2(subword_length): entities to move
// sw2: log_2(word_length): moving area
// 0 <= sw1 < sw2 <= ld_bits
// Example: sw1=2, sw2=5: Shuffle nibbles in dword
// See shuffle_power.

  t_int i,j;

  if (sw2 >= 2) {
    // same code as ibfly but with special non-conforming masks
    for (i = (t_int)(sw2)-2; i >= (t_int)(sw1); --i) {  // UNROLL?
      // x = bit_index_swap(x,i+1,i);
      j = 1 << i;
      x = bit_permute_step(x, a_shuffle_mask[i], j);
      }
    }
  return x;
  }

mycall t_bits unshuffle(t_bits x, t_subword sw1, t_subword sw2) {
// Unshuffle/unzip/uninterlace entities.
// See Hacker's Delight, 7.2 "Shuffling Bits"
// sw1: log_2(subword_length): entities to move
// sw2: log_2(word_length): moving area
// 0 <= sw1 < sw2 <= ld_bits
// Example: sw1=0, sw2=3: Unshuffle bits in bytes
// See unshuffle_power.

  t_int i,j;

  if (sw2 >= 2) {
    // same code as bfly but with special non-conforming masks
    for (i = (t_int)(sw1); i <= (t_int)(sw2)-2; ++i) {  // UNROLL?
      // x = bit_index_swap(x,i+1,i);
      j = 1 << i;
      x = bit_permute_step(x, a_shuffle_mask[i], j);
      }
    }
  return x;
  }


//////
// A "class" for butterfly and other operations


//////
// Compress and expand

//////
// Compress and expand: Compress bit masks

mycall t_bits compress_mask_right(t_bits m, t_subword sw) {
// Move all 1 bits of m to the right of their subwords.
// This is essentially the same code as gen_ce_right.
// We only need the final value of m.
// See Hacker's Delight, 7.4 "Compress, or Generalized Extract"
// UNFOLD all cases of sw
// =compress_right(m,m,sw)

  t_bits mk, mp, mv, mm, m0;
  t_int i,j,s;

  if (sw > 0) {
    m0 = a_bfly_lo[sw];
    mk = ((~m) << 1) & ~m0;             // We will count 0's to right
    for (i = 0; i <= (t_int)(sw)-1; ++i) {  // UNROLL
      mp = mk;                          // Parallel suffix
      for (j = 0; j <= (t_int)(sw)-1; ++j) {  // UNROLL
        s = 1 << j;
        mm = ~(m0 << s) + m0;
        mp = mp ^ ((mp << s) & mm);     // Masking not needed for sw=ld_bits
        }
      mv = mp & m;                      // Bits to move
      m = (m ^ mv) | (mv >> (1 << i));  // Compress m
      mk = mk & ~mp;
      }
    }
  return m;
  }

mycall t_bits compress_mask_left(t_bits m, t_subword sw) {
// Move all 1 bits of m to the left of their subwords.
// This is essentially the same code as gen_ce_left.
// We only need the final value of m.
// See Hacker's Delight, 7.4 "Compress, or Generalized Extract"
// UNFOLD all cases of sw
// =compress_left(m,m,sw)

  t_bits mk, mp, mv, mm, m0, m1;
  t_int i,j,s;

  if (sw > 0) {
    m1 = a_bfly_lo[sw];
    m0 = (m1 >> 1) + hi_bit;            // m0 = a_bfly_hi[sw];
    mk = ((~m) >> 1) & ~m0;             // We will count 0's to right
    for (i = 0; i <= (t_int)(sw)-1; ++i) {  // UNROLL
      mp = mk;                          // Parallel suffix
      for (j = 0; j <= (t_int)(sw)-1; ++j) {  // UNROLL
        s = 1 << j;
        mm = (m0 >> (s-1)) - m1;
        mp = mp ^ ((mp >> s) & mm);     // Masking not needed for sw=ld_bits
        }
      mv = mp & m;                      // Bits to move
      m = (m ^ mv) | (mv << (1 << i));  // Compress m
      mk = mk & ~mp;
      }
    }
  return m;
  }

mycall t_bits compress_mask(t_bits m, t_subword sw, t_direction d) {
// INLINE

  switch (d) {
    case right: return compress_mask_right(m,sw);
    case left:  return compress_mask_left(m,sw);
    default:    return 0;  // this can't happen
    }
  }


//////
// Compress and expand: Generate configuration

mycall void gen_ce_right(tr_bfly* self, t_bits m, t_subword sw) {
// See Hacker's Delight, 7.4 "Compress, or Generalized Extract"
// See Hacker's Delight, 7.5 "Expand, or Generalized Insert"
// To compress use apply_compress_right.
// To expand use apply_expand_right.
// UNFOLD all cases of sw

  t_bits mk, mp, mv, mm, m0;
  t_int i,j,s;

  self->mask = m;                       // Save original mask
  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    self->cfg[i] = 0;
    }
  if (sw > 0) {
    m0 = a_bfly_lo[sw];
    mk = ((~m) << 1) & ~m0;             // We will count 0's to right
    for (i = 0; i <= (t_int)(sw)-1; ++i) {  // UNROLL
      mp = mk;                          // Parallel suffix
      for (j = 0; j <= (t_int)(sw)-1; ++j) {  // UNROLL
        s = 1 << j;
        mm = ~(m0 << s)+m0;
        mp = mp ^ ((mp << s) & mm);     // Masking not needed for sw=ld_bits
        }
      mv = mp & m;                      // Bits to move
      self->cfg[i] = mv;
      m = (m ^ mv) | (mv >> (1 << i));  // Compress m
      mk = mk & ~mp;
      }
    }
  }

mycall void gen_ce_left(tr_bfly* self, t_bits m, t_subword sw) {
// See Hacker's Delight, 7.4 "Compress, or Generalized Extract"
// See Hacker's Delight, 7.5 "Expand, or Generalized Insert"
// To compress use apply_compress_left.
// To expand use apply_expand_left.
// UNFOLD all cases of sw

  t_bits mk, mp, mv, mm, m0, m1;
  t_int i,j,s;

  self->mask = m;                       // Save original mask
  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    self->cfg[i] = 0;
    }
  if (sw > 0) {
    m1 = a_bfly_lo[sw];
    m0 = (m1 >> 1) + hi_bit;            // m0 = a_bfly_hi[sw];
    mk = ((~m) >> 1) & ~m0;             // We will count 0's to right
    for (i = 0; i <= (t_int)(sw)-1; ++i) {  // UNROLL
      mp = mk;                          // Parallel suffix
      for (j = 0; j <= (t_int)(sw)-1; ++j) {  // UNROLL
        s = 1 << j;
        mm = (m0 >> (s-1)) - m1;
        mp = mp ^ ((mp >> s) & mm);     // Masking not needed for sw=ld_bits
        }
      mv = mp & m;                      // Bits to move
      self->cfg[i] = mv;
      m = (m ^ mv) | (mv << (1 << i));  // Compress m
      mk = mk & ~mp;
      }
    }
  }


//////
// Compress and expand: Usage

mycall t_bits apply_compress_right(const tr_bfly* self, t_bits x) {
// See Hacker's Delight, 7.4 "Compress, or Generalized Extract"
// self should be configured by gen_ce_right.

  t_bits t;
  t_int i;

  x = x & self->mask;  // Clear irrelevant bits

  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    t = x & self->cfg[i];
    x = (x ^ t) | (t >> (1 << i));  // Compress x (or ^)
    }

  return x;
  }

mycall t_bits apply_compress_left(const tr_bfly* self, t_bits x) {
// See Hacker's Delight, 7.4 "Compress, or Generalized Extract"
// self should be configured by gen_ce_left.

  t_bits t;
  t_int i;

  x = x & self->mask;  // Clear irrelevant bits

  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    t = x & self->cfg[i];
    x = (x ^ t) | (t << (1 << i));  // Compress x (or ^)
    }

  return x;
  }

mycall t_bits apply_expand_right(const tr_bfly* self, t_bits x) {
// See Hacker's Delight, 7.4 "Compress, or Generalized Extract"
// See Hacker's Delight, 7.5 "Expand, or Generalized Insert"
// self should be configured by gen_ce_right.
// (a & b) | (~a & c) => ((b ^ c) & a) ^ c

  t_int i;

  for (i = ld_bits-1; i >= 0; --i) {  // UNROLL
    x = (((x << (1 << i)) ^ x) & self->cfg[i]) ^ x;
    }

  return x & self->mask;  // Clear out extraneous bits
  }

mycall t_bits apply_expand_left(const tr_bfly* self, t_bits x) {
// See Hacker's Delight, 7.4 "Compress, or Generalized Extract"
// See Hacker's Delight, 7.5 "Expand, or Generalized Insert"
// self should be configured by gen_ce_left.
// (a & b) | (~a & c) => ((b ^ c) & a) ^ c

  t_int i;

  for (i = ld_bits-1; i >= 0; --i) {  // UNROLL
    x = (((x >> (1 << i)) ^ x) & self->cfg[i]) ^ x;
    }

  return x & self->mask;  // Clear out extraneous bits
  }


//////
// Compress and expand: Compound

mycall t_bits compress_right(t_bits x, t_bits m, t_subword sw) {

  tr_bfly ce;

  gen_ce_right(&ce,m,sw);
  return apply_compress_right(&ce,x);
  }

mycall t_bits compress_left(t_bits x, t_bits m, t_subword sw) {

  tr_bfly ce;

  gen_ce_left(&ce,m,sw);
  return apply_compress_left(&ce,x);
  }

mycall t_bits compress(t_bits x, t_bits m, t_subword sw, t_direction d) {
// INLINE

  switch (d) {
    case right: return compress_right(x,m,sw);
    case left:  return compress_left(x,m,sw);
    default:    return 0;  // this can't happen
    }
  }

mycall t_bits expand_right(t_bits x, t_bits m, t_subword sw) {

  tr_bfly ce;

  gen_ce_right(&ce,m,sw);
  return apply_expand_right(&ce,x);
  }

mycall t_bits expand_left(t_bits x, t_bits m, t_subword sw) {

  tr_bfly ce;

  gen_ce_left(&ce,m,sw);
  return apply_expand_left(&ce,x);
  }

mycall t_bits expand(t_bits x, t_bits m, t_subword sw, t_direction d) {
// INLINE

  switch (d) {
    case right: return expand_right(x,m,sw);
    case left:  return expand_left(x,m,sw);
    default:    return 0;  // this can't happen
    }
  }


//////
// Butterfly network

mycall t_bits butterfly(t_bits x, t_bits m, t_subword sw) {
// INLINE
// One butterfly step/stage.
// sw: 0..ld_bits-1
// m & a_bfly_mask[sw] should be == m

  return bit_permute_step(x, m, 1 << sw);
  }

mycall t_bits bfly(const tr_bfly* self, t_bits x) {
// Apply butterfly network on x configured by
//   - gen_frot
//   - gen_vrot
//   - gen_cef_right
//   - gen_cef_left

  t_int stage,j;

  for (stage = ld_bits-1; stage >= 0; --stage) {  // UNROLL
    // x = butterfly(x, self->cfg[stage], sw);
    j = 1 << stage;
    x = bit_permute_step(x, self->cfg[stage], j);
    }

  return x;
  }

mycall t_bits ibfly(const tr_bfly* self, t_bits x) {
// Apply inverse butterfly network on x configured by
//   - gen_frot
//   - gen_vrot
//   - gen_cef_right
//   - gen_cef_left

  t_int stage,j;

  for (stage = 0; stage <= ld_bits-1; ++stage) {  // UNROLL
    // x = butterfly(x, self->cfg[stage], sw);
    j = 1 << stage;
    x = bit_permute_step(x, self->cfg[stage], j);
    }

  return x;
  }

mycall t_bool bfly_parity(const tr_bfly* self) {
// Return the parity of a permutation given by a butterfly network.
// This is false for even parity and true for odd parity.

  t_int stage;
  t_bits x;

  x = 0;
  for (stage = 0; stage <= ld_bits-1; ++stage) {  // UNROLL
    x = x ^ self->cfg[stage];
    }

  return odd(nr_1bits(x));
  }


//////
// Rotate via butterfly

//////
// Rotate via butterfly: Generate configuration

mycall void gen_frot(tr_bfly* self, t_int rot, t_subword sw) {
// Fixed rotate.
// Generate configuration for [inverse] butterfly network.
// To rotate right use bfly.
// To rotate left use ibfly.
// Bit-parallel implementation: 2011-09-14 by Jasper L. Neumann

  t_int i;

  self->mask = rot;
  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    self->cfg[i] = 0;
    }
  if (sw > 0) {
    for (i = 0; i <= (t_int)(sw)-1; ++i) {  // UNROLL
      self->cfg[i] = rolc_lo(0,rot,i) * a_bfly_lo[i+1];
      }
    }
  }

mycall void gen_vrot(tr_bfly* self, t_bits rot, t_subword sw) {
// Field variable rotate.
// Generate configuration for [inverse] butterfly network.
// To rotate right use bfly.
// To rotate left use ibfly.
// Simulate rolc for every subword.
// Bit-parallel implementation: 2011-09-14 by Jasper L. Neumann

  t_bits t,x,y,lo;
  t_int i,s;

  self->mask = rot;
  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    self->cfg[i] = 0;
    }
  switch (sw) {
    case 0: {
      break;  // nothing else to do
      }
    case 1: {
      self->cfg[0] = rot & a_bfly_mask[0];
      break;
      }
    default: {  // UNFOLD all cases of sw
      // this code does not work for sw < 1
      // shift a single 1 to the left...
      lo = a_bfly_lo[sw];
      sw = sw-1;
      x = lo;
      for (i = 0; i <= (t_int)(sw)-1; ++i) {  // UNROLL
        y = (rot >> i) & lo;  // rot bit to #0
        s = 1 << i;
        // (lo_bit << s) - 1 == a_sw_base[i]
        // t = x & (y * ((lo_bit << s) - 1));  // get offending bit(s)
        t = x & ((y << s) - y);  // get offending bit(s)
        x = (x ^ t) ^ (t << s);  // swap subwords if bit set
        }
      // x is e.g. 1000 here (1 << 3), we want 3 ones
      x = x - lo;  // sub 1 to yield 1-string, e.g. 1000-1=111
      y = (rot >> sw) & lo;
      s = 1 << sw;
      // (lo_bit << s) - 1 == a_sw_base[sw]
      // x = x ^ (y * ((lo_bit << s) - 1));  // invert if rot<0
      x = x ^ ((y << s) - y);  // invert if rot<0
      x = x & a_bfly_mask[sw];  // finalize rolc
      self->cfg[sw] = x;
      // and now for the lower stages...
      for (i = (t_int)(sw)-1; i >= 0; --i) {  // sw-1..0; UNROLL
        // xor 2 columns together to get new rolc for the stage...
        s = 1 << i;
        x = (x ^ (x >> s)) & a_bfly_mask[i];
        x = x | (x << (s * 2));  // ...and spread into places
        self->cfg[i] = x;
        }
      break;
      }
    }
  }


//////
// Rotate via butterfly: Compound

mycall t_bits fror_bfly(t_bits x, t_int rot, t_subword sw) {
// Rotate right all subwords by the same amount.
// 2011-09-14 by Jasper L. Neumann

  tr_bfly cef;

  gen_frot(&cef,rot,sw);
  return bfly(&cef,x);
  }

mycall t_bits frol_bfly(t_bits x, t_int rot, t_subword sw) {
// Rotate left all subwords by the same amount.
// 2011-09-14 by Jasper L. Neumann

  tr_bfly cef;

  gen_frot(&cef,rot,sw);
  return ibfly(&cef,x);
  }

mycall t_bits frot_bfly(t_bits x, t_int rot, t_subword sw, t_direction d) {
// Rotate all subwords by the same amount.
// INLINE
// 2011-09-14 by Jasper L. Neumann

  switch (d) {
    case right: return fror_bfly(x,rot,sw);
    case left:  return frol_bfly(x,rot,sw);
    default:    return 0;  // this can't happen
    }
  }

mycall t_bits vror_bfly(t_bits x, t_bits rot, t_subword sw) {
// Rotate right all subwords by a variable amount.
// 2011-09-14 by Jasper L. Neumann

  tr_bfly cef;

  gen_vrot(&cef,rot,sw);
  return bfly(&cef,x);
  }

mycall t_bits vrol_bfly(t_bits x, t_bits rot, t_subword sw) {
// Rotate left all subwords by a variable amount.
// 2011-09-14 by Jasper L. Neumann

  tr_bfly cef;

  gen_vrot(&cef,rot,sw);
  return ibfly(&cef,x);
  }

mycall t_bits vrot_bfly(t_bits x, t_bits rot, t_subword sw, t_direction d) {
// Rotate all subwords by a variable amount.
// INLINE
// 2011-09-14 by Jasper L. Neumann

  switch (d) {
    case right: return vror_bfly(x,rot,sw);
    case left:  return vrol_bfly(x,rot,sw);
    default:    return 0;  // this can't happen
    }
  }

mycall t_bits frol(t_bits x, t_int rot, t_subword sw) {
// Fast alternative to frol_bfly.
// Every (1 << sw) bits are rotated left.
// x: value
// sw: log_2(#bits), must be <=ld_bits
// rot: rotate count
// Bit-parallel implementation: 2011-09-21 by Jasper L. Neumann
// = fror(x,-rot,sw)

  t_int b;  // # affected bits
  t_int r;  // rot % b
  t_bits m;  // mask for affected bits

  b = 1 << sw;

  r = rot & (b - 1);
  if (r == 0) {
    // Prevent shifting by b-r >= bits.
    return x;
    }
  else {
    m = a_bfly_lo[sw];
    m = (m << r) - m;

    return
      ((x << r) & ~m) |
      ((x >> (b - r)) & m);
    }
  }

mycall t_bits fror(t_bits x, t_int rot, t_subword sw) {
// Fast alternative to fror_bfly.
// Every (1 << sw) bits are rotated right.
// x: value
// sw: log_2(#bits), must be <=ld_bits
// rot: rotate count
// Bit-parallel implementation: 2011-09-21 by Jasper L. Neumann
// = frol(x,-rot,sw)

  t_int b;  // # affected bits
  t_int r;  // rot % b
  t_bits m;  // mask for affected bits

  b = 1 << sw;

  r = (b - rot) & (b - 1);
  if (r == 0) {
    // Prevent shifting by b-r >= bits.
    return x;
    }
  else {
    m = a_bfly_lo[sw];
    m = (m << r) - m;

    return
      ((x << r) & ~m) |
      ((x >> (b - r)) & m);
    }
  }

mycall t_bits frot(t_bits x, t_int rot, t_subword sw, t_direction d) {
// Fast alternative to frot_bfly.
// INLINE
// 2011-09-21 by Jasper L. Neumann

  switch (d) {
    case right: return fror(x,rot,sw);
    case left:  return frol(x,rot,sw);
    default:    return 0;  // this can't happen
    }
  }

mycall t_bits frolc(t_bits x, t_int rot, t_subword sw) {
// Every (1 << sw) bits are rotated left.
// x: value
// sw: log_2(#bits), must be <=ld_bits
// rot: rotate count
// Bit-parallel implementation: 2011-09-23 by Jasper L. Neumann

  t_int b;   // # affected bits
  t_int r;   // rot mod b
  t_bits m;  // mask for affected bits

  b = 1 << sw;

  r = rot & (b-1);
  if (r == 0) {
    // Prevent shifting by b-r >= bits.
    }
  else {
    m = a_bfly_lo[sw];
    m = (m << r) - m;

    x =
      ((x << r) & ~m) |
      ((x >> (b-r)) & m);

    // Until here essentially same code as frol.
    x = x ^ m;
    }

  if ((rot & b)!=0) {
    x = ~x;
    }

  return x;
  }

mycall t_bits frorc(t_bits x, t_int rot, t_subword sw) {
// Every (1 << sw) bits are rotated right.
// x: value
// sw: log_2(#bits), must be <=ld_bits
// rot: rotate count
// INLINE

  return frolc(x, -rot, sw);
  }

mycall t_bits frotc(t_bits x, t_int rot, t_subword sw, t_direction d) {
// INLINE

  switch (d) {
    case right: return frorc(x,rot,sw);
    case left:  return frolc(x,rot,sw);
    default:    return 0;  // this can't happen
    }
  }

mycall t_bits vrol(t_bits x, t_bits rot, t_subword sw) {
// Field variable variant of frol.
// Fast alternative to vrol_bfly.
// Gives correct results for all values of rot.
// 2014-07-14 by Jasper L. Neumann
// O(ld_bits)

  t_int i, s;

  s = 1;
  for (i = 0; i <= (t_int)(sw)-1; ++i) {
    // s = 1 << i;
    x = blend(simd_odd(rot, sw), frol(x, s, sw), x);
    rot = rot >> 1;
    s = s << 1;
    }
  return x;
  }

mycall t_bits vror(t_bits x, t_bits rot, t_subword sw) {
// Field variable variant of fror.
// Fast alternative to vror_bfly.
// Gives correct results for all values of rot.
// 2014-07-14 by Jasper L. Neumann
// O(ld_bits)

  t_int i, s;

  s = 1;
  for (i = 0; i <= (t_int)(sw)-1; ++i) {
    // s = 1 << i;
    x = blend(simd_odd(rot, sw), fror(x, s, sw), x);
    rot = rot >> 1;
    s = s << 1;
    }
  return x;
  }

mycall t_bits vrot(t_bits x, t_bits rot, t_subword sw, t_direction d) {
// Fast alternative to vrot_bfly.
// INLINE
// 2011-09-21 by Jasper L. Neumann

  switch (d) {
    case right: return vror(x,rot,sw);
    case left:  return vrol(x,rot,sw);
    default:    return 0;  // this can't happen
    }
  }


//////
// Compress/expand-flip via butterfly

//////
// Compress/expand-flip via butterfly: Generate configuration

mycall void gen_cef_right(tr_bfly* self, t_bits m, t_subword sw) {
// Scatter/gather-flip, compress/expand+flip.
// Generate configuration for [inverse] butterfly network.
// To compress use ibfly.
// To expand use bfly.
// Bit-parallel implementation: 2011-02 by Jasper L. Neumann

  t_bits t,mm,m0;
  t_int i,j,s;

  self->mask = m;
  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    self->cfg[i] = 0;
    }
  if (sw > 0) {  // UNFOLD all cases of sw
    m = ~m;
    m0 = a_bfly_lo[sw];

    for (i = 0; i <= (t_int)(sw)-1; ++i) {  // UNROLL
      t = m;
      for (j = i; j <= (t_int)(sw)-1; ++j) {  // UNROLL
        s = 1 << j;  // j ones; j=2: 1 + ~(1 << 4): 11101111+1=11110000
        mm = ~(m0 << s) + m0;  // mask to hinder shifting into other subwords
        m = m ^ ((m << s) & mm);
        }
      s = 1 << i;
      m = m & a_bfly_mask[i];  // my bfly looks on low bits
      self->cfg[i] = m;
      m = (t ^ (t >> s)) & m;  // do a butterfly op
      m = (t ^ m) ^ (m << s);
      }
    }
  }

mycall void gen_cef_left(tr_bfly* self, t_bits m, t_subword sw) {
// Scatter/gather-flip, compress/expand+flip.
// Generate configuration for [inverse] butterfly network.
// To compress use ibfly.
// To expand use bfly.
// Bit-parallel implementation: 2011-02 by Jasper L. Neumann

  t_bits t,mm,m0,m1;
  t_int i,j,s;

  self->mask = m;
  for (i = 0; i <= ld_bits-1; ++i) {  // UNROLL
    self->cfg[i] = 0;
    }
  if (sw > 0) {  // UNFOLD all cases of sw
    m = ~m;
    m1 = a_bfly_lo[sw];
    m0 = (m1 >> 1) + hi_bit;  // m0 = a_bfly_hi[sw];
    for (i = 0; i <= (t_int)(sw)-1; ++i) {  // UNROLL
      t = m;
      for (j = i; j <= (t_int)(sw)-1; ++j) {  // UNROLL
        s = 1 << j;  // j ones; j=2: 1 + ~(1 << 4): 11101111+1=11110000
        mm = (m0 >> (s-1)) - m1;  // mask to hinder shifting into other subwords
        m = m ^ ((m >> s) & mm);
        }
      s = 1 << i;
      m = (m >> s) & a_bfly_mask[i];  // my bfly looks on low bits
      self->cfg[i] = m;               // so shift into place
      m = (t ^ (t >> s)) & m;  // do a butterfly op
      m = (t ^ m) ^ (m << s);
      }
    }
  }


//////
// Compress/expand-flip via butterfly: Compound

mycall t_bits compress_flip_right(t_bits x, t_bits m, t_subword sw) {
// 2011-02 by Jasper L. Neumann

  tr_bfly cef;

  gen_cef_right(&cef,m,sw);
  return ibfly(&cef,x);
  }

mycall t_bits compress_flip_left(t_bits x, t_bits m, t_subword sw) {
// 2011-02 by Jasper L. Neumann

  tr_bfly cef;

  gen_cef_left(&cef,m,sw);
  return ibfly(&cef,x);
  }

mycall t_bits compress_flip(t_bits x, t_bits m, t_subword sw, t_direction d) {
// INLINE
// 2011-02 by Jasper L. Neumann

  switch (d) {
    case right: return compress_flip_right(x,m,sw);
    case left:  return compress_flip_left(x,m,sw);
    default:    return 0;  // this can't happen
    }
  }

mycall t_bits expand_flip_right(t_bits x, t_bits m, t_subword sw) {
// 2011-02 by Jasper L. Neumann

  tr_bfly cef;

  gen_cef_right(&cef,m,sw);
  return bfly(&cef,x);
  }

mycall t_bits expand_flip_left(t_bits x, t_bits m, t_subword sw) {
// 2011-02 by Jasper L. Neumann

  tr_bfly cef;

  gen_cef_left(&cef,m,sw);
  return bfly(&cef,x);
  }

mycall t_bits expand_flip(t_bits x, t_bits m, t_subword sw, t_direction d) {
// INLINE
// 2011-02 by Jasper L. Neumann

  switch (d) {
    case right: return expand_flip_right(x,m,sw);
    case left:  return expand_flip_left(x,m,sw);
    default:    return 0;  // this can't happen
    }
  }


//////
// Omega/flip

mycall t_bits omega(t_bits x, t_bits m, t_subword sw) {
// Simulate one omega step/stage.
// sw: 0..ld_bits-1
// m: low bits used
// m & a_bfly_mask[sw-1] should be == m
// 2011-07-18 by Jasper L. Neumann

  x = butterfly(x,m,sw-1);
  x = shuffle(x,0,sw);

  return x;
  }

mycall t_bits flip(t_bits x, t_bits m, t_subword sw) {
// Simulate one flip step/stage.
// sw: 0..ld_bits-1
// m: low bits used
// m & a_bfly_mask[sw-1] should be == m
// 2011-07-18 by Jasper L. Neumann

  x = unshuffle(x,0,sw);
  x = butterfly(x,m,sw-1);

  return x;
  }


//////
// Permutations via Benes network

static mycall void exchange_bit_index(t_bit_index* a, t_bit_index* b) {
// INLINE

  t_bit_index q;

  q = *a;
  *a = *b;
  *b = q;
  }

mycall void gen_benes_ex(tr_benes* self, const ta_index c_tgt, const ta_subword a_stage) {
// Generate a configuration for the Benes network with variable stage order.
// Use benes_fwd_ex and benes_bwd_ex.
// Algorithm as sketched by Donal E. Knuth,
//   The art of computer programming, vol. 4, pre-fascicle 1a
// Implemented 2011-12-02 by Jasper L. Neumann
// Modified 2012-08-31 to allow for "don't care" entries.

  ta_index src, inv_src;
  ta_index tgt, inv_tgt;
  t_int stage;
  t_int mask;
  t_bits cfg_src,cfg_tgt;
  t_bits src_set;
  t_int main_idx, aux_idx, src_idx, tgt_idx, idx2;
  t_int stage_idx;
  t_int s;

  for (s = 0; s <= bits-1; ++s) {
    src[s] = no_index;
    tgt[s] = no_index;
    }
  for (s = 0; s <= bits-1; ++s) {
    if (c_tgt[s] != no_index) {
      tgt[s] = s;
      src[c_tgt[s]] = s;
      }
    }
  invert_perm(src,inv_src);
  invert_perm(tgt,inv_tgt);
  for (stage_idx = 0; stage_idx <= ld_bits-1; ++stage_idx) {
    stage = a_stage[stage_idx];
    src_set = 0;
    mask = lo_bit << stage;
    cfg_src = 0;
    cfg_tgt = 0;
    for (main_idx = 0; main_idx <= bits-1; ++main_idx) {  // This order to meet Waksman test
      if ((main_idx & mask) == 0) {  // low only
        for (aux_idx = 0; aux_idx <= 1; ++aux_idx) {
          src_idx = main_idx+(aux_idx << stage);
          if (((lo_bit << src_idx) & src_set) == 0) {  // yet unhandled
            if (src[src_idx] != no_index) {  // not open

              do {
                src_set = src_set | (lo_bit << src_idx);
                tgt_idx = inv_tgt[src[src_idx]];
                if (tgt[tgt_idx] == no_index) {
                  break;  // open end
                  }

                if (((src_idx ^ tgt_idx) & mask) == 0) {
                  // straight
                  tgt_idx = tgt_idx ^ mask;
                  }
                else {
                  // cross
                  cfg_tgt = cfg_tgt | (lo_bit << (tgt_idx & ~mask));
                  idx2 = tgt_idx ^ mask;
                  exchange_bit_index(&tgt[tgt_idx], &tgt[idx2]);
                  inv_tgt[tgt[idx2]] = idx2;
                  if (tgt[tgt_idx] != no_index) {
                    inv_tgt[tgt[tgt_idx]] = tgt_idx;
                    }
                  }

                if (tgt[tgt_idx] == no_index) {
                  break;  // open end
                  }
                src_idx = inv_src[tgt[tgt_idx]];

                if (((src_idx ^ tgt_idx) & mask) == 0) {
                  // straight
                  src_set = src_set | (lo_bit << src_idx);
                  src_idx = src_idx ^ mask;
                  }
                else {
                  // cross
                  cfg_src = cfg_src | (lo_bit << (src_idx & ~mask));
                  idx2 = src_idx ^ mask;
                  src_set = src_set | (lo_bit << idx2);
                  exchange_bit_index(&src[src_idx], &src[idx2]);
                  inv_src[src[idx2]] = idx2;
                  if (src[src_idx] != no_index) {
                    inv_src[src[src_idx]] = src_idx;
                    }
                  }
                if (src[src_idx] == no_index) {
                  break;  // open end
                  }
                if (((lo_bit << src_idx) & src_set) != 0) {  // yet unhandled
                  break;  // already handled
                  }
                } while (true);

              }
            }
          }
        }
      }
    self->b1.cfg[stage] = cfg_src;
    self->b2.cfg[stage] = cfg_tgt;
    }
  // Reduce inner stages to one (not needed)
  // self->b2.cfg[0] = self->b2.cfg[0] ^ self->b1.cfg[0];
  // self->b1.cfg[0] = 0;
  }

mycall void gen_benes(tr_benes* self, const ta_index c_tgt) {
// INLINE
// Generate a configuration for the standard Benes network.

  gen_benes_ex(self,c_tgt,a_stage_bwd);  // standard Benes order
  }

mycall t_bits benes_fwd(const tr_benes* self, t_bits x) {
// Apply Benes network.
// c_tgt of gen_benes selected source indexes.

  return ibfly(&self->b2, bfly(&self->b1,x));
  }

mycall t_bits benes_bwd(const tr_benes* self, t_bits x) {
// Apply Benes network.
// c_tgt of gen_benes selected target indexes.

  return ibfly(&self->b1, bfly(&self->b2,x));
  }

mycall t_bits benes_fwd_ex(const tr_benes* self, t_bits x, const ta_subword a_stage) {
// Apply reordered Benes network.
// c_tgt of gen_benes_ex selected source indexes.

  t_int stage_idx;
  t_int stage;
  t_int j;

  //  benes_fwd = ibfly(self->b2, bfly(self->b1,x));
  for (stage_idx = 0; stage_idx <= ld_bits-1; ++stage_idx) {
    stage = a_stage[stage_idx];
    // x = butterfly(x, self->b1.cfg[stage], sw);
    j = 1 << stage;
    x = bit_permute_step(x, self->b1.cfg[stage], j);
    }
  for (stage_idx = ld_bits-1; stage_idx >= 0; --stage_idx) {
    stage = a_stage[stage_idx];
    // x = butterfly(x, self->b2.cfg[stage], sw);
    j = 1 << stage;
    x = bit_permute_step(x, self->b2.cfg[stage], j);
    }

  return x;
  }

mycall t_bits benes_bwd_ex(const tr_benes* self, t_bits x, const ta_subword a_stage) {
// Apply reordered Benes network.
// c_tgt of gen_benes_ex selected target indexes.

  t_int stage_idx;
  t_int stage;
  t_int j;

  //  benes_fwd = ibfly(self->b1, bfly(self->b2,x));
  for (stage_idx = 0; stage_idx <= ld_bits-1; ++stage_idx) {
    stage = a_stage[stage_idx];
    // x = butterfly(x, self->b2.cfg[stage], sw);
    j = 1 << stage;
    x = bit_permute_step(x, self->b2.cfg[stage], j);
    }
  for (stage_idx = ld_bits-1; stage_idx >= 0; --stage_idx) {
    stage = a_stage[stage_idx];
    // x = butterfly(x, self->b1.cfg[stage], sw);
    j = 1 << stage;
    x = bit_permute_step(x, self->b1.cfg[stage], j);
    }

  return x;
  }

mycall t_bool benes_parity(const tr_benes* self) {
// Return the parity of a permutation given by a Benes network.
// This is false for even parity and true for odd parity.

  t_int stage;
  t_bits x;

  x = 0;
  for (stage = 0; stage <= ld_bits-1; ++stage) {  // UNROLL
    x = x ^ self->b1.cfg[stage] ^ self->b2.cfg[stage];
    }

  return odd(nr_1bits(x));
  }

// eof.
