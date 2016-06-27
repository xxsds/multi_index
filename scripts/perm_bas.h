#pragma once

#include "general.h"

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


//////
// Random replacement

t_bits my_randseed;

mycall t_bits random_bits();


//////
// Auxiliary stuff

typedef enum {right,left} t_direction;

mycall t_bool odd(t_int x);
mycall t_longint gcd(t_longint a, t_longint b);
mycall t_bits mul_inv(t_bits x);
mycall t_bits rol(t_bits x, t_int rot);
mycall t_bits rol_lo(t_bits x, t_int rot, t_subword sw);
mycall t_bits rolc_lo(t_bits x, t_int rot, t_subword sw);
mycall t_bits gray_code(t_bits x);
mycall t_bits inv_gray_code(t_bits x);
mycall t_int nr_1bits(t_bits x);
mycall t_int nr_leading_0bits(t_bits x);
mycall t_int nr_trailing_0bits(t_bits x);
mycall t_bool is_contiguous_1bits(t_bits x);
mycall t_bits tbm(t_bits x, t_int mode);

mycall t_bits blend(t_bits m, t_bits x, t_bits y);
mycall t_bits simd_odd(t_bits x, t_subword sw);

mycall t_bits bit_permute_step(t_bits x, t_bits m, t_uint shift);
mycall t_bits bit_permute_step_simple(t_bits x, t_bits m, t_uint shift);
mycall void bit_permute_step2(t_bits* x1, t_bits* x2, t_bits m, t_uint shift);

mycall void identity_perm(ta_index tgt);
mycall void invert_perm(const ta_index src, ta_index tgt);
mycall void random_perm(ta_index tgt);
mycall t_bits used_source_bits(const ta_index perm);
mycall t_bits used_target_bits(const ta_index perm);


//////
// Bit index operations

mycall t_bits bit_index_complement(t_bits x, t_subword k);
mycall t_bits bit_index_swap(t_bits x, t_subword j, t_subword k);
mycall t_bits bit_index_swap_complement(t_bits x, t_subword j, t_subword k);

mycall t_bits bit_index_ror(t_bits x, t_subword ofs, t_subword field, t_int rot);
mycall t_bits transpose(t_bits x, t_subword ld_fields, t_subword ld_col, t_subword ld_row);
mycall t_bits shuffle_power(t_bits x, t_subword sw1, t_subword sw2, t_int pwr);
mycall t_bits unshuffle_power(t_bits x, t_subword sw1, t_subword sw2, t_int pwr);

mycall t_bits permute_bpc(t_bits x, const ta_subword tgt, t_subword_set k);
mycall void invert_bpc(const ta_subword src, t_subword_set src_k, ta_subword tgt, t_subword_set* tgt_k);


//////
// Generalized Bit Reversal

mycall t_bits general_reverse_bits(t_bits x, t_int k);
mycall t_bits bswap(t_bits x);


//////
// Swap by primitives

mycall t_bits prim_swap(t_bits x, t_bits m);


//////
// Shuffle and unshuffle

mycall t_bits shuffle(t_bits x, t_subword sw1, t_subword sw2);
mycall t_bits unshuffle(t_bits x, t_subword sw1, t_subword sw2);


//////
// A "class" for butterfly and other operations

typedef struct {
  // This structure is used to hold the configuration of
  // butterfly-based operations as well as compress and expand.

  t_bits cfg[ld_bits];  // butterfly configuration
  t_bits mask;  // saved mask, for compress/expand

  // Here is sketched how to convert this to a class:
  // Include all the generator and usage functions as private methods
  // and replace the parameter self by the implicit object pointer this.
  // Add the many compound routines.
  // Remove the name suffix  for all methods.
  // If you want to cache the configuration, add here:
  //   kind: the generator kind
  //     enum (initialized, frot, vrot, ce_right, ce_left, cef_right, cef_left)
  //   sw: the used subword size (t_subword)
  // Add an initializer/constructor which sets kind to initialized.
  // The generator routines must set the keys (kind, mask, sw).
  // The compound routines check the cached keys (kind, mask, sw);
  //   if not equal, call the generator routine and update the configuration;
  // finally they call the usage routine.
  } tr_bfly;


//////
// Compress and expand

//////
// Compress and expand: Compress bit masks

mycall t_bits compress_mask_right(t_bits m, t_subword sw);
mycall t_bits compress_mask_left(t_bits m, t_subword sw);
mycall t_bits compress_mask(t_bits m, t_subword sw, t_direction d);


//////
// Compress and expand: Generate configuration

mycall void gen_ce_right(tr_bfly* self, t_bits m, t_subword sw);
mycall void gen_ce_left(tr_bfly* self, t_bits m, t_subword sw);


//////
// Compress and expand: Usage

mycall t_bits apply_compress_right(const tr_bfly* self, t_bits x);
mycall t_bits apply_compress_left(const tr_bfly* self, t_bits x);
mycall t_bits apply_expand_right(const tr_bfly* self, t_bits x);
mycall t_bits apply_expand_left(const tr_bfly* self, t_bits x);


//////
// Compress and expand: Compound

mycall t_bits compress_right(t_bits x, t_bits m, t_subword sw);
mycall t_bits compress_left(t_bits x, t_bits m, t_subword sw);
mycall t_bits compress(t_bits x, t_bits m, t_subword sw, t_direction d);
mycall t_bits expand_right(t_bits x, t_bits m, t_subword sw);
mycall t_bits expand_left(t_bits x, t_bits m, t_subword sw);
mycall t_bits expand(t_bits x, t_bits m, t_subword sw, t_direction d);


//////
// Butterfly network

mycall t_bits butterfly(t_bits x, t_bits m, t_subword sw);
mycall t_bits bfly(const tr_bfly* self, t_bits x);
mycall t_bits ibfly(const tr_bfly* self, t_bits x);
mycall t_bool bfly_parity(const tr_bfly* self);


//////
// Rotate via butterfly

//////
// Rotate via butterfly: Generate configuration

mycall void gen_frot(tr_bfly* self, t_int rot, t_subword sw);
mycall void gen_vrot(tr_bfly* self, t_bits rot, t_subword sw);


//////
// Rotate via butterfly: Compound

mycall t_bits fror_bfly(t_bits x, t_int rot, t_subword sw);
mycall t_bits frol_bfly(t_bits x, t_int rot, t_subword sw);
mycall t_bits frot_bfly(t_bits x, t_int rot, t_subword sw, t_direction d);

mycall t_bits vror_bfly(t_bits x, t_bits rot, t_subword sw);
mycall t_bits vrol_bfly(t_bits x, t_bits rot, t_subword sw);
mycall t_bits vrot_bfly(t_bits x, t_bits rot, t_subword sw, t_direction d);

mycall t_bits frol(t_bits x, t_int rot, t_subword sw);
mycall t_bits fror(t_bits x, t_int rot, t_subword sw);
mycall t_bits frot(t_bits x, t_int rot, t_subword sw, t_direction d);

mycall t_bits frolc(t_bits x, t_int rot, t_subword sw);
mycall t_bits frorc(t_bits x, t_int rot, t_subword sw);
mycall t_bits frotc(t_bits x, t_int rot, t_subword sw, t_direction d);

mycall t_bits vrol(t_bits x, t_bits rot, t_subword sw);
mycall t_bits vror(t_bits x, t_bits rot, t_subword sw);
mycall t_bits vrot(t_bits x, t_bits rot, t_subword sw, t_direction d);


//////
// Compress/expand-flip via butterfly

//////
// Compress/expand-flip via butterfly: Generate configuration

mycall void gen_cef_right(tr_bfly* self, t_bits m, t_subword sw);
mycall void gen_cef_left(tr_bfly* self, t_bits m, t_subword sw);


//////
// Compress/expand-flip via butterfly: Compound

mycall t_bits compress_flip_right(t_bits x, t_bits m, t_subword sw);
mycall t_bits compress_flip_left(t_bits x, t_bits m, t_subword sw);
mycall t_bits compress_flip(t_bits x, t_bits m, t_subword sw, t_direction d);
mycall t_bits expand_flip_right(t_bits x, t_bits m, t_subword sw);
mycall t_bits expand_flip_left(t_bits x, t_bits m, t_subword sw);
mycall t_bits expand_flip(t_bits x, t_bits m, t_subword sw, t_direction d);


//////
// Omega/flip

mycall t_bits omega(t_bits x, t_bits m, t_subword sw);
mycall t_bits flip(t_bits x, t_bits m, t_subword sw);


//////
// Permutations via Benes network

typedef struct {
  tr_bfly b1,b2;
  } tr_benes;

mycall void gen_benes_ex(tr_benes* self, const ta_index c_tgt, const ta_subword a_stage);
mycall void gen_benes(tr_benes* self, const ta_index c_tgt);
mycall t_bits benes_fwd(const tr_benes* self, t_bits x);
mycall t_bits benes_bwd(const tr_benes* self, t_bits x);
mycall t_bits benes_fwd_ex(const tr_benes* self, t_bits x, const ta_subword a_stage);
mycall t_bits benes_bwd_ex(const tr_benes* self, t_bits x, const ta_subword a_stage);
mycall t_bool benes_parity(const tr_benes* self);

// eof.
