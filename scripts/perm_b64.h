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
// Last change: 2012-09-19

// Here the adaptations are made to enable perm_bas.c
// to be compiled for a word size of 64 bit.
// perm_bas.c must be included afterwards.


//////
// Our base for the bit hacks

#define ld_bits 6
  // log_2 of used bit size (here: 64 bit)
#define ld_bits_factorial (1*2*3*4*5*6)

typedef uint64_t t_bits;  // 1<<ld_bits bits, unsigned

#include "perm_bxx.h"


//////
// Derived stuff

const ta_subword a_stage_fwd = {0,1,2,3,4,5};
const ta_subword a_stage_bwd = {5,4,3,2,1,0};


//////
// Constant masks; must be adapted for other word sizes

const t_bits a_bfly_mask[]={
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

const t_bits a_bfly_lo[]={
  // 0..ld_bits
  // For auxiliary butterfly ops
  // a_bfly_mask with only lowest bit of runs set, index off by 1
  // = all_bits / ((1 << (1 << i)) - 1)
  0xffffffffffffffff,   // 0
  0x5555555555555555,   // 1 => a_bfly_mask[0]
  0x1111111111111111,   // 2
  0x0101010101010101,   // 3
  0x0001000100010001,   // 4
  0x0000000100000001,   // 5
  0x0000000000000001};  // 6

const t_bits a_bfly_hi[]={
  // Inverted a_bfly_mask with only highest bit of runs set, index off by 1.
  // = (a_bfly_lo[] >> 1)+hi_bit
  // = a_bfly_lo[] << ((1 << sw)-1)
  // = a_bfly_lo[] ror 1
  0xffffffffffffffff,   // 0
  0xaaaaaaaaaaaaaaaa,   // 1 => ~a_bfly_mask[0]
  0x8888888888888888,   // 2
  0x8080808080808080,   // 3
  0x8000800080008000,   // 4
  0x8000000080000000,   // 5
  0x8000000000000000};  // 6

const t_bits a_sw_base[]={
  // 0..ld_bits
  // (lo_bit << (1 << sw)) - 1; correct even for sw=ld_bits
  0x0000000000000001,   // 0
  0x0000000000000003,   // 1
  0x000000000000000f,   // 2
  0x00000000000000ff,   // 3
  0x000000000000ffff,   // 4
  0x00000000ffffffff,   // 5
  0xffffffffffffffff};  // 6

const t_bits a_shuffle_mask[]={
  // 0..ld_bits-2
  // For [un]shuffle
  // a_shuffle_mask[i] = a_bfly_mask[i+1] & ~a_bfly_mask[i]
  // => bit_index_swap
  0x2222222222222222,   // 0
  0x0c0c0c0c0c0c0c0c,   // 1
  0x00f000f000f000f0,   // 2
  0x0000ff000000ff00,   // 3
  0x00000000ffff0000};  // 4

const t_bits a_prim_swap[]={
  // 0..ld_bits-1
  // For prim_swap
  // Sum must fill all but highest bit
  // a_prim_swap[i] = a_bfly_lo[i+1] << ((1 << i) - 1)
  0x5555555555555555,   // 0
  0x2222222222222222,   // 1
  0x0808080808080808,   // 2
  0x0080008000800080,   // 3
  0x0000800000008000,   // 5
  0x0000000080000000};  // 6

// eof.
