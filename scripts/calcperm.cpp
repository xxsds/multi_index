//////
// Intro

// Bit permutation code generator

// (c) 2011..2014 by Jasper L. Neumann
// www.sirrida.de / programming.sirrida.de
// E-Mail: info@sirrida.de

// Granted to the public domain
// First version: 2012-08-31 (in Pascal)
// First version: 2013-02-01 (in C++)
// Last change: 2014-02-11

// Compile with
// Gnu C: g++ calcperm.cpp

#include "general.c"
// Choose a working size by including the needed perm_b*.c here:
#include "perm_b64.h"
#include "perm_bas.c"

#include <iostream>
#include <string>
using namespace std;

#define test_count 100000
  // # tests performed

#define nr_hex_digits (bits >> 2)
  // each hexdigit encodes 4 bits
typedef string t_string;


//////
// Porting helpers for Pascal => C++

#define str_ofs (-1)
  // offset for 1-based strings
  // ought to be -1 for C(++), 0 for Pascal
#define maxint 32767
  // ought to be maximum value of t_int
#define assert(x)
  // raise exception if not x
mycall void write(const t_string s) { cout << s; }
mycall void writeln(const t_string s) { cout << s << "\n"; }
mycall void writeln0() { cout << "\n"; }
#define in_set(i,s) (((s>>i)&1)!=0)
  // assuming a simplified set construction
mycall t_int length(const t_string &s) { return s.length(); }
mycall t_string copy(const t_string& s, t_int b, t_int l) { return t_string(s,str_ofs+b,l); }


//////
// Configuration parameter variables

static struct {
  // Output options
  t_bool dump_input;
  t_bool dump_inverse;
  t_bool brief;
  t_bool verbose;

  // Language dependence
  t_string comment_prefix;
  t_string comment_postfix;
  t_string hex_prefix;
  t_string hex_postfix;
  t_string op_assign;
  t_string op_and;
  t_string op_or;
  t_string op_xor;
  t_string op_shl;
  t_string op_shr;

  // Function names
  t_string op_pstep;
  t_string op_pstep_simple;
  t_string op_rol;
  t_string op_gather;
  t_string op_scatter;
  t_string op_bswap;

  // Input options
  t_int in_origin;  // input index origin
  t_int in_base;  // input number base
  t_bool in_indexes_are_target;  // /in_indexes=source or target

  // General costs
  t_int cost_rotate_shift;
  t_int cost_bool;
  t_int cost_bswap;
  t_int cost_mul;
  t_int cost_gs;
  t_int cost_mask;

  // Special costs
  t_int cost_rotate;
  t_int cost_shift;
  t_int cost_and;
  t_int cost_or;
  t_int cost_xor;
  t_int cost_scatter;
  t_int cost_gather;
  t_int cost_bit_permute_step;
  t_int cost_bit_permute_step_simple;

  // Superscalar boni
  t_int bonus_bit_permute_step;
  t_int bonus_bit_permute_step_simple;
  t_int bonus_gs;
  t_int bonus_mask_rol;
  t_int bonus_gs_rol;

  // Calculation options
  t_bool allow_bswap;
  t_bool allow_bmi;
  t_bool test_bpc;
  t_bool test_bfly;
  t_bool test_ibfly;
  t_bool test_benes;
  t_bool test_bit_groups;
  t_bool test_mul;
  t_bool test_gather_scatter;
  t_bool test_gather_shift;
  t_bool test_gather_shift_sloppy;
  t_bool test_shift_scatter;
  t_bool test_shift_scatter_sloppy;
  t_bool test_sag;
  t_bool opt_gs;
  t_bool opt_rol;
  t_bool opt_rol_ex;
  t_bool opt_bswap;

  t_bool self_test;
  } options;

mycall void error_abort() {

  writeln("ERROR");
  exit(1);
  }


//////
// Auxiliary routines

mycall t_longint min(t_longint a, t_longint b) {

  if (a < b) {
    return a;
    }
  else {
    return b;
    }
  }

mycall t_longint max(t_longint a, t_longint b) {

  if (a > b) {
    return a;
    }
  else {
    return b;
    }
  }

mycall t_bits cl_mul(t_bits x, t_bits y) {
// var
  t_bits res;

  res = 0;
  while (x != 0) {
    if (odd(x)) {
      res = res ^ y;
      }
    y = y << 1;
    x = x >> 1;
    }
  return res;
  }

mycall t_bool mul_is_carryless(t_bits x, t_bits y) {
// var
  t_bool res;
  t_bits q;

  res = true;
  q = 0;
  while (x != 0) {
    if (odd(x)) {
      if ((q & y) != 0) {
        res = false;
        break;
        }
      q = q | y;
      }
    y = y << 1;
    x = x >> 1;
    }
  return res;
  }

mycall t_bits carries_of_mul(t_bits x, t_bits y) {
// var
  t_bits res;
  t_bits q;

  res = 0;
  q = 0;
  while (x != 0) {
    if (odd(x)) {
      res = res | (q & y);
      q = q | y;
      }
    y = y << 1;
    x = x >> 1;
    }
  return res;
  }

mycall t_bits mask_hi(t_bits x) {
// Mask all but highest 0 bit.
// var
  t_bits res;

  do {
    res = x;
    x = x | (x+1);  // turn on rightmost bit
    } while (!(x == all_bits));
  return res;
  }

mycall t_bits pdep(t_bits x, t_bits m) {
// scatter

  return expand_right(x,m,ld_bits);
  }

mycall t_bits pext(t_bits x, t_bits m) {
// gather

  return compress_right(x,m,ld_bits);
  }

mycall t_char locase(t_char x) {

  if ( ('A'<=x) &&
       (x<='Z')) {
    x = (t_char)((t_int)(x)-((t_int)('a')-(t_int)('A')));
    }
  return x;
  }

mycall t_char upcase(t_char x) {

  if ( ('a'<=x) &&
       (x<='z')) {
    x = (t_char)((t_int)(x)-((t_int)('A')-(t_int)('a')));
    }
  return x;
  }

mycall t_string lostr(const t_string x) {
// var
  t_int i;
  t_string res;

  res = "";
  for (i = 1; i<=length(x); ++i) {
    res = res + locase(x[str_ofs+i]);
    }
  return res;
  }

mycall t_string upstr(const t_string x) {
// var
  t_int i;
  t_string res;

  res = "";
  for (i = 1; i<=length(x); ++i) {
    res = res + upcase(x[str_ofs+i]);
    }
  return res;
  }

mycall t_string num(t_int x) {
// Convert an integer to a string containing the decimal representation.

  t_string res;
  t_bool neg;

  // num = inttostr(x);
  if (x == 0) {
    res = "0";
    }
  else {
    neg = x < 0;
    if (neg) {
      x = -x;
      }
    res = "";
    while (x != 0) {
      res = (t_char)((x % 10)+(t_int)('0'))+res;
      x = x / 10;
      }
    if (neg) {
      res = "-"+res;
      }
    }
  return res;
  }

mycall t_string hex(t_bits x) {
// Convert a positive integer to a string
// containing the hexadecimal representation.
static const
  t_char digits[16+1] = "0123456789abcdef";  // +1 for 0 char
// var
  t_string res;
  t_int i;

  // hex = hex_prefix+inttohex(x,nr_hex_digits)+hex_postfix;
  res = "";
  for (i = 1; i<=nr_hex_digits; ++i) {
    res = digits[x & 15]+res;
    x = x >> 4;
    }
  return options.hex_prefix+res;
  }

mycall t_string q2s(t_bool x) {
  if (x)  {
    return "true";
    }
  else {
    return "false";
    }
  }

mycall t_longint s2i(const t_string s) {
// Convert a string to an integer.
// Simply halts the program in case of an error
// since we use no exception mechanism for simplicity.
// var
  t_int len;
  t_int i;
  t_int q;
  t_longint v;
  t_bool neg;

  len = length(s);
  // Trim trailing spaces
  while ( (len > 0) &&
          (s[str_ofs+len] == ' ') ) {
    len = len - 1;
    }
  q = 1;
  // Trim leading spaces
  while ( (q <= len) &&
          (s[str_ofs+q] == ' ') ) {
    q = q + 1;
    }
  neg = false;
  // Accept one "-" char
  if ( (q <= len) &&
       (s[str_ofs+q] == '-') ) {
    neg = true;
    q = q + 1;
    // Trim more spaces
    while ( (q <= len) &&
            (s[str_ofs+q] == ' ') ) {
      q = q + 1;
      }
    }
  v = 0;
  for (i = q; i<=len; ++i) {
    if (s[str_ofs+i] == '_') {
      // Skip "_" chars
      }
    else if ( (s[str_ofs+i] >= '0') &&
              (s[str_ofs+i] <= '9') ) {
      v = v*10 + ((t_int)(s[str_ofs+i]) - (t_int)('0'));
      }
    else {
      writeln("ERROR: Number expected in "+s);
      error_abort();
      }
    }
  if (neg) {
    v = -v;
    }
  return v;
  }

#define max_pprim bits
#define max_pstep (ld_bits*2+3)
  // one each extra rol/bswap, one spare
// type

  typedef enum {
    pk_none,
    pk_permute,
    pk_permute_simple,
    pk_mask,
    pk_mask_rol,
    pk_mask_shift,
    pk_rol,
    pk_mul,
    pk_gather_scatter,
    pk_gather,
    pk_gather_shift,
    pk_scatter,
    pk_shift_scatter,
    pk_bswap
    } tq_prim_kind;

  typedef struct {
    // one permutation substep, these are ORed together
    tq_prim_kind kind;
    t_bits mask;   // permute, masked_rol, gather
    t_bits mask2;  // scatter, mul
    t_bits mul;  // mul
    t_int rol;  // rotate/shift count
    t_int rol2;  // mul (not yet used)
    } tr_pprim;

  typedef struct {
    // one permutation step
    t_string description;
    t_bits needed_src_bits;
    t_bits resulting_src_bits;
    t_int nr_pprim;
    tr_pprim a_pprim[max_pprim];
    } tr_pstep;

  typedef struct {
    // complete permutation
    ta_index perm;
    t_string description;
    t_int nr_step;
    tr_pstep a_step[max_pstep];
    t_bits needed_src_bits;
    t_bool is_bfly;  // used in route_benes
    } tr_imp;

  typedef struct {
    // for (routing gather/scatter and derivatives
    t_bits mask_src;
    t_bits mask_tgt;
    } tr_gather_scatter;

  typedef tr_gather_scatter tar_gather_scatter[bits];

  typedef struct {
    t_int cost;
    t_int cmask;
    } tr_performance;

// var
  tr_imp pre_imp;
  tr_imp post_imp;

mycall void imp_dump(const tr_imp &self);  // forward
mycall void dump_perm(const ta_index &perm);  // forward


//////
// Evaluating routines

mycall t_bits pprim_eval(const tr_pprim &self, t_bits x) {
// var
  t_bits res;

  switch (self.kind) {
    case pk_none: res = x;  break;
    case pk_permute: res = bit_permute_step(x, self.mask, self.rol);  break;
    case pk_permute_simple: res = bit_permute_step_simple(x, self.mask, self.rol);  break;
    case pk_rol: res = rol(x, self.rol);  break;
    case pk_mask: res = x & self.mask;  break;
    case pk_mask_rol: res = rol(x & self.mask, self.rol);  break;
    case pk_mask_shift: {
      if (self.rol >= 0) {
        res = (x & self.mask) << self.rol;
        }
      else {
        res = (x & self.mask) >> (-self.rol);
        }
      break;
      }
    case pk_mul: res = rol((t_bits)((rol(x, self.rol) & self.mask)*self.mul) & self.mask2, self.rol2);  break;
    case pk_gather_scatter: res = pdep(pext(x, self.mask), self.mask2);  break;
    case pk_gather: res = pext(x, self.mask);  break;
    case pk_gather_shift: res = (pext(x, self.mask) << self.rol);  break;
    case pk_scatter: res = pdep(x, self.mask2);  break;
    case pk_shift_scatter: res = pdep(x >> ((-self.rol) & (bits-1)), self.mask2);  break;
    case pk_bswap: res = bswap(x);  break;
    default: {
      res = 0;
      writeln("ERROR: pprim_eval: unknown kind"+num((t_int)(self.kind)));
      error_abort();
      break;
      }
    }
  return res;
  }

mycall t_bits pstep_eval(const tr_pstep &self, t_bits x) {
// var
  t_int i;
  t_bits res;

  if (self.nr_pprim == 0) {
    res = x;
    }
  else {
    res = pprim_eval(self.a_pprim[0], x);
    for (i = 1; i<=self.nr_pprim-1; ++i) {
      res = res | pprim_eval(self.a_pprim[i], x);
      }
    }
  return res;
  }

mycall t_bits imp_eval(const tr_imp &self, t_bits x) {
// var
  t_int i;
  t_bits res;

  res = x;
  for (i = 0; i<=self.nr_step-1; ++i) {
    res = pstep_eval(self.a_step[i], res);
    }
  return res;
  }

mycall void imp_check(const tr_imp &imp) {
// var
  t_int q;
  t_bits mask;
  t_bits x;
  t_bool error;
  ta_index inv_perm;

  error = false;
  mask = 1;
  for (q = 0; q<=bits-1; ++q) {
    if (imp.perm[q] != no_index) {
      x = lo_bit << imp.perm[q];
      x = imp_eval(imp, x);
      if (x != mask) {
        error = true;
        break;
        }
      }
    mask = mask << 1;
    }
  if (error) {
    writeln("ERROR (mask)");
    dump_perm(imp.perm);
    imp_dump(imp);
    invert_perm(imp.perm,inv_perm);
    mask = 1;
    for (q = 0; q<=bits-1; ++q) {
      x = imp_eval(imp, mask);
      write(hex(mask)+" => "+hex(x));
      if ( (inv_perm[q] != no_index) &&
           (lo_bit << inv_perm[q] != x) ) {
        writeln(" ERROR, expected: " + hex(lo_bit << inv_perm[q]));
        }
      else {
        writeln(" OK");
        }
      mask = mask << 1;
      }
    error_abort();
    }
  if (imp.nr_step != 0) {
    if (imp_eval(imp, imp.needed_src_bits) !=
        imp.a_step[imp.nr_step-1].resulting_src_bits) {
      writeln("ERROR (needed bits defective)");
      dump_perm(imp.perm);
      imp_dump(imp);
      error_abort();
      }
    }
  }


//////
// Init routines

mycall void init_gather_scatter(tr_gather_scatter &self) {

  self.mask_src = 0;
  self.mask_tgt = 0;
  }

mycall void init_pprim(tr_pprim &self) {

  self.kind = pk_none;
  self.mask = 0;
  self.mask2 = 0;
  self.mul = 0;
  self.rol = 0;
  self.rol2 = 0;
  }

mycall void init_pstep(tr_pstep &self) {
// var
  t_int i;

  self.description = "";
  self.needed_src_bits = 0;
  self.resulting_src_bits = 0;
  self.nr_pprim = 0;
  for (i = 0; i<=max_pprim-1; ++i) {
    init_pprim(self.a_pprim[i]);
    }
  }

mycall void init_imp0(tr_imp &self, const ta_index &perm) {
// var
  t_bits needed;
  t_int i;

  needed = used_source_bits(perm);
  for (i=0; i<=bits-1; ++i) {
    self.perm[i] = perm[i];
    }
  self.description = "";
  self.nr_step = 0;
  for (i = 0; i<=max_pstep-1; ++i) {
    init_pstep(self.a_step[i]);
    }
  self.needed_src_bits = needed;
  self.is_bfly = false;
  self.a_step[0].needed_src_bits = needed;
  }

mycall void init_imp(tr_imp &self, ta_index &perm) {
// var
  ta_index perm1;
  t_int q;
  t_bits x;

  self = pre_imp;
  for (q = 0; q<=bits-1; ++q) {
    perm[q] = no_index;
    perm1[q] = no_index;
    }
  for (q = 0; q<=bits-1; ++q) {
    if (pre_imp.perm[q] != no_index) {
      x = lo_bit << pre_imp.perm[q];
      x = imp_eval(pre_imp, x);
      perm1[q] = nr_trailing_0bits(x);
      assert(x == lo_bit << perm1[q]);
      }
    }
  for (q = 0; q<=bits-1; ++q) {
    x = lo_bit << q;
    x = imp_eval(post_imp, x);
    perm[q] = perm1[nr_trailing_0bits(x)];
    }
  }


//////
// Count mask routines

mycall t_int pprim_cmask(const tr_pprim &self) {
// var
  t_int res;

  switch (self.kind) {
    case pk_none:           res = 0;  break;
    case pk_permute:        res = 1;  break;
    case pk_permute_simple: res = 1;  break;
    case pk_rol:            res = 0;  break;
    case pk_mask:           res = 1;  break;
    case pk_mask_rol:       res = 1;  break;
    case pk_mask_shift:     res = 1;  break;
    case pk_mul:            res = 3;  break;
    case pk_gather_scatter: res = 2;  break;
    case pk_gather:         res = 1;  break;
    case pk_gather_shift:   res = 1;  break;
    case pk_scatter:        res = 1;  break;
    case pk_shift_scatter:  res = 1;  break;
    case pk_bswap:          res = 0;  break;
    default: {
      writeln("ERROR: pprim_cmask: unknown kind");
      res = 0;
      error_abort();
      break;
      }
    }
  return res;
  }

mycall t_int pstep_cmask(const tr_pstep &self) {
// var
  t_int i;
  t_int res;

  res = 0;
  for (i = 0; i<=self.nr_pprim-1; ++i) {
    res = res + pprim_cmask(self.a_pprim[i]);
    }
  return res;
  }

mycall t_int imp_cmask(const tr_imp &self) {
// var
  t_int i;
  t_int res;

  res = 0;
  for (i = 0; i<=self.nr_step-1; ++i) {
    res = res + pstep_cmask(self.a_step[i]);
    }
  return res;
  }


//////
// Cost routines

mycall t_int pprim_cost(const tr_pprim &self) {
// var
  t_int res;

  switch (self.kind) {
    case pk_none:           res = 0;                                         break;
    case pk_permute:        res = options.cost_bit_permute_step;             break;
    case pk_permute_simple: res = options.cost_bit_permute_step_simple;      break;
    case pk_rol:            res = options.cost_rotate;                       break;
    case pk_mask:           res = options.cost_and;                          break;
    case pk_mask_rol:       res = options.cost_and+options.cost_rotate;      break;
    case pk_mask_shift:     res = options.cost_and+options.cost_shift;       break;
    case pk_mul: {
      res = options.cost_mul + options.cost_and * 2;
      if (self.rol != 0) {
        res = res + options.cost_rotate;
        }
      if (self.rol2 != 0) {
        res = res + options.cost_rotate;
        }
      break;
      }
    case pk_gather_scatter: res = options.cost_scatter+options.cost_gather;  break;
    case pk_gather:         res = options.cost_gather;                       break;
    case pk_gather_shift:   res = options.cost_gather+options.cost_shift;    break;
    case pk_scatter:        res = options.cost_scatter;                      break;
    case pk_shift_scatter:  res = options.cost_shift+options.cost_scatter;   break;
    case pk_bswap:          res = options.cost_shift+options.cost_bswap;     break;
    default: {
      writeln("ERROR: pprim_cost: unknown kind");
      res = 0;
      error_abort();
      break;
      }
    }
  return res + pprim_cmask(self) * options.cost_mask;
  }

mycall t_int pstep_cost(const tr_pstep &self) {
const
  t_int s_gs = (
    (1 << pk_gather_scatter) |
    (1 << pk_gather) |
    (1 << pk_gather_shift) |
    (1 << pk_scatter) |
    (1 << pk_shift_scatter) );
const
  t_int s_mr = (
    (1 << pk_rol) |
    (1 << pk_mask) |
    (1 << pk_mask_rol) |
    (1 << pk_mask_shift) |
    (1 << pk_mul) );
// var
  t_int i;
  t_int res;

  if (self.nr_pprim == 0) {
    res = 0;
    }
  else {
    res = pprim_cost(self.a_pprim[0]);
    for (i = 1; i<=self.nr_pprim-1; ++i) {
      res = res + options.cost_or + pprim_cost(self.a_pprim[i]);
      if ( in_set(self.a_pprim[i-1].kind, s_gs) &&
           in_set(self.a_pprim[i].kind, s_gs)) {
        res = res - options.bonus_gs;
        }
      else if ( in_set(self.a_pprim[i-1].kind, s_mr) &&
                in_set(self.a_pprim[i].kind, s_mr)) {
        res = res - options.bonus_mask_rol;
        }
      else if ( in_set(self.a_pprim[i-1].kind, s_mr) &&
                in_set(self.a_pprim[i].kind, s_gs)) {
        res = res - options.bonus_gs_rol;
        }
      else if ( in_set(self.a_pprim[i-1].kind, s_gs) &&
                in_set(self.a_pprim[i].kind, s_mr)) {
        res = res - options.bonus_gs_rol;
        }
      }
    switch (self.a_pprim[0].kind) {
      case pk_permute:        res = res - options.bonus_bit_permute_step;  break;
      case pk_permute_simple: res = res - options.bonus_bit_permute_step_simple;  break;
      default: break;
      }
    }
  return res;
  }

mycall t_int imp_cost(const tr_imp &self) {
// var
  t_int i;
  t_int res;

  res = 0;
  for (i = 0; i<=self.nr_step-1; ++i) {
    res = res + pstep_cost(self.a_step[i]);
    }
  return res;
  }


//////
// Performance routines

mycall void init_performance(tr_performance &self) {

  self.cost = maxint;
  self.cmask = maxint;
  }

mycall t_int cmp_performance(const tr_performance &a, const tr_performance &b) {
// -1: a<b
// 0: a==b
// 1: a>b
// var
  t_int res;

  if (a.cost < b.cost) {
    res = -1;
    }
  else if (a.cost > b.cost) {
    res = 1;
    }
  else if (a.cmask < b.cmask) {
    res = -1;
    }
  else if (a.cmask > b.cmask) {
    res = 1;
    }
  else {
    res = 0;
    }
  return res;
  }

mycall void imp_performance(const tr_imp &imp, tr_performance &perf) {

  perf.cost = imp_cost(imp);
  perf.cmask = imp_cmask(imp);
  }


//////
// Dumping routines

// var
  t_int progress_pos;

mycall void progress() {
static const
  t_string chars[4]={"-","\\","|","/"};

  progress_pos = progress_pos+1;
  write(chars[progress_pos & 3]);
  write("\x0d");
  }

mycall void dump_perm(const ta_index &perm) {
// var
  t_int i;

  for (i = 0; i<=bits-1; ++i) {
    if (perm[i] == no_index) {
      write("* ");
      }
    else {
      write(num(perm[i])+" ");
      }
    }
  writeln0();
  }

mycall t_string dump_performance(const tr_performance &perf) {

  return
    num(perf.cost)+
    " cycles, "+
    num(perf.cmask)+
    " masks";
  }

mycall t_string pprim_dump(const tr_pprim &self) {
// var
  t_string res;

  switch (self.kind) {
    case pk_none: {
      res = "";
      break;
      }
    case pk_permute: {
      res = options.op_pstep+"(x, "+hex(self.mask)+", "+num(self.rol)+')';
      break;
      }
    case pk_permute_simple: {
      res = options.op_pstep_simple+
        "(x, "+hex(self.mask)+", "+num(self.rol)+')';
      break;
      }
    case pk_rol: {
      res = options.op_rol+"(x, "+num(self.rol)+')';
      break;
      }
    case pk_mask: {
      res = "(x "+options.op_and+' '+hex(self.mask)+')';
      break;
      }
    case pk_mask_rol: {
      res = options.op_rol+
        "(x "+options.op_and+' '+hex(self.mask)+", "+num(self.rol)+')';
      break;
      }
    case pk_mask_shift: {
      if (self.rol >= 0) {
        res = "((x "+options.op_and+' '+hex(self.mask)+") "+
          options.op_shl+' '+num(self.rol)+')';
        }
      else {
        res = "((x "+options.op_and+' '+hex(self.mask)+") "+
          options.op_shr+' '+num(-self.rol)+')';
        }
      break;
      }
    case pk_mul: {
      if (self.rol==0) {
        if (self.rol2==0) {
          res = "(((x "+options.op_and+' '+hex(self.mask)+") * "+
            hex(self.mul)+") "+options.op_and+' '+hex(self.mask2)+')';
          }
        else {
          res = options.op_rol+
            "(((x "+options.op_and+' '+hex(self.mask)+") * "+
            hex(self.mul)+") "+options.op_and+' '+hex(self.mask2)+", "+
            num(self.rol2 & (bits-1))+')';
          }
        }
      else {
        if (self.rol2==0) {
          res = "((("+options.op_rol+"(x, "+num(self.rol)+") "+
            options.op_and+' '+hex(self.mask)+") * "+hex(self.mul)+") "+
            options.op_and+' '+hex(self.mask2)+')';
          }
        else {
          res = options.op_rol+
            "((("+options.op_rol+"(x, "+num(self.rol)+") "+
            options.op_and+' '+hex(self.mask)+") * "+hex(self.mul)+") "+
            options.op_and+' '+hex(self.mask2)+", "+
            num(self.rol2 & (bits-1))+')';
          }
        }
      break;
      }
    case pk_gather_scatter: {
      res = options.op_scatter+
        '('+options.op_gather+"(x, "+hex(self.mask)+"), "+hex(self.mask2)+')';
      break;
      }
    case pk_gather: {
      res = options.op_gather+"(x, "+hex(self.mask)+')';
      break;
      }
    case pk_gather_shift: {
      res = '('+options.op_gather+
        "(x, "+hex(self.mask)+") "+options.op_shl+' '+num(self.rol)+')';
      break;
      }
    case pk_scatter: {
      res = options.op_scatter+"(x, "+hex(self.mask2)+')';
      break;
      }
    case pk_shift_scatter: {
      res = options.op_scatter+"(x "+options.op_shr+' '+
        num((-self.rol) & (bits-1))+", "+hex(self.mask2)+')';
      break;
      }
    case pk_bswap: {
      res = options.op_bswap+"(x)";
      break;
      }
    default: {
      res = "";
      writeln("ERROR: pprim_dump: unknown kind");
      error_abort();
      break;
      }
    }
  return res;
  }

mycall void pstep_dump(const tr_pstep &self) {
// var
  t_int i;

  if (self.nr_pprim != 0) {
    std::cerr << "x " << options.op_assign << " " << pprim_dump(self.a_pprim[0]);
    for (i = 1; i<=self.nr_pprim-1; ++i) {
      std::cerr << std::endl;
      std::cerr << "  " << options.op_or << " " << pprim_dump(self.a_pprim[i]);
      }
    std::cerr << ";";
    }
  }

mycall void imp_dump_simple(const tr_imp &self) {
// var
  tr_performance perf;

  imp_performance(self, perf);
  writeln(
    self.description+
    ": "+
    dump_performance(perf));
  }

mycall void imp_dump(const tr_imp &self) {
// var
  t_int i;

  imp_dump_simple(self);
  
  for (i = 0; i<=self.nr_step-1; ++i) {
    pstep_dump(self.a_step[i]);
    if (self.a_step[i].description != "") {
      std::cerr << "  " << options.comment_prefix << self.a_step[i].description << options.comment_postfix;
      }
      std::cerr << std::endl;
    }
  }


//////
// Code construction routines

mycall void concat_description(t_string &tgt, const t_string &src) {

  if (tgt == "") {
    tgt = src;
    }
  else if (src != "") {
    tgt = tgt + " + " + src;
    }
  }

mycall void finish_step(tr_imp &imp) {
// var
  t_bits needed;
  t_bits resulting;

  if (imp.nr_step == 0) {
    needed = imp.needed_src_bits;
    }
  else {
    needed = imp.a_step[imp.nr_step-1].resulting_src_bits;
    }
  resulting = pstep_eval(imp.a_step[imp.nr_step], needed);
  imp.a_step[imp.nr_step].resulting_src_bits = resulting;
  imp.nr_step = imp.nr_step+1;
  imp.a_step[imp.nr_step].needed_src_bits = resulting;
  }

mycall void finish_perm(tr_imp &imp) {
// var
  t_int i;

  for (i = 0; i<=post_imp.nr_step-1; ++i) {
    imp.a_step[imp.nr_step] = post_imp.a_step[i];
    finish_step(imp);
    }
  concat_description(imp.description, post_imp.description);
  }

mycall void add_permute(tr_imp &imp, t_bits mask, t_int rol) {
// Add a permute step
// var
  t_bits n;

  if (mask != 0) {
    rol = rol & (bits-1);
    imp.a_step[imp.nr_step].a_pprim[0].rol = rol;
    imp.a_step[imp.nr_step].a_pprim[0].mask = mask;
    n = imp.a_step[imp.nr_step].needed_src_bits;
    if (nr_1bits(bit_permute_step_simple(n,mask,rol)) == nr_1bits(n)) {
      // bit_permute_step_simple does not kill needed bits
      imp.a_step[imp.nr_step].a_pprim[0].kind = pk_permute_simple;
      }
    else {
      // we need a full permute step
      imp.a_step[imp.nr_step].a_pprim[0].kind = pk_permute;
      }
    imp.a_step[imp.nr_step].nr_pprim = 1;
    finish_step(imp);
    }
  }

mycall void make_rol_step(tr_pstep &pstep, t_bits mask, t_int rol) {
// Add a mask/rol step

  rol = rol & (bits-1);
  pstep.a_pprim[pstep.nr_pprim].rol = rol;
  pstep.a_pprim[pstep.nr_pprim].mask = mask;

  if (mask == all_bits) {
    // nothing to mask
    pstep.a_pprim[pstep.nr_pprim].kind = pk_rol;
    }
  else if (rol == 0) {
    // nothing to rotate
    pstep.a_pprim[pstep.nr_pprim].kind = pk_mask;
    }
  else if (((t_bits)(mask << rol) >> rol) == mask) {
    // shift left is sufficient
    pstep.a_pprim[pstep.nr_pprim].kind = pk_mask_shift;
    }
  else if (((t_bits)(mask >> (bits-rol)) << (bits-rol)) == mask) {
    // shift right is sufficient
    pstep.a_pprim[pstep.nr_pprim].kind = pk_mask_shift;
    pstep.a_pprim[pstep.nr_pprim].rol = rol-bits;  // < 0
    }
  else {
    // we need a rotate; a shift is not enough
    pstep.a_pprim[pstep.nr_pprim].kind = pk_mask_rol;
    }
  }

mycall void add_rol_step(tr_pstep &pstep, t_bits mask, t_int rol) {

  rol = rol & (bits-1);
  if ( (mask != 0) &&
       ( (mask != all_bits) ||
         (rol != 0) ) ) {
    make_rol_step(pstep, mask, rol);
    pstep.nr_pprim = pstep.nr_pprim+1;
    }
  }

mycall void add_rol(tr_imp &imp, t_bits mask, t_int rol) {

  add_rol_step(imp.a_step[imp.nr_step], mask, rol);
  }

mycall void add_bswap_step(tr_pstep &pstep) {

  pstep.a_pprim[pstep.nr_pprim].kind = pk_bswap;
  pstep.nr_pprim = pstep.nr_pprim+1;
  }

mycall void add_bswap(tr_imp &imp) {

  add_bswap_step(imp.a_step[imp.nr_step]);
  }

mycall void add_gs_step(tr_pstep &pstep, t_bits mask, t_bits mask2) {
// Add a gather/scatter step
// var
  t_int t1,t2;

  if (mask != 0) {
    pstep.a_pprim[pstep.nr_pprim].rol = 0;
    pstep.a_pprim[pstep.nr_pprim].mask = mask;
    pstep.a_pprim[pstep.nr_pprim].mask2 = mask2;
    t1 = nr_trailing_0bits(mask);
    t2 = nr_trailing_0bits(mask2);
    if (mask == mask2) {
      // we scatter what we gather, only masking needed
      pstep.a_pprim[pstep.nr_pprim].kind = pk_mask;
      }
    else if ((mask >> t1) == (mask2 >> t2)) {
      // the masks differ by a shift, so mask and rotate
      make_rol_step(pstep, mask, t2-t1);
      }
    else if (is_contiguous_1bits(mask2)) {
      // the target bits are contiguous
      if (t2 == 0) {
        // single gather, no shift necessary
        pstep.a_pprim[pstep.nr_pprim].kind = pk_gather;
        }
      else {
        // gather and shift
        pstep.a_pprim[pstep.nr_pprim].kind = pk_gather_shift;
        pstep.a_pprim[pstep.nr_pprim].rol = t2;
        }
      }
    else if (is_contiguous_1bits(mask)) {
      // the source bits are contiguous
      if (t1 == 0) {
        // single scatter, no shift necessary
        pstep.a_pprim[pstep.nr_pprim].kind = pk_scatter;
        }
      else {
        // shift to 0 position and scatter
        pstep.a_pprim[pstep.nr_pprim].kind = pk_shift_scatter;
        pstep.a_pprim[pstep.nr_pprim].rol = (-t1) & (bits-1);
        }
      }
    else {
      // last resort: gather and scatter
      pstep.a_pprim[pstep.nr_pprim].kind = pk_gather_scatter;
      }
    pstep.nr_pprim = pstep.nr_pprim+1;
    }
  }

mycall void add_gs(tr_imp &imp, t_bits mask, t_bits mask2) {

  add_gs_step(imp.a_step[imp.nr_step], mask, mask2);
  }

mycall void add_sag_step(tr_pstep &pstep, t_bits mask, t_int rol) {
// Add a sheep and goats step
// var
  t_bits lo_mask;

  rol = rol & (bits-1);
  lo_mask = (lo_bit << rol)-1;
  add_gs_step(pstep, mask, lo_mask << rol);
  add_gs_step(pstep, (~mask) & (lo_mask+(lo_mask << rol)), lo_mask);
  }

mycall void add_sag(tr_imp &imp, t_bits mask, t_int rol) {

  add_sag_step(imp.a_step[imp.nr_step], mask, rol);
  finish_step(imp);
  }

mycall void add_mul_step(tr_pstep &pstep, t_bits mask, t_bits mul, t_bits mask2, t_int rol, t_int rol2) {
// Add a mul step

  if (mask != 0) {
    pstep.a_pprim[pstep.nr_pprim].kind = pk_mul;
    pstep.a_pprim[pstep.nr_pprim].rol = rol & (bits-1);
    pstep.a_pprim[pstep.nr_pprim].rol2 = rol2 & (bits-1);
    pstep.a_pprim[pstep.nr_pprim].mask = mask;
    pstep.a_pprim[pstep.nr_pprim].mask2 = mask2;
    pstep.a_pprim[pstep.nr_pprim].mul = mul;
    pstep.nr_pprim = pstep.nr_pprim+1;
    }
  }


//////
// Routing routines

mycall t_bool is_identity() {
// var
  tr_imp imp;
  ta_index perm;
  t_int i;
  t_bool res;

  res = true;
  init_imp(imp,perm);
  for (i = 0; i<=bits-1; ++i) {
    if (perm[i] != no_index) {
      if (i != perm[i]) {
        res = false;
        }
      }
    }
  return res;
  }

mycall void route_benes(tr_imp &imp, const ta_subword &a_stage) {
// var
  ta_index perm;
  tr_benes benes;
  t_int stage_idx,stage;

  init_imp(imp,perm);
  gen_benes_ex(&benes,perm,a_stage);
  concat_description(imp.description, "Benes ");
  imp.is_bfly = true;
  for (stage_idx = 0; stage_idx<=ld_bits-1; ++stage_idx) {
    stage = a_stage[stage_idx];
    imp.description = imp.description+num(stage);
    if (benes.b1.cfg[stage] != 0) {
      imp.is_bfly = false;
      add_permute(imp, benes.b1.cfg[stage], 1 << stage);
      }
    }
  for (stage_idx = ld_bits-1; stage_idx>=0; --stage_idx) {
    stage = a_stage[stage_idx];
    if (benes.b2.cfg[stage] != 0) {
      add_permute(imp, benes.b2.cfg[stage], 1 << stage);
      }
    }
  finish_perm(imp);
  }

mycall void route_bit_groups(tr_imp &imp) {
// var
  ta_index perm;
  t_int i,j,k,count,min_pos;
  t_bits s;  // set) { differences
  t_bits masks[bits];

  init_imp(imp,perm);
  concat_description(imp.description, "Bit group moving");

  s = 0;
  for (i = 0; i<=bits-1; ++i) {
    if (perm[i] != no_index) {
      s = s | (lo_bit << ((i-perm[i]) & (bits-1)));
      }
    }

  count = 0;
  min_pos = 0;
  for (i = bits-1; i>=0; --i) {
    if (((lo_bit << i) & s) != 0) {
      count = count+1;
      min_pos = i;
      }
    }

  switch (count) {
    case 0: {
      finish_perm(imp);
      return;
      break;
      }
    case 1: {
      if (min_pos == 0) {
        finish_perm(imp);
        return;
        }
      add_rol(imp, all_bits, min_pos);
      finish_step(imp);
      finish_perm(imp);
      return;
      break;
      }
    default: break;
    }

  for (i = 0; i<=bits-1; ++i) {
    masks[i] = 0;
    }
  for (i = 0; i<=bits-1; ++i) {
    if (perm[i] != no_index) {
      j = (i-perm[i]) & (bits-1);
      k = (i-j) & (bits-1);
      masks[j] = masks[j] | (lo_bit << k);
      }
    }

  for (i = 0; i<=bits-1; ++i) {
    add_rol(imp, masks[i], i);
    }

  finish_step(imp);
  finish_perm(imp);
  }

mycall void route_mul1(tr_imp &imp) {
// 2013-03-19  J. Neumann

// var
  ta_index perm;
  t_int i;
  t_int idx;
  t_int best_idx;
  t_int diff;
  t_int rol;
  t_bits u,use,m,mul,tgt;
  t_bool empty;
  t_int pop, best_pop;
  ta_index save_perm;
  t_int nr_pprim;
  t_bool redo;
  t_int my_cost_mul, my_cost_rol;

  init_imp(imp,perm);
  concat_description(imp.description, "Shift by mul (1)");

  nr_pprim = imp.a_step[imp.nr_step].nr_pprim;
  for (i = 0; i<=bits-1; ++i) {
    save_perm[i] = perm[i];
    }

  do {
    redo = false;
    rol = 0;
    while (true) {

      empty = true;
      for (i = 0; i<=bits-1; ++i) {
        if (perm[i] != no_index) {
          empty = false;
          break;
          }
        }
      if (empty) {
        break;
        }

      best_idx = 0;
      best_pop = 0;
      for (idx = 0; idx<=bits-1; ++idx) {

        // Pseudo-calc
        mul = 0;
        use = 0;
        for (i = 0; i<=bits-1; ++i) {
          if (perm[i] != no_index) {
            diff = i-((perm[i]-idx) & (bits-1));
            if (diff < 0) {
              continue;  // Can not be reached by mul
              }
            u = lo_bit << ((perm[i]-idx) & (bits-1));
            m = lo_bit << diff;
            if (((t_bits)(mul*use) & (t_bits)(mul*u)) != 0) {
              continue;  // Conflict
              }
            if ((m & mul) != 0) {
              // already handled
              }
            else {
              // new entry
              if (((t_bits)(mul*(use | u)) & (t_bits)(m*(use | u))) != 0) {
                continue;  // Conflict
                }
              mul = mul | m;
              }
            use = use | u;
            }
          }

        pop = nr_1bits(mul) * 10 + nr_1bits(use);
        if (((idx+rol) & (bits-1)) == 0) {
          pop = pop + 5;  // prefer not to rotate
          }
        if (pop > best_pop) {
          best_idx = idx;
          best_pop = pop;
          }
        }

      for (i = 0; i<=bits-1; ++i) {
        if (perm[i] != no_index) {
          perm[i] = (perm[i]-best_idx) & (bits-1);
          }
        }
      rol = rol + best_idx;

      mul = 0;
      use = 0;
      tgt = 0;
      for (i = 0; i<=bits-1; ++i) {
        if (perm[i] != no_index) {
          diff = i-perm[i];
          if (diff<0) {
            continue;  // Can not be reached by mul
            }
          u = lo_bit << perm[i];
          m = lo_bit << diff;
          if (((t_bits)(mul*use) & (t_bits)(mul*u)) != 0) {
            continue;  // Conflict
            }
          if ((m & mul) != 0) {
            // already handled
            }
          else {
            // new entry
            if (((t_bits)(mul*(use | u)) & (t_bits)(m*(use | u))) != 0) {
              continue;  // Conflict
              }
            mul = mul | m;
            }
          use = use | u;
          tgt = tgt | (lo_bit << i);
          perm[i] = no_index;
          }
        }

      my_cost_mul =
          // +options.cost_rotate*(t_int)(rol!=0)
          +options.cost_rotate
          +options.cost_mul
          +options.cost_and*2
          +options.cost_mask*2
          +options.cost_or
          -options.bonus_mask_rol;
      my_cost_rol =
          // -(t_int)((mul and 1)!=0)*options.cost_rotate
          +nr_1bits(mul)*(
            +options.cost_and
            +options.cost_rotate
            +options.cost_or
            +options.cost_mask
            -options.bonus_mask_rol);
      if (my_cost_mul <= my_cost_rol) {
        add_mul_step(imp.a_step[imp.nr_step], use,mul,tgt,-rol,0);
        }
      else {

        diff = (rol-nr_trailing_0bits(mul)) & (bits-1);
        use = 0;
        for (i = 0; i<=bits-1; ++i) {
          if (save_perm[i] != no_index) {
            if (((save_perm[i]-i) & (bits-1)) == diff) {
              use = use | (lo_bit << ((i+diff) & (bits-1)));
              save_perm[i] = no_index;
              }
            }
          }
        imp.a_step[imp.nr_step].nr_pprim = nr_pprim;
        add_rol(imp, use, -diff);
        nr_pprim = imp.a_step[imp.nr_step].nr_pprim;
        for (i = 0; i<=bits-1; ++i) {
          perm[i] = save_perm[i];
          }

        redo = true;
        break;
        }
      }
    } while (!(!redo));
  finish_step(imp);
  finish_perm(imp);
  }

mycall void route_mul2(tr_imp &imp) {
// 2013-03-22  J. Neumann

// var
  ta_index perm;
  t_int i,j,best_count;
  t_int rotate;
  t_bits mask,u,use,m,mul,used,tgt,use0,u1;
  t_bool found;
  t_int my_cost_mul, my_cost_rol;
  t_bits masks[bits];
  t_int counts[bits];  // nr_1bits(masks[i])

  init_imp(imp,perm);
  concat_description(imp.description, "Shift by mul (2)");

  mask = 0;
  for (i = 0; i<=bits-1; ++i) {
    masks[i] = 0;
    counts[i] = 0;
    }
  for (i = 0; i<=bits-1; ++i) {
    if (perm[i] != no_index) {
      j = (i-perm[i]) & (bits-1);
      u = lo_bit << ((i-j) & (bits-1));
      mask = mask | u;
      masks[j] = masks[j] | u;
      counts[j] = counts[j] + 1;
      }
    }

  while (mask!=0) {

    best_count = 0;
    i = 0;

    for (j = 0; j<=bits-1; ++j) {
      if (masks[j]!=0) {
        if (counts[j]>best_count) {
          best_count = counts[j];
          i = j;
          }
        }
      }

    use = rol(masks[i],i);  // used bits
    rotate = i;
    used = lo_bit << i;  // set of indexes which masks are used

    mul = 1;
    tgt = use;
    use0 = use;

    my_cost_rol = 1;

    do {
      found = false;
      for (i = 0; i<=bits-1; ++i) {
        if ( (masks[i] != 0) &&
             ((used & (lo_bit << i)) == 0) ) {
          u = rol(masks[i],rotate);
          j = ((i-rotate) & (bits-1));
          m = lo_bit << j;
          if (nr_1bits(m*u) == nr_1bits(u)) {
            // if (carries_of_mul(mul | m, use | u) == 0) {
            if ((carries_of_mul(mul | m, use | u) &
                 mask_hi(tgt | (u << j))) == 0) {
              my_cost_rol = my_cost_rol + 1;
              use = use | u;
              mul = mul | m;
              used = used | (lo_bit << i);
              tgt = tgt | (u << j);
              found = true;
              }
            }
          }
        }
      } while (!(!found));

    my_cost_mul =
        +options.cost_rotate*(t_int)(rotate!=0)
        +options.cost_mul
        +options.cost_and*2
        +options.cost_mask*2
        +options.cost_or
        -options.bonus_mask_rol;
    my_cost_rol =
         my_cost_rol*(
          +options.cost_and
          +options.cost_rotate
          +options.cost_or
          +options.cost_mask
          -options.bonus_mask_rol);
    if (my_cost_mul >= my_cost_rol) {
      mask = mask & ~masks[rotate];
      masks[rotate] = 0;
      add_rol(imp, rol(use0,-rotate), rotate);
      }
    else {
      for (i = 0; i<=bits-1; ++i) {
        if ((used & (lo_bit << i)) != 0) {
          mask = mask & ~masks[i];
          masks[i] = 0;
          counts[i] = 0;
          }
        }

      do {
        found = false;
        for (i = bits-1; i>=0; --i) {
          if (masks[i] != 0) {
            u = rol(masks[i],rotate);
            j = ((i-rotate) & (bits-1));
            m = lo_bit << j;
            while (u!=0) {
              u1 = u & (-u);
              if (nr_1bits(m*u1) == nr_1bits(u1)) {
                if ((carries_of_mul(mul | m, use | u1) &
                     mask_hi(tgt | (u1 << j))) == 0) {
                  // my_cost_rol = my_cost_rol + 1;
                  use = use | u1;
                  mul = mul | m;
                  // used = used | (lo_bit << i);
                  tgt = tgt | (u1 << j);
                  mask = mask & ~rol(u1,-rotate);
                  masks[i] = masks[i] & ~rol(u1,-rotate);
                  counts[i] = counts[i] - 1;
                  found = true;
                  }
                }
              u = u & ~u1;
              }
            }
          }
        } while (!(!found));

      add_mul_step(imp.a_step[imp.nr_step], use,mul,tgt,rotate,0);
      }
    }
  finish_step(imp);
  finish_perm(imp);
  }

mycall void route_mul3(tr_imp &imp) {
// 2013-03-19  J. Neumann

// var
  ta_index perm;
  t_int i,j;
  t_int idx;
  t_int best_idx;
  t_int diff;
  t_int rotate;
  t_bits u,use,m,mul,tgt;
  t_bool empty;
  t_int pop, best_pop;
  ta_index save_perm;
  t_int nr_pprim;
  t_bool redo;
  t_int hibit;
  t_bits himask;
  t_int my_cost_mul, my_cost_rol;

  init_imp(imp,perm);
  concat_description(imp.description, "Shift by mul (3)");

  nr_pprim = imp.a_step[imp.nr_step].nr_pprim;
  for (i = 0; i<=bits-1; ++i) {
    save_perm[i] = perm[i];
    }

  do {
    redo = false;
    rotate = 0;
    while (true) {

      empty = true;
      for (i = 0; i<=bits-1; ++i) {
        if (perm[i] != no_index) {
          empty = false;
          break;
          }
        }
      if (empty) {
        break;
        }

      best_idx = 0;
      best_pop = 0;
      for (idx = 0; idx<=bits-1; ++idx) {

        // Pseudo-calc
        mul = 0;
        use = 0;
        hibit = 0;
        for (i = bits-1; i>=0; --i) {
          if (perm[i] != no_index) {
            j = ((perm[i]-idx) & (bits-1));
            diff = i-j;
            if (diff < 0) {
              continue;  // Can not be reached by mul
              }
            u = lo_bit << j;
            m = lo_bit << diff;
            if (i > hibit) {
              j = i;
              }
            else {
              j = hibit;
              }
            himask = (lo_bit << j)*2-1;
            if ((t_bits(mul*use) & t_bits(mul*u) & himask) != 0) {
              continue;  // Conflict
              }
            if ((m & mul) != 0) {
              // already handled
              }
            else {
              // new entry
              if (((t_bits)(mul*(use | u)) & (t_bits)(m*(use | u)) & himask) != 0) {
                continue;  // Conflict
                }
              mul = mul | m;
              }
            use = use | u;
            if (i > hibit) {
              hibit = i;
              }
            }
          }

        pop = nr_1bits(mul) + nr_1bits(use) * 10;
        if (((idx+rotate) & (bits-1)) == 0) {
          pop = pop + 1;  // prefer not to rotate
          }
        if (pop > best_pop) {
          best_idx = idx;
          best_pop = pop;
          }
        }

      for (i = 0; i<=bits-1; ++i) {
        if (perm[i] != no_index) {
          perm[i] = (perm[i]-best_idx) & (bits-1);
          }
        }
      rotate = rotate + best_idx;

      mul = 0;
      use = 0;
      tgt = 0;
      hibit = 0;

      for (i = bits-1; i>=0; --i) {
        if (perm[i] != no_index) {
          diff = i-perm[i];
          if (diff<0) {
            continue;  // Can not be reached by mul
            }
          u = lo_bit << perm[i];
          m = lo_bit << diff;
          if (i > hibit) {
            j = i;
            }
          else {
            j = hibit;
            }
          himask = (lo_bit << j)*2-1;
          if (((t_bits)(mul*use) & (t_bits)(mul*u) & himask) != 0) {
            continue;  // Conflict
            }
          if ((m & mul) != 0) {
            // already handled
            }
          else {
            // new entry
            if (((t_bits)(mul*(use | u)) & (t_bits)(m*(use | u)) & himask) != 0) {
              continue;  // Conflict
              }
            mul = mul | m;
            }
          use = use | u;
          tgt = tgt | (lo_bit << i);
          if (i > hibit) {
            hibit = i;
            }
          perm[i] = no_index;
          }
        }

      my_cost_mul =
          +options.cost_rotate*(t_int)(rotate!=0)
          +options.cost_mul
          +options.cost_and*2
          +options.cost_mask*2
          +options.cost_or
          -options.bonus_mask_rol;
      my_cost_rol =
          +nr_1bits(mul)*(
            +options.cost_and
            +options.cost_rotate
            +options.cost_or
            +options.cost_mask
            -options.bonus_mask_rol);
      if (my_cost_mul <= my_cost_rol) {
        add_mul_step(imp.a_step[imp.nr_step], use,mul,tgt,-rotate,0);
        }
      else {

        best_idx = 0;
        best_pop = 0;
        for (idx = 0; idx<=bits-1; ++idx) {
          if (((lo_bit << idx) & mul) != 0) {
            diff = (rotate-idx) & (bits-1);

            // Pseudo-calc
            pop = 0;
            for (i = 0; i<=bits-1; ++i) {
              if (save_perm[i] != no_index) {
                if (((save_perm[i]-i) & (bits-1)) == diff) {
                  pop = pop + 1;
                  }
                }
              }

            if (pop > best_pop) {
              best_idx = idx;
              best_pop = pop;
              }
            }
          }

        diff = (rotate-best_idx) & (bits-1);
        use = 0;
        for (i = 0; i<=bits-1; ++i) {
          if (save_perm[i] != no_index) {
            if (((save_perm[i]-i) & (bits-1)) == diff) {
              use = use | (lo_bit << ((i+diff) & (bits-1)));
              save_perm[i] = no_index;
              }
            }
          }

        imp.a_step[imp.nr_step].nr_pprim = nr_pprim;
        add_rol(imp, use, -diff);
        nr_pprim = imp.a_step[imp.nr_step].nr_pprim;
        for (i = 0; i<=bits-1; ++i) {
          perm[i] = save_perm[i];
          }

        redo = true;
        break;
        }
      }
    } while (!(!redo));
  finish_step(imp);
  finish_perm(imp);
  }

mycall void route_sag(tr_imp &imp) {
// 2012-09-14 (Knuth)

// var
  ta_index perm;
  struct {
    t_bits mask;
    t_bits needed;
    } a[ld_bits];
  t_int i,j;
  t_bits mask_src;
  t_int max_src,max_tgt;
  t_int my_ld_bits;
  t_int my_bits;

  t_bits needed_src;
  ta_index inv_perm;
  ta_index perm1;
  ta_index perm2;

  init_imp(imp,perm);
  concat_description(imp.description, "Sheep-and-goats method");

  for (i = 0; i<=ld_bits-1; ++i) {
    a[i].mask = 0;
    a[i].needed = 0;
    }
  invert_perm(perm,inv_perm);

  needed_src = 0;  // set of needed source bits
  mask_src = 1;
  max_src = -1;
  max_tgt = -1;
  for (i = 0; i<=bits-1; ++i) {
    if (perm[i] > max_tgt) {
      max_tgt = perm[i];
      }
    if (inv_perm[i] != no_index) {
      needed_src = needed_src | mask_src;
      if (inv_perm[i] > max_src) {
        max_src = inv_perm[i];
        }
      }
    mask_src = mask_src << 1;
    }

  // fill don't care entries ascending with hole indexes
  j = 0;
  for (i = 0; i<=bits-1; ++i) {
    if (perm[i] != no_index) {
      perm1[i] = perm[i];
      }
    else {
      while (((lo_bit << j) & needed_src) != 0) {
        j = j+1;
        }
      perm1[i] = j;
      j = j+1;
      }
    }

  if (max_tgt < 0) {
    finish_perm(imp);
    return;  // nothing to do
    }
  if (max_src > max_tgt) {
    max_tgt = max_src;
    }
  my_ld_bits = 0;
  while (max_tgt != 0) {
    my_ld_bits = my_ld_bits+1;
    max_tgt = max_tgt >> 1;
    }
  my_bits = 1 << my_ld_bits;

  invert_perm(perm1,perm2);
  for (i = 0; i<=my_ld_bits-1; ++i) {
    for (j = 0; j<=my_bits-1; ++j) {
      a[i].mask = a[i].mask | ((t_bits)(((perm2[j] >> i) & 1)) << j);
      }
    a[i].needed = needed_src;
    }

  for (i = 0; i<=my_ld_bits-1; ++i) {
    for (j = i+1; j<=my_ld_bits-1; ++j) {
      a[j].mask =
        compress_right(a[j].mask,~a[i].mask,my_ld_bits) |
        (compress_right(a[j].mask,a[i].mask,my_ld_bits) << (my_bits >> 1));
      a[j].needed =
        compress_right(a[j].needed,~a[i].mask,my_ld_bits) |
        (compress_right(a[j].needed,a[i].mask,my_ld_bits) << (my_bits >> 1));
      }
    add_sag(imp, a[i].mask, my_bits >> 1);
    }

  finish_perm(imp);
  }

typedef
  mycall t_int (*tf_aux_route_gs)(
    tar_gather_scatter &,  // a
    const ta_index &,  // perm
    t_bits,  // needed_src
    t_bits);  // needed_tgt

mycall t_int aux_route_gather_scatter(
  tar_gather_scatter &a,
  const ta_index &perm,
  t_bits needed_src,
  t_bits needed_tgt) {
// FAR;
// 2012-09-27  J. Neumann
// var
  t_int a_count;
  t_int s,t;
  t_bits mask_src;
  t_bits mask_tgt;

  a_count = 0;
  while (needed_src != 0) {
    s = 0;
    t = -1;
    while (true) {
      t = t+1;
      while ( (t < bits) &&
              (((lo_bit << t) & needed_tgt) == 0) ) {
        t = t+1;
        }
      if (t >= bits) {
        break;
        }

      if (perm[t] < s) {
        continue;
        }
      s = perm[t];
      // assert(s != no_index);

      mask_src = lo_bit << s;
      mask_tgt = lo_bit << t;
      a[a_count].mask_src = a[a_count].mask_src | mask_src;
      a[a_count].mask_tgt = a[a_count].mask_tgt | mask_tgt;
      needed_src = needed_src & ~mask_src;
      needed_tgt = needed_tgt & ~mask_tgt;
      }
    a_count = a_count+1;
    }
  return a_count;
  }

mycall t_int aux_route_gather_shift_sloppy(
  tar_gather_scatter &a,
  const ta_index &perm,
  t_bits needed_src,
  t_bits needed_tgt) {
// FAR;
// 2012-09-25  J. Neumann
// var
  t_int a_count;
  t_int s,t;
  t_bits mask_src;
  t_bits mask_tgt;

  a_count = 0;
  while (needed_src != 0) {
    s = 0;  // -bits;
    t = -1;
    while (true) {
      t = t+1;
      while ( (t < bits) &&
              (((lo_bit << t) & needed_tgt) == 0) ) {
        t = t+1;
        // s = s+1;
        }
      if (t >= bits) {
        break;
        }

      if (perm[t] < s) {
        break;
        }
      s = perm[t];
      // assert(s != no_index);

      mask_src = lo_bit << s;
      mask_tgt = lo_bit << t;
      a[a_count].mask_src = a[a_count].mask_src | mask_src;
      a[a_count].mask_tgt = a[a_count].mask_tgt | mask_tgt;
      needed_src = needed_src & ~mask_src;
      needed_tgt = needed_tgt & ~mask_tgt;
      }
    a_count = a_count+1;
    }
  return a_count;
  }

mycall t_int aux_route_gather_shift(
  tar_gather_scatter &a,
  const ta_index &perm,
  t_bits needed_src,
  t_bits needed_tgt) {
// FAR;
// 2012-09-26  J. Neumann
// var
  t_int a_count;
  t_int s,t;
  t_bits mask_src;
  t_bits mask_tgt;

  a_count = 0;
  while (needed_src != 0) {
    s = 0;  // -bits;
    t = 0;
    while ( (t < bits) &&
            (((lo_bit << t) & needed_tgt) == 0) ) {
      t = t+1;
      }
    t = t-1;
    while (true) {
      t = t+1;
      if (t >= bits) {
        break;
        }

      if (perm[t] == no_index) {
        break;
        }
      if (perm[t] < s) {
        break;
        }
      s = perm[t];

      mask_src = lo_bit << s;
      mask_tgt = lo_bit << t;
      a[a_count].mask_src = a[a_count].mask_src | mask_src;
      a[a_count].mask_tgt = a[a_count].mask_tgt | mask_tgt;
      needed_src = needed_src & ~mask_src;
      needed_tgt = needed_tgt & ~mask_tgt;
      }
    a_count = a_count+1;
    }
  return a_count;
  }

mycall t_int aux_route_shift_scatter_sloppy(
  tar_gather_scatter &a,
  const ta_index &perm,
  t_bits needed_src,
  t_bits needed_tgt) {
// FAR;
// 2012-09-25  J. Neumann
// var
  t_int a_count;
  t_int s,t;
  t_bits mask_src;
  t_bits mask_tgt;
  ta_index inv_perm;

  invert_perm(perm,inv_perm);
  a_count = 0;
  while (needed_tgt != 0) {
    t = 0;  // -bits;
    s = -1;
    while (true) {
      s = s+1;
      while ( (s < bits) &&
              (((lo_bit << s) & needed_src) == 0) ) {
        s = s+1;
        // t = t+1;
        }
      if (s >= bits) {
        break;
        }

      if (inv_perm[s] < t) {
        break;
        }
      t = inv_perm[s];
      // assert(t != no_index);

      mask_tgt = lo_bit << t;
      mask_src = lo_bit << s;
      a[a_count].mask_tgt = a[a_count].mask_tgt | mask_tgt;
      a[a_count].mask_src = a[a_count].mask_src | mask_src;
      needed_tgt = needed_tgt & ~mask_tgt;
      needed_src = needed_src & ~mask_src;
      }
    a_count = a_count+1;
    }
  return a_count;
  }

mycall t_int aux_route_shift_scatter(
  tar_gather_scatter &a,
  const ta_index &perm,
  t_bits needed_src,
  t_bits needed_tgt) {
// FAR;
// 2012-09-26  J. Neumann
// var
  t_int a_count;
  t_int s,t;
  t_bits mask_src;

  t_bits mask_tgt;
  ta_index inv_perm;

  invert_perm(perm,inv_perm);
  a_count = 0;
  while (needed_tgt != 0) {
    t = 0;  // -bits;
    s = 0;
    while ( (s < bits) &&
            (((lo_bit << s) & needed_src) == 0) ) {
      s = s+1;
      }
    s = s-1;
    while (true) {
      s = s+1;
      if (s >= bits) {
        break;
        }

      if (inv_perm[s] == no_index) {
        break;
        }
      if (inv_perm[s] < t) {
        break;
        }
      t = inv_perm[s];

      mask_tgt = lo_bit << t;
      mask_src = lo_bit << s;
      a[a_count].mask_tgt = a[a_count].mask_tgt | mask_tgt;
      a[a_count].mask_src = a[a_count].mask_src | mask_src;
      needed_tgt = needed_tgt & ~mask_tgt;
      needed_src = needed_src & ~mask_src;
      }
    a_count = a_count+1;
    }
  return a_count;
  }

mycall void decorate_route_gs(tr_imp &imp, t_bool opt) {

  if (opt) {
    imp.description = imp.description+" opt";
    }
  else {
    imp.description = imp.description+" pure";
    }
  }

mycall void frame_route_gs(tr_imp &imp, t_bool opt, tf_aux_route_gs fn, const t_string &des) {
// 2012-09-27  J. Neumann
// var
  ta_index perm;
  tar_gather_scatter a;
  t_int a_count;
  t_int i,j,s,t;
  t_bits mask_src;
  t_bits mask_tgt;
  t_int skip;
  t_int diff;
  t_bool ok;

  init_imp(imp,perm);
  concat_description(imp.description, des);
  decorate_route_gs(imp,opt);
  skip = 0;
  do {
    ok = true;
    imp.a_step[imp.nr_step].nr_pprim = skip;
    for (i = 0; i<=bits-1; ++i) {
      init_gather_scatter(a[i]);
      }
    a_count = (*fn)(a,perm,used_source_bits(perm),used_target_bits(perm));
    for (i = 0; i<=a_count-1; ++i) {
      mask_src = a[i].mask_src;
      mask_tgt = a[i].mask_tgt;

      if (opt) {
        s = nr_trailing_0bits(mask_src);
        t = nr_trailing_0bits(mask_tgt);
        if (mask_src >> s == mask_tgt >> t) {
          imp.a_step[imp.nr_step].nr_pprim = skip;
          skip = skip+1;
          diff = (t-s) & (bits-1);
          mask_src = 0;
          for (j = 0; j<=bits-1; ++j) {
            if ( (perm[j] != no_index) &&
                 (((j-perm[j]) & (bits-1)) == diff) ) {
              mask_src = mask_src | (lo_bit << j);
              perm[j] = no_index;
              }
            }
          add_rol(imp, rol(mask_src,-diff), diff);
          ok = false;
          break;
          }
        }
      add_gs(imp,mask_src,mask_tgt);
      }
    } while (!ok);
  finish_step(imp);
  finish_perm(imp);
  }

mycall void route_gather_scatter(tr_imp &imp, t_bool opt) {

  frame_route_gs(imp,opt,aux_route_gather_scatter, "Gather/scatter");
  }

mycall void route_gather_shift_sloppy(tr_imp &imp, t_bool opt) {

  frame_route_gs(imp,opt,aux_route_gather_shift_sloppy, "Gather/shift sloppy");
  }

mycall void route_gather_shift(tr_imp &imp, t_bool opt) {

  frame_route_gs(imp,opt,aux_route_gather_shift, "Gather/shift");
  }

mycall void route_shift_scatter_sloppy(tr_imp &imp, t_bool opt) {

  frame_route_gs(imp,opt,aux_route_shift_scatter_sloppy, "Shift/scatter sloppy");
  }

mycall void route_shift_scatter(tr_imp &imp, t_bool opt) {

  frame_route_gs(imp,opt,aux_route_shift_scatter, "Shift/scatter");
  }


//////
// Best performance

// var
  tr_performance best_performance;
  tr_imp best_imp;

mycall void use_best_imp(const tr_imp &imp) {
// var
  tr_performance performance;

  imp_performance(imp, performance);
  if (cmp_performance(performance, best_performance) < 0) {
    best_performance = performance;
    best_imp = imp;
    }
  }

mycall void check_best(const tr_imp &imp) {

  if (options.verbose) {
    imp_dump_simple(imp);
    }
  imp_check(imp);
  use_best_imp(imp);
  }


//////
// Test routing variations (BPC / Benes)

typedef
  mycall void (*tf_test_route) (const ta_subword &idx);

mycall void permute_sub(t_int x, ta_subword &idx, tf_test_route fn) {
// var
  t_int i;
  t_int q;

  if (x == ld_bits-1) {
    (*fn)(idx);
    }
  else {
    for (i = x+1; i<=ld_bits-1; ++i) {
      permute_sub(x+1,idx,fn);
      q = idx[x];
      idx[x] = idx[i];
      idx[i] = q;
      }
    permute_sub(x+1,idx,fn);
    for (i = ld_bits-1; i>=x+1; --i) {
      q = idx[x];
      idx[x] = idx[i];
      idx[i] = q;
      }
    }
  }

mycall void permute(tf_test_route fn) {
// var
  ta_subword idx;
  t_int i;

  for (i = 0; i<=ld_bits-1; ++i) {
    idx[i] = i;
    }
  permute_sub(0, idx, fn);
  }

// var
  t_bool is_bpc_perm;
  tr_performance bpc_performance;
  tr_performance benes_performance;

mycall void route_bpc(tr_imp &imp, const ta_subword &tgt, t_bit_index k) {
// 2012-08-31
// var
  ta_index perm;
  ta_subword src,inv_src;
  t_int n,m;
  t_subword i,j,ii,kk;

  init_imp(imp,perm);
  concat_description(imp.description, "BPC permutation ");

  for (i = 0; i<=ld_bits-1; ++i) {
    imp.description = imp.description+num(tgt[i]);
    src[i] = i;
    inv_src[i] = i;
    }
  imp.description = imp.description+'/'+hex(k);

  kk = 0;  // k for generated

  for (i = 0; i<=ld_bits-1; ++i) {  // any order
    n = 1 << i;
    ii = src[i];
    if (tgt[i] == ii) {
      if (((k ^ kk) & n) != 0) {
        // x = bit_index_complement(x,i);
        imp.a_step[imp.nr_step].needed_src_bits = all_bits;
        add_permute(imp, a_bfly_mask[i], n);  // 1 << i;
        }
      }
    else {
      j = inv_src[tgt[i]];
      m = 1 << j;
      if (((k & n) != 0) ^  // boolean xor
          ((kk & m) != 0)) {
        // x = bit_index_swap_complement(x,i,j);
        imp.a_step[imp.nr_step].needed_src_bits = all_bits;
        add_permute(imp, a_bfly_mask[j] & a_bfly_mask[i], m+n);
        if ((kk & n) == 0) {
          kk = kk | m;
          }
        else {
          kk = kk & ~m;
          }
        }
      else {
        // x = bit_index_swap(x,i,j);
        imp.a_step[imp.nr_step].needed_src_bits = all_bits;
        add_permute(imp, a_bfly_mask[j] & ~a_bfly_mask[i], m-n);
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
  finish_perm(imp);
  }

mycall void test_route_bpc(const ta_subword &idx) {
// FAR;
// var
  ta_index perm;
  t_bit_index inv;  // 0..bits-1, i.e. set of [0..ld_bits-1]
  t_int i,j;
  t_bool found;
  t_bit_index v;  // 0..bits-1, i.e. set of [0..ld_bits-1]
  tr_imp imp;
  tr_performance performance;

  init_imp(imp,perm);
  for (inv = 0; inv<=bits-1; ++inv) {
    found = true;
    for (j = 0; j<=bits-1; ++j) {
      if (perm[j] != no_index) {
        v = 0;
        for (i = 0; i<=ld_bits-1; ++i) {
          if (((j ^ inv) & (1 << i)) != 0) {
            v = v | (1 << idx[i]);
            }
          }
        if (perm[j] != v) {
          found = false;
          break;
          }
        }
      }
    if (found) {
      is_bpc_perm = true;
      route_bpc(imp,idx,inv);
      imp_check(imp);
      // imp_dump_simple(imp);
      use_best_imp(imp);
      imp_performance(imp, performance);
      if (cmp_performance(performance, bpc_performance) < 0) {
        bpc_performance = performance;
        }
      }
    }
  }

mycall void test_route_benes(const ta_subword &idx) {
// FAR;
// var
  tr_imp imp;
  tr_performance performance;

  route_benes(imp,idx);
  imp_check(imp);
  // imp_dump_simple(imp);
  use_best_imp(imp);
  imp_performance(imp, performance);
  if (cmp_performance(performance, benes_performance) < 0) {
    benes_performance = performance;
    }
  }


//////
// Permutation decorators

mycall void test_all(t_bool ex) {
// needs pre_imp and post_imp correctly initialized
// var
  tr_imp imp;
  ta_index perm;
  t_bool save_is_bpc_perm;

  save_is_bpc_perm=is_bpc_perm;

  if (options.verbose) {
    writeln0();
    }

  if (is_identity()) {
    init_imp(imp,perm);
    // concat_description(imp.description, "Identity");
    finish_perm(imp);
    check_best(imp);
    }
  else {
    init_performance(bpc_performance);
    init_performance(benes_performance);

    if (options.test_bpc && ex) {
      permute(test_route_bpc);
      if (options.verbose) {
        if (is_bpc_perm) {
          writeln(
            "Altered best BPC permutation: "+dump_performance(bpc_performance));
          }
        else {
          writeln("Altered BPC: Not possible");
          }
        }
      }

    if (options.test_benes && ex) {
      permute(test_route_benes);
      if (options.verbose) {
        writeln("Altered best Benes: "+dump_performance(benes_performance));
        }
      }

    if (options.test_bit_groups) {
      route_bit_groups(imp);
      check_best(imp);
      }

    if (options.test_mul) {
      route_mul1(imp);
      check_best(imp);

      route_mul2(imp);
      check_best(imp);

      route_mul3(imp);
      check_best(imp);
      }

    if (options.allow_bmi) {
      if (options.test_gather_scatter) {
        route_gather_scatter(imp,options.opt_gs);
        check_best(imp);
        }

      if (options.test_gather_shift) {
        route_gather_shift(imp,options.opt_gs);
        check_best(imp);
        }

      if (options.test_gather_shift_sloppy) {
        route_gather_shift_sloppy(imp,options.opt_gs);
        check_best(imp);
        }

      if (options.test_shift_scatter) {
        route_shift_scatter(imp,options.opt_gs);
        check_best(imp);
        }

      if (options.test_shift_scatter_sloppy) {
        route_shift_scatter_sloppy(imp,options.opt_gs);
        check_best(imp);
        }

      if (options.test_sag) {
        route_sag(imp);
        check_best(imp);
        }
      }
    }

  is_bpc_perm=save_is_bpc_perm;
  }

mycall void try_pre_rol(const ta_index &perm, t_bool opt_bswap) {
// var
  t_int i,j,xswap;

  if (opt_bswap) {
    xswap = 7;
    }
  else {
    xswap = 0;
    }
  for (i=bits-1; i>=1; --i) {  // all rotates but 0
    // progress();
    write(num(i)+" \x0d");
    for (j=0; j<=xswap; ++j) {
      init_imp0(pre_imp,perm);
      if ((j & 1)!=0) {
        concat_description(pre_imp.description, "Bswap");
        add_bswap(pre_imp);
        finish_step(pre_imp);
        }
      concat_description(pre_imp.description, "Rol "+num(i));
      add_rol(pre_imp,all_bits,i);
      finish_step(pre_imp);
      if ((j & 2)!=0) {
        concat_description(pre_imp.description, "Bswap");
        add_bswap(pre_imp);
        finish_step(pre_imp);
        }

      init_imp0(post_imp,perm);
      if ((j & 4)!=0) {
        concat_description(post_imp.description, "Bswap");
        add_bswap(post_imp);
        finish_step(post_imp);
        }

      test_all(options.opt_rol_ex);
      }
    }
  }

mycall void try_post_rol(const ta_index &perm, t_bool opt_bswap) {
// var
  t_int i,j,xswap;

  if (opt_bswap) {
    xswap = 7;
    }
  else {
    xswap = 0;
    }
  for (i=bits-1; i>=1; --i) {  // all rotates but 0
    // progress();
    write(num(i)+" \x0d");
    for (j=0; j<=xswap; ++j) {
      init_imp0(pre_imp,perm);
      if ((j & 1)!=0) {
        concat_description(pre_imp.description, "Bswap");
        add_bswap(pre_imp);
        finish_step(pre_imp);
        }

      init_imp0(post_imp,perm);
      if ((j & 2)!=0) {
        concat_description(post_imp.description, "Bswap");
        add_bswap(post_imp);
        finish_step(post_imp);
        }
      concat_description(post_imp.description, "Rol "+num(i));
      add_rol(post_imp,all_bits,i);
      finish_step(post_imp);
      if ((j & 4)!=0) {
        concat_description(post_imp.description, "Bswap");
        add_bswap(post_imp);
        finish_step(post_imp);
        }

      test_all(options.opt_rol_ex);
      }
    }
  }

mycall void try_bswap(const ta_index &perm) {
// var
  t_int j;

  for (j=1; j<=3; ++j) {
    init_imp0(pre_imp,perm);
    if ((j & 1)!=0) {
      concat_description(pre_imp.description, "Bswap");
      add_bswap(pre_imp);
      finish_step(pre_imp);
      }

    init_imp0(post_imp,perm);
    if ((j & 2)!=0) {
      concat_description(post_imp.description, "Bswap");
      add_bswap(post_imp);
      finish_step(post_imp);
      }

    test_all(true);
    }
  }

mycall void my_random_perm(ta_index &perm) {
// var
  t_int i;
  t_int hi;
  t_int r;

  random_perm(perm);
  switch (random_int(4)) {
    case 0: {
      r = random_int(bits);
      for (i = 1; i<=r; ++i) {
        perm[random_int(bits)] = no_index;
        }
      break;
      }
    default: break;
    }
  switch (random_int(4)) {
    case 0: {
      hi = random_int(bits);
      for (i = 0; i<=bits-1; ++i) {
        if (perm[i] >= hi) {
          perm[i] = no_index;
          }
        }
      break;
      }
    default: break;
    }
  switch (random_int(4)) {
    case 0: {
      hi = random_int(bits);
      for (i = hi; i<=bits-1; ++i) {
        perm[i] = no_index;
        }
      break;
      }
    default: break;
    }
  }

mycall void self_test() {
// var
  ta_index perm;
  tr_imp imp;
  t_longint loop;

  for (loop = 1030; loop>=0; --loop) {
    if ((loop & 0x3f) == 0) {
      write(num(loop)+"  \x0d");
      fflush(0);
      }

    my_random_perm(perm);
    init_imp0(pre_imp,perm);
    init_imp0(post_imp,perm);

    route_mul1(imp);
    imp_check(imp);

    route_mul2(imp);
    imp_check(imp);

    route_mul3(imp);
    imp_check(imp);

    route_bit_groups(imp);
    imp_check(imp);

    route_benes(imp,a_stage_bwd);
    imp_check(imp);

    route_benes(imp,a_stage_fwd);
    imp_check(imp);

    route_gather_scatter(imp,true);
    imp_check(imp);

    route_gather_shift_sloppy(imp,true);
    imp_check(imp);

    route_gather_shift(imp,true);
    imp_check(imp);

    route_shift_scatter_sloppy(imp,true);
    imp_check(imp);

    route_shift_scatter(imp,true);
    imp_check(imp);

    route_sag(imp);
    imp_check(imp);
    }
  write("        \x0d");

  for (loop = 260; loop>=0; --loop) {
    if ((loop & 0x0f) == 0) {
      write(num(loop)+"  \x0d");
      fflush(0);
      }

    my_random_perm(perm);
    init_imp0(pre_imp,perm);
    init_imp0(post_imp,perm);

    permute(test_route_bpc);
    permute(test_route_benes);
    }
  write("        \x0d");
  writeln0();
  }

// var
  ta_index perm;


//////
// Startup/parameter routines

mycall void split_opt(const t_string &opt, t_string &name, t_string &value) {
// var
  t_int i,len;

  name = "";
  value = "";

  len = length(opt);

  // search for (line comment sign (#)
  i = 1;
  while ( (i <= len) &&
          (opt[str_ofs+i] != '#') ) {
    i = i + 1;
    }
  if (i <= len) {  // we found the #, now trim
    len = i - 1;
    }

  // trim right
  while ( (len > 0) &&
          (opt[str_ofs+len] >= '\x00') &&  // char might be signed
          (opt[str_ofs+len] <= ' ') ) {
    len = len - 1;
    }
  if (len == 0) {
    return;  // nothing to do
    }

  // trim left
  i = 1;
  while ( (i <= len) &&
          (opt[str_ofs+i] >= '\x00') &&  // char might be signed
          (opt[str_ofs+i] <= ' ') ) {
    i = i + 1;
    }

  // skip lead in character if present
  if ( (i <= len) &&
       ( (opt[str_ofs+i] == '/') ||
         (opt[str_ofs+i] == '-') ) ) {
    i = i + 1;
    }

  // extract name
  while ( (i <= len) &&
          (opt[str_ofs+i] != '=') &&
          (opt[str_ofs+i] != ':') ) {
    name = name + locase(opt[str_ofs+i]);
    i = i + 1;
    }

  // skip delimiter
  i = i + 1;

  // extract value
  while (i <= len) {
    value = value + opt[str_ofs+i];
    i = i + 1;
    }
  }

mycall void setup_default_pascal() {

  options.comment_prefix="// ";
  options.comment_postfix="";
  options.hex_prefix = '$';
  options.hex_postfix = "";
  options.op_assign = ":=";
  options.op_and = "and";
  options.op_or = "or";
  options.op_xor = "xor";
  options.op_shl = "shl";
  options.op_shr = "shr";
  }

mycall void setup_default_c() {

  options.comment_prefix="// ";
  options.comment_postfix="";
  options.hex_prefix = "0x";
  options.hex_postfix = "";
  options.op_assign = '=';
  options.op_and = '&';
  options.op_or = '|';
  options.op_xor = '^';
  options.op_shl = "<<";
  options.op_shr = ">>";
  }

mycall void recalc_cost_permute() {

  options.cost_bit_permute_step =
    options.cost_shift*2+options.cost_xor*3+options.cost_and;
  options.cost_bit_permute_step_simple =
    options.cost_shift*2+options.cost_and*2+options.cost_or;
  }

mycall void recalc_cost_opt() {

  options.cost_rotate = options.cost_rotate_shift;
  options.cost_shift = options.cost_rotate_shift;
  options.cost_and = options.cost_bool;
  options.cost_or = options.cost_bool;
  options.cost_xor = options.cost_bool;
  options.cost_scatter = options.cost_gs;
  options.cost_gather = options.cost_gs;

  recalc_cost_permute();
  }

mycall void setup_defaults() {

  options.dump_input = true;
  options.dump_inverse = true;
  options.verbose = true;
  options.brief = true;

  // setup_default_pascal();
  setup_default_c();

  options.op_pstep = "bit_permute_step";
  options.op_pstep_simple = "bit_permute_step_simple";
  options.op_rol = "rol";
  options.op_gather = "pext";
  options.op_scatter = "pdep";
  options.op_bswap = "bswap";

  options.in_origin = 0;  // input index origin
  options.in_base = 10;  // input number base
  options.in_indexes_are_target = false;  // /in_indexes=source or target

  options.cost_rotate_shift = 1;
  options.cost_bool = 1;
  options.cost_bswap = 1;
  options.cost_mul = 4;
  options.cost_gs = 3;
  options.cost_mask = 0;

  recalc_cost_opt();

  options.bonus_bit_permute_step = 1;  // implicitely parallel
  options.bonus_bit_permute_step_simple = 1;  // implicitely parallel
  options.bonus_gs = 3;  // 2 parallel gs
  options.bonus_mask_rol = 2;  // 2 parallel mask_rol/mask_shift ops
  options.bonus_gs_rol = 1;  // parallel mask_rol/mask_shift and gs ops

  options.allow_bswap = true;
  options.allow_bmi = false;
  options.test_bpc = true;
  options.test_bfly = true;
  options.test_ibfly = true;
  options.test_benes = true;
  options.test_bit_groups = true;
  options.test_mul = true;
  options.test_gather_scatter = true;
  options.test_gather_shift = true;
  options.test_gather_shift_sloppy = true;
  options.test_shift_scatter = true;
  options.test_shift_scatter_sloppy = true;
  options.test_sag = true;
  options.opt_gs = true;
  options.opt_rol = true;
  options.opt_rol_ex = false;
  options.opt_bswap = true;

  options.self_test = false;
  }

mycall void setup(const t_string &opt);  // forward

mycall void setup_file(const t_string &n) {
// var
  FILE* f;
  t_string s;
  t_int i;
  t_string n1;
  t_char buf[256];

  n1 = "";
  for (i = 1; i<=length(n); ++i) {
    if (n[str_ofs+i]!='"') {
      n1 = n1 + n[str_ofs+i];
      }
    }

  f = fopen(n1.c_str(),"r");
  while (fgets(buf,sizeof(buf)-1,f)) {
    s = buf;  // char array => string
    while ( (s != "") &&
            (s[str_ofs+length(s)] == '\x0a') ) {
      s = copy(s,1,length(s)-1);
      }
    setup(s);
    }
  fclose(f);
  }

mycall void setup(const t_string &opt) {
// var
  t_string n,v;

  if (opt == "") {
    // ignore
    }
  else if (opt[str_ofs+1] == '@') {
    setup_file(copy(opt,2,length(opt)-1));
    }
  else {
    split_opt(opt,n,v);

    if (n == "") {
      // ignore
      }

    else if (n == "dump_input") {
      options.dump_input = s2i(v)!=0;
      }
    else if (n == "dump_inverse") {
      options.dump_inverse = s2i(v)!=0;
      }
    else if (n == "verbose") {
      options.verbose = s2i(v)!=0;
      }
    else if (n == "brief") {
      options.brief = s2i(v)!=0;
      }

    else if (n == "output_pas") {
      setup_default_pascal();
      }
    else if (n == "output_c") {
      setup_default_c();
      }

    else if (n == "comment_prefix") {
      if (v == "") {
        options.comment_prefix = "";
        }
      else {
        options.comment_prefix = v+' ';
        }
      }
    else if (n == "comment_postfix") {
      if (v == "") {
        options.comment_postfix = "";
        }
      else {
        options.comment_postfix = ' '+v;
        }
      }
    else if (n == "hex_prefix") {
      options.hex_prefix = v;
      }
    else if (n == "hex_postfix") {
      options.hex_postfix = v;
      }
    else if (n == "op_assign") {
      options.op_assign = v;
      }
    else if (n == "op_and") {
      options.op_and = v;
      }
    else if (n == "op_or") {
      options.op_or = v;
      }
    else if (n == "op_xor") {
      options.op_xor = v;
      }
    else if (n == "op_shl") {
      options.op_shl = v;
      }
    else if (n == "op_shr") {
      options.op_shr = v;
      }

    else if (n == "in_origin") {
      options.in_origin = s2i(v);
      }
    else if (n == "in_base") {
      options.in_base = s2i(v);
      }
    else if (n == "in_indexes") {
      v = lostr(v);
      if (v == "source") {
        options.in_indexes_are_target = false;
        }
      else if (v == "target") {
        options.in_indexes_are_target = true;
        }
      else {
        writeln("ERROR: Unknown option to /in_indexes "+v);
        error_abort();
        }
      }

    else if (n == "op_pstep") {
      options.op_pstep = v;
      }
    else if (n == "op_pstep_simple") {
      options.op_pstep_simple = v;
      }
    else if (n == "op_rol") {
      options.op_rol = v;
      }
    else if (n == "op_gather") {
      options.op_gather = v;
      }
    else if (n == "op_scatter") {
      options.op_scatter = v;
      }
    else if (n == "op_bswap") {
      options.op_bswap = v;
      }

    else if (n == "cost_rotate_shift") {
      options.cost_rotate_shift = s2i(v);
      recalc_cost_opt();
      }
    else if (n == "cost_bool") {
      options.cost_bool = s2i(v);
      recalc_cost_opt();
      }
    else if (n == "cost_bswap") {
      options.cost_bswap = s2i(v);
      recalc_cost_opt();
      }
    else if (n == "cost_mul") {
      options.cost_mul = s2i(v);
      recalc_cost_opt();
      }
    else if (n == "cost_gs") {
      options.cost_gs = s2i(v);
      recalc_cost_opt();
      }
    else if (n == "cost_mask") {
      options.cost_mask = s2i(v);
      }

    else if (n == "cost_rotate") {
      options.cost_rotate = s2i(v);
      recalc_cost_permute();
      }
    else if (n == "cost_shift") {
      options.cost_shift = s2i(v);
      recalc_cost_permute();
      }
    else if (n == "cost_and") {
      options.cost_and = s2i(v);
      recalc_cost_permute();
      }
    else if (n == "cost_or") {
      options.cost_or = s2i(v);
      recalc_cost_permute();
      }
    else if (n == "cost_xor") {
      options.cost_xor = s2i(v);
      recalc_cost_permute();
      }
    else if (n == "cost_scatter") {
      options.cost_scatter = s2i(v);
      }
    else if (n == "cost_gather") {
      options.cost_gather = s2i(v);
      }

    else if (n == "cost_bit_permute_step") {
      options.cost_bit_permute_step = s2i(v);
      }
    else if (n == "cost_bit_permute_step_simple") {
      options.cost_bit_permute_step_simple = s2i(v);
      }

    else if (n == "bonus_bit_permute_step") {
      options.bonus_bit_permute_step = s2i(v);
      }
    else if (n == "bonus_bit_permute_step_simple") {
      options.bonus_bit_permute_step_simple = s2i(v);
      }

    else if (n == "bonus_gs") {
      options.bonus_gs = s2i(v);
      }
    else if (n == "bonus_mask_rol") {
      options.bonus_mask_rol = s2i(v);
      }
    else if (n == "bonus_gs_rol") {
      options.bonus_gs_rol = s2i(v);
      }

    else if (n == "allow_bswap") {
      options.allow_bswap = s2i(v)!=0;
      }
    else if (n == "allow_bmi") {
      options.allow_bmi = s2i(v)!=0;
      }
    else if (n == "test_bpc") {
      options.test_bpc = s2i(v)!=0;
      }
    else if (n == "test_bfly") {
      options.test_bfly = s2i(v)!=0;
      }
    else if (n == "test_ibfly") {
      options.test_ibfly = s2i(v)!=0;
      }
    else if (n == "test_benes") {
      options.test_benes = s2i(v)!=0;
      }
    else if (n == "test_bit_groups") {
      options.test_bit_groups = s2i(v)!=0;
      }
    else if (n == "test_mul") {
      options.test_mul = s2i(v)!=0;
      }
    else if (n == "test_gather_scatter") {
      options.test_gather_scatter = s2i(v)!=0;
      }
    else if (n == "test_gather_shift") {
      options.test_gather_shift = s2i(v)!=0;
      }
    else if (n == "test_gather_shift_sloppy") {
      options.test_gather_shift_sloppy = s2i(v)!=0;
      }
    else if (n == "test_shift_scatter") {
      options.test_shift_scatter = s2i(v)!=0;
      }
    else if (n == "test_shift_scatter_sloppy") {
      options.test_shift_scatter_sloppy = s2i(v)!=0;
      }
    else if (n == "test_sag") {
      options.test_sag = s2i(v)!=0;
      }
    else if (n == "opt_gs") {
      options.opt_gs = s2i(v)!=0;
      }
    else if (n == "opt_rol") {
      options.opt_rol = s2i(v)!=0;
      }
    else if (n == "opt_rol_ex") {
      options.opt_rol_ex = s2i(v)!=0;
      }
    else if (n == "opt_bswap") {
      options.opt_bswap = s2i(v)!=0;
      }

    else if (n == "self_test") {
      options.self_test = s2i(v)!=0;
      }

    else {
      writeln("ERROR: Unknown option "+n);
      error_abort();
      }
    }
  }

mycall void get_perm(const t_string &opt, ta_index &perm) {
// var
  t_int i;
  t_int idx;
  t_bits got;
  t_char c;
  t_int my_base;
  t_int val;
  t_int letter_ofs;
  t_int digit;
  t_string s;

  for (i = 0; i<=bits-1; ++i) {
    perm[i] = no_index;
    }

  my_base = options.in_base;
  got = 0;
  idx = -1;
  i = 1;
  while (i <= length(opt)) {
    c = opt[str_ofs+i];
    val = -2;
    switch (c) {

      case '\f':
      case '\n':
      case '\r':
      case '\t':
      case '\v':
      case ' ':
      case ',': {
        i = i+1;
        continue;
        break;
        }

      case '*': {
        val = -1;
        break;
        }

      case '$': {
        my_base = 16;
        i = i+1;
        continue;
        break;
        }

      case '@': {
        i = i+1;  // skip @
        s = "";
        while ( (i <= length(opt)) &&
                (opt[str_ofs+i] != ' ') ) {
          if (opt[str_ofs+i] != '"') {
            s = s+opt[str_ofs+i];
            i = i+1;
            }
          else {
            while ( (i <= length(opt)) &&
                    (opt[str_ofs+i] != '"') ) {
              s = s+opt[str_ofs+i];
              i = i+1;
              }
            i = i+1;
            }
          }
        setup_file(s);
        my_base = options.in_base;  // in_base might have changed
        continue;
        break;
        }

      case '/':
      case '-': {
        s = "";
        while ( (i <= length(opt)) &&
                (opt[str_ofs+i] != ' ') ) {
          s = s+opt[str_ofs+i];
          i = i+1;
          }
        setup(s);
        my_base = options.in_base;  // in_base might have changed
        continue;
        break;
        }

      case '#': {
        goto BREAK_1;  // comment until end of line
        break;
        }

      default: {
        if ( (('0'<=c) && (c<='9')) ||
             (('a'<=c) && (c<='z')) ||
             (('A'<=c) && (c<='Z')) ) {
          if (my_base == 26) {
            letter_ofs = 0;
            }
          else {
            letter_ofs = 10;
            }
          val = 0;
          digit = 0;  // Compiler, shut up!
          while (i <= length(opt)) {
            c = opt[str_ofs+i];
            if (false) {}
            else if (('0'<=c) && (c<='9')) {
              digit = (t_int)(c)-(t_int)('0');
              }
            else if (('a'<=c) && (c<='z')) {
              digit = (t_int)(c)-(t_int)('a')+letter_ofs;
              }
            else if (('A'<=c) && (c<='Z')) {
              digit = (t_int)(c)-(t_int)('A')+letter_ofs;
              }
            else {
              goto BREAK_2;
              }
            if ( (my_base != 26) &&
                 (digit >= my_base) ) {
              writeln("ERROR: Invalid digit");
              error_abort();
              }
            val = val*my_base + digit;
            i = i+1;
            }
        BREAK_2:
          i = i-1;
          }
        break;
        }

      }
    if (val == -2) {
      writeln("ERROR: Illegal character");
      i = i+1;
      continue;
      }
    idx = idx+1;
    if (val == -1) {
      // ignore
      }
    else if ( (val < options.in_origin) ||
              (val-options.in_origin > bits-1) ) {
      writeln("ERROR: Out of range");
      error_abort();
      }
    else {
      val = val - options.in_origin;
      if (((lo_bit << val) & got) != 0) {
        writeln("ERROR: Dupe: "+num(val));
        error_abort();
        }
      got = got | (lo_bit << val);
      }
    perm[idx] = val;
    my_base = options.in_base;
    i = i+1;
    }

BREAK_1:
  if (idx < 0) {
    // random_perm(perm);
    // val = random_int(bits);
    // for (i = 1; i<=val; ++i) {
    //   perm[random_int(bits)] = no_index;
    //   }
    my_random_perm(perm);
    }
  }

mycall void startup(int argc, char* argv[]) {
// var
  t_string opt;
  t_int i;
  ta_index inv_perm;

  setup_defaults();

  setup_file("calcperm.ini");

  opt = "";
  for (i = 1; i<=argc-1; ++i) {
    opt = opt+" "+argv[i];
    }

  get_perm(opt, perm);

  if (options.in_indexes_are_target) {
    invert_perm(perm,inv_perm);
    for (i=0; i<=bits-1; ++i) {
      perm[i] = inv_perm[i];
      }
    }

  init_imp0(pre_imp,perm);
  init_imp0(post_imp,perm);
  }

// var
  ta_index inv_perm;
  tr_imp imp;
  t_bool is_bfly_perm;
  t_bool is_ibfly_perm;


//////
// Main program

int main(int argc, char* argv[]) {

  if (!init_general()) {
    return 1;
    }

  printf("Permutation code generator (in C++, bits=%i)...\n", bits);
  if (sizeof(t_bits)*8 != bits) {
    printf("Sizes wrong! sizeof(t_bits)=%i\n", (int)sizeof(t_bits));
    return 1;
    }

  srand((unsigned int)(time(NULL)));
  my_randseed = (t_bits)(rand());

  startup(argc, argv);

  if (options.self_test) {
    self_test();
    writeln("OK");
    }
  else {
    if (options.dump_input) {
      write("Permutation vector: ");
      dump_perm(perm);
      }

    if (options.dump_inverse) {
      write("Inverse vector: ");
      invert_perm(perm,inv_perm);
      dump_perm(inv_perm);
      }

    init_performance(best_performance);
    init_performance(bpc_performance);
    init_performance(benes_performance);
    is_bpc_perm = false;
    is_bfly_perm = false;
    is_ibfly_perm = false;

    if (options.verbose) {
      writeln0();
      }

    if (is_identity()) {
      is_bpc_perm = true;
      is_bfly_perm = true;
      is_ibfly_perm = true;
      init_imp(imp,perm);
      concat_description(imp.description, "Identity");
      finish_perm(imp);
      check_best(imp);
      }
    else {

      if (options.test_bpc) {
        permute(test_route_bpc);
        if (options.verbose) {
          if (is_bpc_perm) {
            writeln("Best BPC permutation: "+dump_performance(bpc_performance));
            }
          }
        }

      if (options.test_bfly) {
        route_benes(imp,a_stage_bwd);
        if (imp.is_bfly) {
          check_best(imp);
          is_bfly_perm = true;
          }
        }

      if (options.test_ibfly) {
        route_benes(imp,a_stage_fwd);
        if (imp.is_bfly) {
          check_best(imp);
          is_ibfly_perm = true;
          }
        }

      if (options.test_benes) {
        permute(test_route_benes);
        if (options.verbose) {
          writeln("Best Benes: "+dump_performance(benes_performance));
          }
        }

      if (options.test_bit_groups) {
        route_bit_groups(imp);
        check_best(imp);
        }

      if (options.test_mul) {
        route_mul1(imp);
        check_best(imp);

        route_mul2(imp);
        check_best(imp);

        route_mul3(imp);
        check_best(imp);
        }

      if (options.allow_bmi) {
        if (options.test_gather_scatter) {
          route_gather_scatter(imp,options.opt_gs);
          check_best(imp);
          }

        if (options.test_gather_shift) {
          route_gather_shift(imp,options.opt_gs);
          check_best(imp);
          }

        if (options.test_gather_shift_sloppy) {
          route_gather_shift_sloppy(imp,options.opt_gs);
          check_best(imp);
          }

        if (options.test_shift_scatter) {
          route_shift_scatter(imp,options.opt_gs);
          check_best(imp);
          }

        if (options.test_shift_scatter_sloppy) {
          route_shift_scatter_sloppy(imp,options.opt_gs);
          check_best(imp);
          }

        if (options.test_sag) {
          route_sag(imp);
          check_best(imp);
          }
        }

      if (options.opt_rol) {
        try_pre_rol(perm, options.opt_bswap && options.allow_bswap);
        try_post_rol(perm, options.opt_bswap && options.allow_bswap);
        write(" \x0d");  // get rid of progress display
        }
      if (options.opt_bswap && options.allow_bswap) {
        try_bswap(perm);
        }
      }

    if (options.brief || options.verbose) {
      writeln0();
      if (options.test_bpc) {
        writeln("BPC permutation: "+q2s(is_bpc_perm));
        }
      if (options.test_bfly) {
        writeln("Butterfly: "+q2s(is_bfly_perm));
        }
      if (options.test_ibfly) {
        writeln("Inverse Butterfly: "+q2s(is_ibfly_perm));
        }

      writeln0();
      writeln("=== Best method ===");
      }

    if (best_performance.cost == maxint) {
      writeln("NOT ROUTABLE WITH ALLOWED METHODS");
      }
    else {
      imp_dump(best_imp);
      }
    }
  }

// eof.
