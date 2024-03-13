/* This file was automatically generated by CasADi 3.6.3.
 *  It consists of: 
 *   1) content generated by CasADi runtime: not copyrighted
 *   2) template code copied from CasADi source: permissively licensed (MIT-0)
 *   3) user code: owned by the user
 *
 */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) tricycle_LF_cost_y_fun_ ## ID
#endif

#include <math.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int int
#endif

/* Add prefix to internal symbols */
#define casadi_copy CASADI_PREFIX(copy)
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)
#define casadi_s3 CASADI_PREFIX(s3)
#define casadi_s4 CASADI_PREFIX(s4)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

void casadi_copy(const casadi_real* x, casadi_int n, casadi_real* y) {
  casadi_int i;
  if (y) {
    if (x) {
      for (i=0; i<n; ++i) *y++ = *x++;
    } else {
      for (i=0; i<n; ++i) *y++ = 0.;
    }
  }
}

static const casadi_int casadi_s0[12] = {8, 1, 0, 8, 0, 1, 2, 3, 4, 5, 6, 7};
static const casadi_int casadi_s1[6] = {2, 1, 0, 2, 0, 1};
static const casadi_int casadi_s2[3] = {0, 0, 0};
static const casadi_int casadi_s3[20] = {16, 1, 0, 16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
static const casadi_int casadi_s4[11] = {7, 1, 0, 7, 0, 1, 2, 3, 4, 5, 6};

/* tricycle_LF_cost_y_fun:(i0[8],i1[2],i2[],i3[],i4[16])->(o0[7]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real *rr, *ss;
  casadi_real w0, w1, w2, w3, w4, w5, w6, w7, *w8=w+8, *w9=w+16;
  /* #0: @0 = input[0][0] */
  w0 = arg[0] ? arg[0][0] : 0;
  /* #1: @1 = input[0][1] */
  w1 = arg[0] ? arg[0][1] : 0;
  /* #2: @2 = input[0][2] */
  w2 = arg[0] ? arg[0][2] : 0;
  /* #3: @3 = input[0][3] */
  w3 = arg[0] ? arg[0][3] : 0;
  /* #4: @4 = input[0][4] */
  w4 = arg[0] ? arg[0][4] : 0;
  /* #5: @5 = input[0][5] */
  w5 = arg[0] ? arg[0][5] : 0;
  /* #6: @6 = input[0][6] */
  w6 = arg[0] ? arg[0][6] : 0;
  /* #7: @7 = input[0][7] */
  w7 = arg[0] ? arg[0][7] : 0;
  /* #8: @8 = vertcat(@0, @1, @2, @3, @4, @5, @6, @7) */
  rr=w8;
  *rr++ = w0;
  *rr++ = w1;
  *rr++ = w2;
  *rr++ = w3;
  *rr++ = w4;
  *rr++ = w5;
  *rr++ = w6;
  *rr++ = w7;
  /* #9: @9 = @8[3:8] */
  for (rr=w9, ss=w8+3; ss!=w8+8; ss+=1) *rr++ = *ss;
  /* #10: output[0][0] = @9 */
  casadi_copy(w9, 5, res[0]);
  /* #11: @0 = input[1][0] */
  w0 = arg[1] ? arg[1][0] : 0;
  /* #12: output[0][1] = @0 */
  if (res[0]) res[0][5] = w0;
  /* #13: @0 = input[1][1] */
  w0 = arg[1] ? arg[1][1] : 0;
  /* #14: output[0][2] = @0 */
  if (res[0]) res[0][6] = w0;
  return 0;
}

CASADI_SYMBOL_EXPORT int tricycle_LF_cost_y_fun(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int tricycle_LF_cost_y_fun_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int tricycle_LF_cost_y_fun_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tricycle_LF_cost_y_fun_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int tricycle_LF_cost_y_fun_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tricycle_LF_cost_y_fun_release(int mem) {
}

CASADI_SYMBOL_EXPORT void tricycle_LF_cost_y_fun_incref(void) {
}

CASADI_SYMBOL_EXPORT void tricycle_LF_cost_y_fun_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int tricycle_LF_cost_y_fun_n_in(void) { return 5;}

CASADI_SYMBOL_EXPORT casadi_int tricycle_LF_cost_y_fun_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real tricycle_LF_cost_y_fun_default_in(casadi_int i) {
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tricycle_LF_cost_y_fun_name_in(casadi_int i) {
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    case 4: return "i4";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tricycle_LF_cost_y_fun_name_out(casadi_int i) {
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tricycle_LF_cost_y_fun_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    case 2: return casadi_s2;
    case 3: return casadi_s2;
    case 4: return casadi_s3;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tricycle_LF_cost_y_fun_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s4;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int tricycle_LF_cost_y_fun_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 13;
  if (sz_res) *sz_res = 2;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 21;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
