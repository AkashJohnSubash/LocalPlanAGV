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
  #define CASADI_PREFIX(ID) tricycle_frenet_impl_dae_fun_ ## ID
#endif

#include <math.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int int
#endif

/* Add prefix to internal symbols */
#define casadi_c0 CASADI_PREFIX(c0)
#define casadi_c1 CASADI_PREFIX(c1)
#define casadi_clear CASADI_PREFIX(clear)
#define casadi_clear_casadi_int CASADI_PREFIX(clear_casadi_int)
#define casadi_copy CASADI_PREFIX(copy)
#define casadi_de_boor CASADI_PREFIX(de_boor)
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_f1 CASADI_PREFIX(f1)
#define casadi_fill CASADI_PREFIX(fill)
#define casadi_fill_casadi_int CASADI_PREFIX(fill_casadi_int)
#define casadi_low CASADI_PREFIX(low)
#define casadi_nd_boor_eval CASADI_PREFIX(nd_boor_eval)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)
#define casadi_s3 CASADI_PREFIX(s3)
#define casadi_s4 CASADI_PREFIX(s4)
#define casadi_s5 CASADI_PREFIX(s5)
#define casadi_s6 CASADI_PREFIX(s6)
#define casadi_s7 CASADI_PREFIX(s7)

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

void casadi_de_boor(casadi_real x, const casadi_real* knots, casadi_int n_knots, casadi_int degree, casadi_real* boor) {
  casadi_int d, i;
  for (d=1;d<degree+1;++d) {
    for (i=0;i<n_knots-d-1;++i) {
      casadi_real b, bottom;
      b = 0;
      bottom = knots[i + d] - knots[i];
      if (bottom) b = (x - knots[i]) * boor[i] / bottom;
      bottom = knots[i + d + 1] - knots[i + 1];
      if (bottom) b += (knots[i + d + 1] - x) * boor[i + 1] / bottom;
      boor[i] = b;
    }
  }
}

void casadi_fill(casadi_real* x, casadi_int n, casadi_real alpha) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = alpha;
  }
}

void casadi_fill_casadi_int(casadi_int* x, casadi_int n, casadi_int alpha) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = alpha;
  }
}

void casadi_clear(casadi_real* x, casadi_int n) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = 0;
  }
}

void casadi_clear_casadi_int(casadi_int* x, casadi_int n) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = 0;
  }
}

casadi_int casadi_low(casadi_real x, const casadi_real* grid, casadi_int ng, casadi_int lookup_mode) {
  switch (lookup_mode) {
    case 1:
      {
        casadi_real g0, dg;
        casadi_int ret;
        g0 = grid[0];
        dg = grid[ng-1]-g0;
        ret = (casadi_int) ((x-g0)*(ng-1)/dg);
        if (ret<0) ret=0;
        if (ret>ng-2) ret=ng-2;
        return ret;
      }
    case 2:
      {
        casadi_int start, stop, pivot;
        if (ng<2 || x<grid[1]) return 0;
        if (x>grid[ng-1]) return ng-2;
        start = 0;
        stop  = ng-1;
        while (1) {
          pivot = (stop+start)/2;
          if (x < grid[pivot]) {
            if (pivot==stop) return pivot;
            stop = pivot;
          } else {
            if (pivot==start) return pivot;
            start = pivot;
          }
        }
      }
    default:
      {
        casadi_int i;
        for (i=0; i<ng-2; ++i) {
          if (x < grid[i+1]) break;
        }
        return i;
      }
  }
}

void casadi_nd_boor_eval(casadi_real* ret, casadi_int n_dims, const casadi_real* all_knots, const casadi_int* offset, const casadi_int* all_degree, const casadi_int* strides, const casadi_real* c, casadi_int m, const casadi_real* all_x, const casadi_int* lookup_mode, casadi_int* iw, casadi_real* w) {
  casadi_int n_iter, k, i, pivot;
  casadi_int *boor_offset, *starts, *index, *coeff_offset;
  casadi_real *cumprod, *all_boor;
  boor_offset = iw; iw+=n_dims+1;
  starts = iw; iw+=n_dims;
  index = iw; iw+=n_dims;
  coeff_offset = iw;
  cumprod = w; w+= n_dims+1;
  all_boor = w;
  boor_offset[0] = 0;
  cumprod[n_dims] = 1;
  coeff_offset[n_dims] = 0;
  n_iter = 1;
  for (k=0;k<n_dims;++k) {
    casadi_real *boor;
    const casadi_real* knots;
    casadi_real x;
    casadi_int degree, n_knots, n_b, L, start;
    boor = all_boor+boor_offset[k];
    degree = all_degree[k];
    knots = all_knots + offset[k];
    n_knots = offset[k+1]-offset[k];
    n_b = n_knots-degree-1;
    x = all_x[k];
    L = casadi_low(x, knots+degree, n_knots-2*degree, lookup_mode[k]);
    start = L;
    if (start>n_b-degree-1) start = n_b-degree-1;
    starts[k] = start;
    casadi_clear(boor, 2*degree+1);
    if (x>=knots[0] && x<=knots[n_knots-1]) {
      if (x==knots[1]) {
        casadi_fill(boor, degree+1, 1.0);
      } else if (x==knots[n_knots-1]) {
        boor[degree] = 1;
      } else if (knots[L+degree]==x) {
        boor[degree-1] = 1;
      } else {
        boor[degree] = 1;
      }
    }
    casadi_de_boor(x, knots+start, 2*degree+2, degree, boor);
    boor+= degree+1;
    n_iter*= degree+1;
    boor_offset[k+1] = boor_offset[k] + degree+1;
  }
  casadi_clear_casadi_int(index, n_dims);
  for (pivot=n_dims-1;pivot>=0;--pivot) {
    cumprod[pivot] = (*(all_boor+boor_offset[pivot]))*cumprod[pivot+1];
    coeff_offset[pivot] = starts[pivot]*strides[pivot]+coeff_offset[pivot+1];
  }
  for (k=0;k<n_iter;++k) {
    casadi_int pivot = 0;
    for (i=0;i<m;++i) ret[i] += c[coeff_offset[0]+i]*cumprod[0];
    index[0]++;
    {
      while (index[pivot]==boor_offset[pivot+1]-boor_offset[pivot]) {
        index[pivot] = 0;
        if (pivot==n_dims-1) break;
        index[++pivot]++;
      }
      while (pivot>0) {
        cumprod[pivot] = (*(all_boor+boor_offset[pivot]+index[pivot]))*cumprod[pivot+1];
        coeff_offset[pivot] = (starts[pivot]+index[pivot])*strides[pivot]+coeff_offset[pivot+1];
        pivot--;
      }
    }
    cumprod[0] = (*(all_boor+index[0]))*cumprod[1];
    coeff_offset[0] = (starts[0]+index[0])*m+coeff_offset[1];
  }
}

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

static const casadi_int casadi_s0[2] = {0, 227};
static const casadi_int casadi_s1[1] = {3};
static const casadi_int casadi_s2[1] = {1};
static const casadi_int casadi_s3[1] = {2};
static const casadi_int casadi_s4[12] = {8, 1, 0, 8, 0, 1, 2, 3, 4, 5, 6, 7};
static const casadi_int casadi_s5[6] = {2, 1, 0, 2, 0, 1};
static const casadi_int casadi_s6[3] = {0, 0, 0};
static const casadi_int casadi_s7[20] = {16, 1, 0, 16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

static const casadi_real casadi_c0[227] = {0., 0., 0., 0., 1.2500000000000000e+00, 1.8750000000000000e+00, 2.5000000000000000e+00, 3.1250000000000000e+00, 3.7500000000000000e+00, 4.3750000000000000e+00, 5., 5.6250000000000000e+00, 6.2500000000000000e+00, 6.8750000000000000e+00, 7.5000000000000000e+00, 8.1250000000000000e+00, 8.7500000000000000e+00, 9.3750000000000000e+00, 10., 1.0245338000000000e+01, 1.0490677000000000e+01, 1.0736015000000000e+01, 1.0981353000000000e+01, 1.1226692000000000e+01, 1.1472030000000000e+01, 1.1717369000000000e+01, 1.1962707000000000e+01, 1.2208045000000000e+01, 1.2453384000000000e+01, 1.2698722000000000e+01, 1.2944060000000000e+01, 1.3189399000000000e+01, 1.3434737000000000e+01, 1.3680076000000000e+01, 1.3925414000000000e+01, 1.4170752000000000e+01, 1.4416091000000000e+01, 1.4661429000000000e+01, 1.4906767000000000e+01, 1.5152106000000000e+01, 1.5397444000000000e+01, 1.5642783000000000e+01, 1.5888121000000000e+01, 1.6133458999999998e+01, 1.6378798000000000e+01, 1.6624136000000000e+01, 1.6869474000000000e+01, 1.7114813000000002e+01, 1.7360150999999998e+01, 1.7605490000000000e+01, 1.7850828000000000e+01, 1.8096166000000000e+01, 1.8341505000000002e+01, 1.8586842999999998e+01, 1.8832180999999999e+01, 1.9077520000000000e+01, 1.9322858000000000e+01, 1.9568196000000000e+01, 1.9813535000000002e+01, 2.0058872999999998e+01, 2.0304212000000000e+01, 2.0549550000000000e+01, 2.0794888000000000e+01, 2.1040227000000002e+01, 2.1285564999999998e+01, 2.1530902999999999e+01, 2.1776242000000000e+01, 2.2021580000000000e+01, 2.2266919000000001e+01, 2.2512257000000002e+01, 2.2757594999999998e+01, 2.3002934000000000e+01, 2.3248272000000000e+01, 2.3493610000000000e+01, 2.3738949000000002e+01, 2.3984286999999998e+01, 2.4229626000000000e+01, 2.4474964000000000e+01, 2.4720302000000000e+01, 2.4965641000000002e+01, 2.5210978999999998e+01, 2.5456316999999999e+01, 2.5701656000000000e+01, 2.6326656000000000e+01, 2.6951656000000000e+01, 2.7576656000000000e+01, 2.8201656000000000e+01, 2.8826656000000000e+01, 2.9451656000000000e+01, 3.0076656000000000e+01, 3.0701656000000000e+01, 3.1326656000000000e+01, 3.1951656000000000e+01, 3.2576656000000000e+01, 3.3201656000000000e+01, 3.3826656000000000e+01, 3.4451656000000000e+01, 3.5076656000000000e+01, 3.5701656000000000e+01, 3.5946993999999997e+01, 3.6192332999999998e+01, 3.6437671000000002e+01, 3.6683008999999998e+01, 3.6928348000000000e+01, 3.7173685999999996e+01, 3.7419024000000000e+01, 3.7664363000000002e+01, 3.7909700999999998e+01, 3.8155040000000000e+01, 3.8400378000000003e+01, 3.8645716000000000e+01, 3.8891055000000001e+01, 3.9136392999999998e+01, 3.9381731000000002e+01, 3.9627070000000003e+01, 3.9872408000000000e+01, 4.0117745999999997e+01, 4.0363084999999998e+01, 4.0608423000000002e+01, 4.0853762000000003e+01, 4.1099100000000000e+01, 4.1344437999999997e+01, 4.1589776999999998e+01, 4.1835115000000002e+01, 4.2080452999999999e+01, 4.2325792000000000e+01, 4.2571129999999997e+01, 4.2816468999999998e+01, 4.3061807000000002e+01, 4.3307144999999998e+01, 4.3552484000000000e+01, 4.4177484000000000e+01, 4.4802484000000000e+01, 4.5427484000000000e+01, 4.6052484000000000e+01, 4.6677484000000000e+01, 4.7302484000000000e+01, 4.7927484000000000e+01, 4.8552484000000000e+01, 4.9177484000000000e+01, 4.9802484000000000e+01, 5.0427484000000000e+01, 5.1052484000000000e+01, 5.1677484000000000e+01, 5.2302484000000000e+01, 5.2927484000000000e+01, 5.3552484000000000e+01, 5.4177484000000000e+01, 5.4802484000000000e+01, 5.5427484000000000e+01, 5.6052484000000000e+01, 5.6677484000000000e+01, 5.7302484000000000e+01, 5.7927484000000000e+01, 5.8552484000000000e+01, 5.9177484000000000e+01, 5.9802484000000000e+01, 6.0427484000000000e+01, 6.1052484000000000e+01, 6.1677484000000000e+01, 6.2302484000000000e+01, 6.2927484000000000e+01, 6.3552484000000000e+01, 6.4177484000000007e+01, 6.4802484000000007e+01, 6.5427484000000007e+01, 6.6052484000000007e+01, 6.6677484000000007e+01, 6.7302484000000007e+01, 6.7927484000000007e+01, 6.8552484000000007e+01, 6.9177484000000007e+01, 6.9802484000000007e+01, 7.0427484000000007e+01, 7.1052484000000007e+01, 7.1677484000000007e+01, 7.2302484000000007e+01, 7.2927484000000007e+01, 7.3552484000000007e+01, 7.4177484000000007e+01, 7.4802484000000007e+01, 7.5427484000000007e+01, 7.6052484000000007e+01, 7.6677484000000007e+01, 7.7302484000000007e+01, 7.7927484000000007e+01, 7.8552484000000007e+01, 7.9177484000000007e+01, 7.9802484000000007e+01, 8.0427484000000007e+01, 8.1052484000000007e+01, 8.1677484000000007e+01, 8.2302484000000007e+01, 8.2927484000000007e+01, 8.3552484000000007e+01, 8.4177484000000007e+01, 8.4802484000000007e+01, 8.6052484000000007e+01, 8.6677484000000007e+01, 8.7302484000000007e+01, 8.7927484000000007e+01, 8.8552484000000007e+01, 8.9177484000000007e+01, 8.9802484000000007e+01, 9.0427484000000007e+01, 9.1052484000000007e+01, 9.1677484000000007e+01, 9.2302484000000007e+01, 9.3552484000000007e+01, 9.4177484000000007e+01, 9.4802484000000007e+01, 9.5427484000000007e+01, 9.6052484000000007e+01, 9.6677484000000007e+01, 9.7302484000000007e+01, 9.7927484000000007e+01, 9.8552484000000007e+01, 9.9177484000000007e+01, 9.9802484000000007e+01, 1.0042748000000000e+02, 1.0105248000000000e+02, 1.0167748000000000e+02, 1.0230248000000000e+02, 1.0355248000000000e+02, 1.0355248000000000e+02, 1.0355248000000000e+02, 1.0355248000000000e+02};
static const casadi_real casadi_c1[223] = {-5.1238904394406278e-15, 7.4275199132374975e-09, -1.4855033738910648e-08, 2.2282544046342485e-08, -8.3559568618237473e-08, 3.1195572150899659e-07, -1.1642633093926858e-06, 4.3450974961510678e-06, -1.6216126682332229e-05, 6.0519409239803897e-05, -2.2586151028722404e-04, 8.4292663192242706e-04, -3.1458450173964905e-03, 1.1740453437652678e-02, -4.3815968733205271e-02, 1.6352342149516136e-01, -4.5359351034754869e-01, -3.9131804876800746e-01, -4.0232631525185575e-01, -3.9937666685974438e-01, -4.0016702239059498e-01, -3.9995524661694448e-01, -4.0001199166822327e-01, -3.9999678685126078e-01, -4.0000086095900106e-01, -3.9999976930571990e-01, -4.0000006181391473e-01, -3.9999998343802717e-01, -4.0000000443409739e-01, -3.9999999882566389e-01, -4.0000000026324867e-01, -4.0000000012135994e-01, -3.9999999925129132e-01, -4.0000000287349596e-01, -3.9999998925469682e-01, -4.0000004010751783e-01, -3.9999985031490054e-01, -4.0000055863439532e-01, -3.9999791515410094e-01, -4.0000778077371996e-01, -3.9997096189257714e-01, -4.0010837189269122e-01, -3.9959554944820674e-01, -4.0150942295524999e-01, -3.9436674642595926e-01, -4.2102364791962510e-01, -3.2153890926772066e-01, -6.9282163821605802e-01, 6.9282013405473020e-01, 3.2153871492188607e-01, 4.2102357860442297e-01, 3.9436675950587374e-01, 4.0150942929436528e-01, 3.9959555078169856e-01, 4.0010837164099622e-01, 3.9997096177087371e-01, 4.0000778074715215e-01, 3.9999791514896760e-01, 4.0000055863245354e-01, 3.9999985031559465e-01, 4.0000004010637291e-01, 3.9999998925964902e-01, 4.0000000285515397e-01, 3.9999999931969654e-01, 3.9999999986604334e-01, 4.0000000121612117e-01, 3.9999999526946967e-01, 4.0000001770588650e-01, 3.9999993390686034e-01, 4.0000024666732964e-01, 3.9999907942830787e-01, 4.0000343562694191e-01, 3.9998717802943845e-01, 4.0004785210441557e-01, 3.9982141298985413e-01, 4.0066649268667265e-01, 3.9751261083017769e-01, 4.0928308897522375e-01, 3.6535520218066680e-01, 5.2929638472715690e-01, -3.9814946604157186e-01, 1.4353535259800015e-01, -3.8460187384590776e-02, 1.0305396940353718e-02, -2.7614003768241502e-03, 7.4020456694596983e-04, -1.9941789096464954e-04, 5.7466996925705177e-05, -3.0450096743826066e-05, 6.4333390042746550e-05, -2.2688346342198990e-04, 8.4320046364711799e-04, -3.1459183911704404e-03, 1.1740473101041229e-02, -4.3815974012996084e-02, 1.6352342295092923e-01, -4.5359351048603291e-01, -3.9131804874554177e-01, -4.0232631525798240e-01, -3.9937666685772283e-01, -4.0016702239254326e-01, -3.9995524661118537e-01, -4.0001199161221368e-01, -3.9999678684220574e-01, -4.0000086096049309e-01, -3.9999976930567971e-01, -4.0000006181406683e-01, -3.9999998343744025e-01, -4.0000000443631517e-01, -3.9999999881738058e-01, -4.0000000029417604e-01, -4.0000000000591063e-01, -3.9999999968217537e-01, -4.0000000126538282e-01, -3.9999999525629182e-01, -4.0000001770949861e-01, -3.9999993390592736e-01, -4.0000024666757789e-01, -3.9999907942824720e-01, -4.0000343562694834e-01, -3.9998717802944023e-01, -4.0004785202197013e-01, -3.9982141349258721e-01, -4.0066649580135333e-01, -3.9751261114423736e-01, -4.0928308888938225e-01, -3.6535520220997919e-01, -5.2929638469575291e-01, 3.9814946587943661e-01, -1.4353535093793945e-01, 3.8460181369143379e-02, -1.0305374538616557e-02, 2.7613167853173417e-03, -7.3989260267002686e-04, 1.9825362536985924e-04, -5.3121898813468310e-05, 1.4233969885254099e-05, -3.8139807372009740e-06, 1.0219530629562317e-06, -2.7383150339603397e-07, 7.3372929220691309e-08, -1.9660213976388415e-08, 5.2679373432304501e-09, -1.4115372153278691e-09, 3.7822349639064066e-10, -1.0135543794193616e-10, 2.7189827517958113e-11, -7.4135964444517040e-12, 2.4787545695412166e-12, -2.5049148447897441e-12, 7.5550282879918693e-12, -2.7735153225941954e-11, 1.0338453598152428e-10, -3.8579728678878728e-10, 1.4398130052653584e-09, -5.3734527440964276e-09, 2.0053988447366531e-08, -7.4842487182859004e-08, 2.7931595579650842e-07, -1.0424213373942076e-06, 3.8903693978511740e-06, -1.4519056257585493e-05, 5.4185855637822188e-05, -2.0222436631364785e-04, 7.5471160963031078e-04, -2.8166220722168002e-03, 1.0511776679233713e-02, -3.9230484644701642e-02, 1.4641016189956638e-01, -5.4641016295357725e-01, -3.6076951008523178e-01, -4.1051179670550830e-01, -3.9718330309272454e-01, -4.0075499092359929e-01, -3.9979673321287934e-01, -4.0005807622488815e-01, -3.9997096188756054e-01, -4.0005807622486100e-01, -3.9979673321300219e-01, -4.0075499092312650e-01, -3.9718330309449817e-01, -4.1051179669886173e-01, -3.6076951011005748e-01, -5.4641016286091459e-01, 1.4641016155371589e-01, -3.9230483353952371e-02, 1.0511771862101676e-02, -2.8166040944650013e-03, 7.5464451577083919e-04, -2.0197396862252245e-04, 5.3251358703281429e-05, -1.1031466178164046e-05, -9.1254940015087092e-06, 4.7533442194959205e-05, -2.5718884710954261e-04, 1.6101881609976076e-03, -2.7535320699654973e-03, 1.0494870176602210e-02, -3.9225948636444720e-02, 1.4640892436916397e-01, -5.4640974884019611e-01, -3.6076992900838972e-01, -4.1051053512622504e-01, -3.9718793048672096e-01, -4.0073774292688086e-01, -3.9956888276539260e-01, -4.0006684893472683e-01, -3.9999079367288243e-01, -3.9998899018919476e-01, -4.0005324557034966e-01, -3.9979802752941418e-01, -4.0075464431196800e-01, -3.9718339522272178e-01, -4.1051177479714318e-01, -3.6076950558870946e-01, -5.4641020284802300e-01, 1.4640883896368467e-01, -3.9229335985231324e-02, 1.0507886274517108e-02, -2.8021030065367557e-03, 1.8680686710260523e-03, -9.3403433551443149e-04, 4.0294332838309812e-16};

/* kapparef_s:(x)->(f) */
static int casadi_f1(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real w0, w1;
  /* #0: @0 = input[0][0] */
  w0 = arg[0] ? arg[0][0] : 0;
  /* #1: @1 = BSpline(@0) */
  casadi_clear((&w1), 1);
  CASADI_PREFIX(nd_boor_eval)((&w1),1,casadi_c0,casadi_s0,casadi_s1,casadi_s2,casadi_c1,1,(&w0),casadi_s3, iw, w);
  /* #2: output[0][0] = @1 */
  if (res[0]) res[0][0] = w1;
  return 0;
}

/* tricycle_frenet_impl_dae_fun:(i0[8],i1[8],i2[2],i3[],i4[],i5[16])->(o0[8]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_int i;
  casadi_real **res1=res+1, *rr;
  const casadi_real **arg1=arg+6, *cs;
  casadi_real w0, w1, w2, w3, w4, w5, w6, w7, *w8=w+20, w9, w10, *w11=w+30;
  /* #0: @0 = input[1][0] */
  w0 = arg[1] ? arg[1][0] : 0;
  /* #1: @1 = input[1][1] */
  w1 = arg[1] ? arg[1][1] : 0;
  /* #2: @2 = input[1][2] */
  w2 = arg[1] ? arg[1][2] : 0;
  /* #3: @3 = input[1][3] */
  w3 = arg[1] ? arg[1][3] : 0;
  /* #4: @4 = input[1][4] */
  w4 = arg[1] ? arg[1][4] : 0;
  /* #5: @5 = input[1][5] */
  w5 = arg[1] ? arg[1][5] : 0;
  /* #6: @6 = input[1][6] */
  w6 = arg[1] ? arg[1][6] : 0;
  /* #7: @7 = input[1][7] */
  w7 = arg[1] ? arg[1][7] : 0;
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
  /* #9: @0 = input[0][6] */
  w0 = arg[0] ? arg[0][6] : 0;
  /* #10: @1 = input[0][7] */
  w1 = arg[0] ? arg[0][7] : 0;
  /* #11: @2 = cos(@1) */
  w2 = cos( w1 );
  /* #12: @2 = (@0*@2) */
  w2  = (w0*w2);
  /* #13: @3 = input[0][2] */
  w3 = arg[0] ? arg[0][2] : 0;
  /* #14: @4 = cos(@3) */
  w4 = cos( w3 );
  /* #15: @4 = (@2*@4) */
  w4  = (w2*w4);
  /* #16: @3 = sin(@3) */
  w3 = sin( w3 );
  /* #17: @3 = (@2*@3) */
  w3  = (w2*w3);
  /* #18: @1 = sin(@1) */
  w1 = sin( w1 );
  /* #19: @0 = (@0*@1) */
  w0 *= w1;
  /* #20: @1 = 1.03 */
  w1 = 1.0300000000000000e+00;
  /* #21: @0 = (@0/@1) */
  w0 /= w1;
  /* #22: @1 = input[0][5] */
  w1 = arg[0] ? arg[0][5] : 0;
  /* #23: @5 = cos(@1) */
  w5 = cos( w1 );
  /* #24: @5 = (@2*@5) */
  w5  = (w2*w5);
  /* #25: @6 = 1 */
  w6 = 1.;
  /* #26: @7 = input[0][3] */
  w7 = arg[0] ? arg[0][3] : 0;
  /* #27: @9 = kapparef_s(@7) */
  arg1[0]=(&w7);
  res1[0]=(&w9);
  if (casadi_f1(arg1, res1, iw, w, 0)) return 1;
  /* #28: @10 = input[0][4] */
  w10 = arg[0] ? arg[0][4] : 0;
  /* #29: @9 = (@9*@10) */
  w9 *= w10;
  /* #30: @6 = (@6-@9) */
  w6 -= w9;
  /* #31: @5 = (@5/@6) */
  w5 /= w6;
  /* #32: @1 = sin(@1) */
  w1 = sin( w1 );
  /* #33: @2 = (@2*@1) */
  w2 *= w1;
  /* #34: @1 = kapparef_s(@7) */
  arg1[0]=(&w7);
  res1[0]=(&w1);
  if (casadi_f1(arg1, res1, iw, w, 0)) return 1;
  /* #35: @1 = (@1*@5) */
  w1 *= w5;
  /* #36: @1 = (@0-@1) */
  w1  = (w0-w1);
  /* #37: @7 = input[2][0] */
  w7 = arg[2] ? arg[2][0] : 0;
  /* #38: @6 = input[2][1] */
  w6 = arg[2] ? arg[2][1] : 0;
  /* #39: @11 = vertcat(@4, @3, @0, @5, @2, @1, @7, @6) */
  rr=w11;
  *rr++ = w4;
  *rr++ = w3;
  *rr++ = w0;
  *rr++ = w5;
  *rr++ = w2;
  *rr++ = w1;
  *rr++ = w7;
  *rr++ = w6;
  /* #40: @8 = (@8-@11) */
  for (i=0, rr=w8, cs=w11; i<8; ++i) (*rr++) -= (*cs++);
  /* #41: output[0][0] = @8 */
  casadi_copy(w8, 8, res[0]);
  return 0;
}

CASADI_SYMBOL_EXPORT int tricycle_frenet_impl_dae_fun(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int tricycle_frenet_impl_dae_fun_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int tricycle_frenet_impl_dae_fun_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tricycle_frenet_impl_dae_fun_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int tricycle_frenet_impl_dae_fun_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tricycle_frenet_impl_dae_fun_release(int mem) {
}

CASADI_SYMBOL_EXPORT void tricycle_frenet_impl_dae_fun_incref(void) {
}

CASADI_SYMBOL_EXPORT void tricycle_frenet_impl_dae_fun_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int tricycle_frenet_impl_dae_fun_n_in(void) { return 6;}

CASADI_SYMBOL_EXPORT casadi_int tricycle_frenet_impl_dae_fun_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real tricycle_frenet_impl_dae_fun_default_in(casadi_int i) {
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tricycle_frenet_impl_dae_fun_name_in(casadi_int i) {
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    case 4: return "i4";
    case 5: return "i5";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tricycle_frenet_impl_dae_fun_name_out(casadi_int i) {
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tricycle_frenet_impl_dae_fun_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s4;
    case 1: return casadi_s4;
    case 2: return casadi_s5;
    case 3: return casadi_s6;
    case 4: return casadi_s6;
    case 5: return casadi_s7;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tricycle_frenet_impl_dae_fun_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s4;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int tricycle_frenet_impl_dae_fun_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 14;
  if (sz_res) *sz_res = 4;
  if (sz_iw) *sz_iw = 8;
  if (sz_w) *sz_w = 38;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
