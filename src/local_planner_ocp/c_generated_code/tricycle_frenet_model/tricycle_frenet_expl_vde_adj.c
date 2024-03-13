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
  #define CASADI_PREFIX(ID) tricycle_frenet_expl_vde_adj_ ## ID
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
#define casadi_c2 CASADI_PREFIX(c2)
#define casadi_c3 CASADI_PREFIX(c3)
#define casadi_clear CASADI_PREFIX(clear)
#define casadi_clear_casadi_int CASADI_PREFIX(clear_casadi_int)
#define casadi_copy CASADI_PREFIX(copy)
#define casadi_de_boor CASADI_PREFIX(de_boor)
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_f1 CASADI_PREFIX(f1)
#define casadi_f2 CASADI_PREFIX(f2)
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
#define casadi_s8 CASADI_PREFIX(s8)

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

static const casadi_int casadi_s0[2] = {0, 254};
static const casadi_int casadi_s1[1] = {2};
static const casadi_int casadi_s2[1] = {1};
static const casadi_int casadi_s3[2] = {0, 256};
static const casadi_int casadi_s4[1] = {3};
static const casadi_int casadi_s5[9] = {5, 1, 0, 5, 0, 1, 2, 3, 4};
static const casadi_int casadi_s6[6] = {2, 1, 0, 2, 0, 1};
static const casadi_int casadi_s7[20] = {16, 1, 0, 16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
static const casadi_int casadi_s8[11] = {7, 1, 0, 7, 0, 1, 2, 3, 4, 5, 6};

static const casadi_real casadi_c0[254] = {0., 0., 0., 5.0000000000000000e-01, 7.5000000000000000e-01, 1., 1.2500000000000000e+00, 1.5000000000000000e+00, 1.7500000000000000e+00, 2., 2.2500000000000000e+00, 2.5000000000000000e+00, 2.7500000000000000e+00, 3., 3.2500000000000000e+00, 3.5000000000000000e+00, 3.7500000000000000e+00, 4., 4.2500000000000000e+00, 4.5000000000000000e+00, 4.7500000000000000e+00, 5., 5.2500000000000000e+00, 5.5000000000000000e+00, 5.7500000000000000e+00, 6., 6.2500000000000000e+00, 6.5000000000000000e+00, 6.7500000000000000e+00, 7., 7.2500000000000000e+00, 7.5000000000000000e+00, 7.7500000000000000e+00, 8., 8.2500000000000000e+00, 8.5000000000000000e+00, 8.7500000000000000e+00, 9., 9.2500000000000000e+00, 9.5000000000000000e+00, 9.7500000000000000e+00, 10., 1.0250000000000000e+01, 1.0500000000000000e+01, 1.0750000000000000e+01, 11., 1.1250000000000000e+01, 1.1500000000000000e+01, 1.1750000000000000e+01, 12., 1.2250000000000000e+01, 1.2500000000000000e+01, 1.2750000000000000e+01, 13., 1.3250000000000000e+01, 1.3500000000000000e+01, 1.3750000000000000e+01, 14., 1.4250000000000000e+01, 1.4500000000000000e+01, 1.4750000000000000e+01, 15., 1.5250000000000000e+01, 1.5500000000000000e+01, 1.5750000000000000e+01, 16., 1.6250000000000000e+01, 1.6500000000000000e+01, 1.6750000000000000e+01, 17., 1.7250000000000000e+01, 1.7500000000000000e+01, 1.7750000000000000e+01, 18., 1.8250000000000000e+01, 1.8500000000000000e+01, 1.8750000000000000e+01, 19., 1.9250000000000000e+01, 1.9500000000000000e+01, 1.9750000000000000e+01, 20., 2.0250000000000000e+01, 2.0500000000000000e+01, 2.0750000000000000e+01, 21., 2.1250000000000000e+01, 2.1500000000000000e+01, 2.1750000000000000e+01, 22., 2.2250000000000000e+01, 2.2500000000000000e+01, 2.2750000000000000e+01, 23., 2.3250000000000000e+01, 2.3500000000000000e+01, 2.3750000000000000e+01, 24., 2.4250000000000000e+01, 2.4500000000000000e+01, 2.4750000000000000e+01, 25., 2.5250000000000000e+01, 2.5500000000000000e+01, 2.5750000000000000e+01, 26., 2.6250000000000000e+01, 2.6500000000000000e+01, 2.6750000000000000e+01, 27., 2.7250000000000000e+01, 2.7500000000000000e+01, 2.7750000000000000e+01, 28., 2.8250000000000000e+01, 2.8500000000000000e+01, 2.8750000000000000e+01, 29., 2.9250000000000000e+01, 2.9500000000000000e+01, 2.9750000000000000e+01, 30., 3.0250000000000000e+01, 3.0500000000000000e+01, 3.0750000000000000e+01, 31., 3.1250000000000000e+01, 3.1500000000000000e+01, 3.1750000000000000e+01, 32., 3.2250000000000000e+01, 3.2500000000000000e+01, 3.2750000000000000e+01, 33., 3.3250000000000000e+01, 3.3500000000000000e+01, 3.3750000000000000e+01, 34., 3.4250000000000000e+01, 3.4500000000000000e+01, 3.4750000000000000e+01, 35., 3.5250000000000000e+01, 3.5500000000000000e+01, 3.5750000000000000e+01, 36., 3.6250000000000000e+01, 3.6500000000000000e+01, 3.6750000000000000e+01, 37., 3.7250000000000000e+01, 3.7500000000000000e+01, 3.7750000000000000e+01, 38., 3.8250000000000000e+01, 3.8500000000000000e+01, 3.8750000000000000e+01, 39., 3.9250000000000000e+01, 3.9500000000000000e+01, 3.9750000000000000e+01, 40., 4.0250000000000000e+01, 4.0500000000000000e+01, 4.0750000000000000e+01, 41., 4.1250000000000000e+01, 4.1500000000000000e+01, 4.1750000000000000e+01, 42., 4.2250000000000000e+01, 4.2500000000000000e+01, 4.2750000000000000e+01, 43., 4.3250000000000000e+01, 4.3500000000000000e+01, 4.3750000000000000e+01, 44., 4.4250000000000000e+01, 4.4500000000000000e+01, 4.4750000000000000e+01, 45., 4.5250000000000000e+01, 4.5500000000000000e+01, 4.5750000000000000e+01, 46., 4.6250000000000000e+01, 4.6500000000000000e+01, 4.6750000000000000e+01, 47., 4.7250000000000000e+01, 4.7500000000000000e+01, 4.7750000000000000e+01, 48., 4.8250000000000000e+01, 4.8500000000000000e+01, 4.8750000000000000e+01, 49., 4.9250000000000000e+01, 4.9500000000000000e+01, 4.9750000000000000e+01, 50., 5.0250000000000000e+01, 5.0500000000000000e+01, 5.0750000000000000e+01, 51., 5.1250000000000000e+01, 5.1500000000000000e+01, 5.1750000000000000e+01, 52., 5.2250000000000000e+01, 5.2500000000000000e+01, 5.2750000000000000e+01, 53., 5.3250000000000000e+01, 5.3500000000000000e+01, 5.3750000000000000e+01, 54., 5.4250000000000000e+01, 5.4500000000000000e+01, 5.4750000000000000e+01, 55., 5.5250000000000000e+01, 5.5500000000000000e+01, 5.5750000000000000e+01, 56., 5.6250000000000000e+01, 5.6500000000000000e+01, 5.6750000000000000e+01, 57., 5.7250000000000000e+01, 5.7500000000000000e+01, 5.7750000000000000e+01, 58., 5.8250000000000000e+01, 5.8500000000000000e+01, 5.8750000000000000e+01, 59., 5.9250000000000000e+01, 5.9500000000000000e+01, 5.9750000000000000e+01, 60., 6.0250000000000000e+01, 6.0500000000000000e+01, 6.0750000000000000e+01, 61., 6.1250000000000000e+01, 6.1500000000000000e+01, 6.1750000000000000e+01, 62., 6.2250000000000000e+01, 6.2750000000000000e+01, 6.2750000000000000e+01, 6.2750000000000000e+01};
static const casadi_real casadi_c1[251] = {1.1833456234363895e-17, -1.0555468703052595e-15, 1.7843280782391468e-17, 4.0910660920777408e-15, 4.7549193836826801e-15, 1.9022390207010821e-15, -4.9706646798780596e-15, -9.0642166970977375e-15, -1.9477472667329016e-15, 9.7963930430839043e-15, 5.6503933313518111e-15, -1.1632699454868932e-14, -7.0056370334316835e-15, 1.6171971959966925e-14, 4.4421239557011827e-15, -2.0838502805252252e-14, 6.2522169536411763e-15, 1.7457757755592623e-14, -2.0952303964202297e-14, 4.3222575999042431e-15, 1.5075092659740806e-14, -2.6199706268286423e-14, 2.4168655697398291e-14, -6.1293627189561108e-15, -2.0990288014001185e-14, 4.1168884686160216e-14, -4.3595139912726591e-14, 3.2153842550134029e-14, -1.7091251904655817e-14, 5.0474846670271758e-15, -2.2537086279636764e-17, 6.2209970384576567e-15, -1.9587727726583707e-14, 2.8558726781067208e-14, -4.0323618078056844e-14, 1.1901121981725430e-13, -4.5575756121987154e-13, 1.6982227181027287e-12, -6.3027906255953023e-12, 2.3515598596831988e-11, -8.7800563108583474e-11, 3.2769239277043139e-10, -1.2229269659368150e-09, 4.5640107163259931e-09, -1.7033150496497642e-08, 6.3568583429406067e-08, -2.3724116919554574e-07, 8.8539609774134944e-07, -3.3043432291068531e-06, 1.2331976817437176e-05, -4.6023564030257351e-05, 1.7176227931390835e-04, -6.4102555322483891e-04, 2.3923399335793018e-03, -8.9283341810976937e-03, 3.3320996790809701e-02, -1.2435565298214085e-01, 4.6410161513775472e-01, -1.7320508075688774e+00, 4.6410161513775461e-01, -1.2435565298214057e-01, 3.3320996790809798e-02, -8.9283341810978412e-03, 2.3923399335794571e-03, -6.4102555322453902e-04, 1.7176227931359200e-04, -4.6023564029162856e-05, 1.2331976817048229e-05, -3.3043432297041875e-06, 8.8539609455207113e-07, -2.3724117348411511e-07, 6.3568581065709395e-08, -1.7033144894540442e-08, 4.5640201618013521e-09, -1.2229256496354424e-09, 3.2768265878502234e-10, -8.7805429593856843e-11, 2.3526958159436617e-11, -6.2951865942295626e-12, 1.6810997038874120e-12, -4.6007642140466487e-13, 1.3955503419538218e-13, -4.5408121707168903e-14, 1.0547118733938987e-14, 1.4432899320127035e-15, 1.6653345369377348e-15, -1.4543921622589551e-14, 3.0642155479654321e-14, -4.0412118096355698e-14, 3.7858605139717838e-14, -2.2759572004815709e-14, 0., 2.2759572004815709e-14, -3.7858605139717838e-14, 4.0412118096355698e-14, -3.0642155479654321e-14, 1.4543921622589551e-14, -1.6653345369377348e-15, -1.4432899320127035e-15, -1.0547118733938987e-14, 4.5408121707168903e-14, -1.3955503419538218e-13, 4.6007642140466487e-13, -1.6810997038874120e-12, 6.2951865942295626e-12, -2.3526958159436617e-11, 8.7805429593856843e-11, -3.2768265878502234e-10, 1.2229256496354424e-09, -4.5640201618013521e-09, 1.7033144894540442e-08, -6.3568581065709395e-08, 2.3724117348411511e-07, -8.8539609455207113e-07, 3.3043432297041875e-06, -1.2331976817048229e-05, 4.6023564029162856e-05, -1.7176227931359200e-04, 6.4102555322453902e-04, -2.3923399335794571e-03, 8.9283341810978412e-03, -3.3320996790810020e-02, 1.2435565298214113e-01, -4.6410161513775472e-01, 1.7320508075688774e+00, -4.6410161513775444e-01, 1.2435565298214132e-01, -3.3320996790809687e-02, 8.9283341810966685e-03, -2.3923399335823233e-03, 6.4102555322067428e-04, -1.7176227931550131e-04, 4.6023564035358129e-05, -1.2331976808325675e-05, 3.3043432310690201e-06, -8.8539610753385373e-07, 2.3724116354613130e-07, -6.3568571796328268e-08, 1.7033157502053923e-08, -4.5640268882640606e-09, 1.2229225238164480e-09, -3.2767155427608324e-10, 8.7794310904760468e-11, -2.3533056353879196e-11, 6.3237429245900984e-12, -1.7025449751610884e-12, 4.4068246784999029e-13, -9.2811512612451441e-14, 1.6154961917642753e-14, -2.2429363892086403e-14, 4.0578015689480637e-14, -4.7389881724444369e-14, 4.3617676999211918e-14, -3.7201327210543949e-14, 3.4182503809954172e-14, -3.7201327210052129e-14, 4.3617677003597442e-14, -4.7389881761103013e-14, 4.0578015798147937e-14, -2.2429364090282549e-14, 1.6154962449669868e-14, -9.2811513309836414e-14, 4.4068246824184949e-13, -1.7025449757524050e-12, 6.3237429296281332e-12, -2.3533056354544382e-11, 8.7794310891577458e-11, -3.2767155426741717e-10, 1.2229225238141044e-09, -4.5640268882903500e-09, 1.7033157502177119e-08, -6.3568571796491997e-08, 2.3724116354613075e-07, -8.8539610753374149e-07, 3.3043432310691641e-06, -1.2331976808326161e-05, 4.6023564035358427e-05, -1.7176227931549809e-04, 6.4102555322066799e-04, -2.3923399335823146e-03, 8.9283341810966390e-03, -3.3320996790809618e-02, 1.2435565298214125e-01, -4.6410161513775439e-01, 1.7320508075688774e+00, -4.6410161513775439e-01, 1.2435565298214057e-01, -3.3320996790809798e-02, 8.9283341810978412e-03, -2.3923399335786799e-03, 6.4102555322576027e-04, -1.7176227931237076e-04, 4.6023564030162056e-05, -1.2331976817492318e-05, 3.3043432259294292e-06, -8.8539609910398553e-07, 2.3724116959833452e-07, -6.3568579400374858e-08, 1.7033149002365633e-08, -4.5640188295337225e-09, 1.2229254275908374e-09, -3.2768054936127555e-10, 8.7806428794579006e-11, -2.3532287229954818e-11, 6.2938543266000124e-12, -1.6718848527830232e-12, 4.5818904226280210e-13, -1.5043521983670871e-13, 5.3623772089395061e-14, -4.4408920985006262e-15, -9.2148511043887993e-15, -2.5535129566378600e-15, 1.5765166949677223e-14, -2.3203661214665772e-14, 3.1419311596891930e-14, -4.0301095793893182e-14, 4.0079051188968151e-14, -2.6645352591003757e-14, 6.6613381477509392e-15, 1.0880185641326534e-14, -2.0983215165415459e-14, 1.9428902930940239e-14, -6.6613381477509392e-15, -6.5503158452884236e-15, 1.0103029524088925e-14, -4.8849813083506888e-15, -2.3314683517128287e-15, 8.1046280797636427e-15, -7.1054273576010019e-15, -1.2212453270876722e-15, 7.4384942649885488e-15, -5.5511151231257827e-15, -3.6637359812630166e-15, 5.3290705182007514e-15, -2.3314683517128287e-15, -2.5535129566378600e-15, 1.0325074129013956e-14, 5.7731597280508140e-15, -6.4392935428259079e-15, -7.4384942649885488e-15, 5.5511151231257827e-16, 2.4424906541753444e-15, -2.2204460492503131e-16, -1.6653345369377348e-15, -1.1102230246251565e-16, 8.8817841970012523e-16, 7.7715611723760958e-16, -5.5511151231257827e-16, 2.7755575615628914e-16, 6.6613381477509392e-16, 1.9984014443252818e-15};
static const casadi_real casadi_c2[256] = {0., 0., 0., 0., 5.0000000000000000e-01, 7.5000000000000000e-01, 1., 1.2500000000000000e+00, 1.5000000000000000e+00, 1.7500000000000000e+00, 2., 2.2500000000000000e+00, 2.5000000000000000e+00, 2.7500000000000000e+00, 3., 3.2500000000000000e+00, 3.5000000000000000e+00, 3.7500000000000000e+00, 4., 4.2500000000000000e+00, 4.5000000000000000e+00, 4.7500000000000000e+00, 5., 5.2500000000000000e+00, 5.5000000000000000e+00, 5.7500000000000000e+00, 6., 6.2500000000000000e+00, 6.5000000000000000e+00, 6.7500000000000000e+00, 7., 7.2500000000000000e+00, 7.5000000000000000e+00, 7.7500000000000000e+00, 8., 8.2500000000000000e+00, 8.5000000000000000e+00, 8.7500000000000000e+00, 9., 9.2500000000000000e+00, 9.5000000000000000e+00, 9.7500000000000000e+00, 10., 1.0250000000000000e+01, 1.0500000000000000e+01, 1.0750000000000000e+01, 11., 1.1250000000000000e+01, 1.1500000000000000e+01, 1.1750000000000000e+01, 12., 1.2250000000000000e+01, 1.2500000000000000e+01, 1.2750000000000000e+01, 13., 1.3250000000000000e+01, 1.3500000000000000e+01, 1.3750000000000000e+01, 14., 1.4250000000000000e+01, 1.4500000000000000e+01, 1.4750000000000000e+01, 15., 1.5250000000000000e+01, 1.5500000000000000e+01, 1.5750000000000000e+01, 16., 1.6250000000000000e+01, 1.6500000000000000e+01, 1.6750000000000000e+01, 17., 1.7250000000000000e+01, 1.7500000000000000e+01, 1.7750000000000000e+01, 18., 1.8250000000000000e+01, 1.8500000000000000e+01, 1.8750000000000000e+01, 19., 1.9250000000000000e+01, 1.9500000000000000e+01, 1.9750000000000000e+01, 20., 2.0250000000000000e+01, 2.0500000000000000e+01, 2.0750000000000000e+01, 21., 2.1250000000000000e+01, 2.1500000000000000e+01, 2.1750000000000000e+01, 22., 2.2250000000000000e+01, 2.2500000000000000e+01, 2.2750000000000000e+01, 23., 2.3250000000000000e+01, 2.3500000000000000e+01, 2.3750000000000000e+01, 24., 2.4250000000000000e+01, 2.4500000000000000e+01, 2.4750000000000000e+01, 25., 2.5250000000000000e+01, 2.5500000000000000e+01, 2.5750000000000000e+01, 26., 2.6250000000000000e+01, 2.6500000000000000e+01, 2.6750000000000000e+01, 27., 2.7250000000000000e+01, 2.7500000000000000e+01, 2.7750000000000000e+01, 28., 2.8250000000000000e+01, 2.8500000000000000e+01, 2.8750000000000000e+01, 29., 2.9250000000000000e+01, 2.9500000000000000e+01, 2.9750000000000000e+01, 30., 3.0250000000000000e+01, 3.0500000000000000e+01, 3.0750000000000000e+01, 31., 3.1250000000000000e+01, 3.1500000000000000e+01, 3.1750000000000000e+01, 32., 3.2250000000000000e+01, 3.2500000000000000e+01, 3.2750000000000000e+01, 33., 3.3250000000000000e+01, 3.3500000000000000e+01, 3.3750000000000000e+01, 34., 3.4250000000000000e+01, 3.4500000000000000e+01, 3.4750000000000000e+01, 35., 3.5250000000000000e+01, 3.5500000000000000e+01, 3.5750000000000000e+01, 36., 3.6250000000000000e+01, 3.6500000000000000e+01, 3.6750000000000000e+01, 37., 3.7250000000000000e+01, 3.7500000000000000e+01, 3.7750000000000000e+01, 38., 3.8250000000000000e+01, 3.8500000000000000e+01, 3.8750000000000000e+01, 39., 3.9250000000000000e+01, 3.9500000000000000e+01, 3.9750000000000000e+01, 40., 4.0250000000000000e+01, 4.0500000000000000e+01, 4.0750000000000000e+01, 41., 4.1250000000000000e+01, 4.1500000000000000e+01, 4.1750000000000000e+01, 42., 4.2250000000000000e+01, 4.2500000000000000e+01, 4.2750000000000000e+01, 43., 4.3250000000000000e+01, 4.3500000000000000e+01, 4.3750000000000000e+01, 44., 4.4250000000000000e+01, 4.4500000000000000e+01, 4.4750000000000000e+01, 45., 4.5250000000000000e+01, 4.5500000000000000e+01, 4.5750000000000000e+01, 46., 4.6250000000000000e+01, 4.6500000000000000e+01, 4.6750000000000000e+01, 47., 4.7250000000000000e+01, 4.7500000000000000e+01, 4.7750000000000000e+01, 48., 4.8250000000000000e+01, 4.8500000000000000e+01, 4.8750000000000000e+01, 49., 4.9250000000000000e+01, 4.9500000000000000e+01, 4.9750000000000000e+01, 50., 5.0250000000000000e+01, 5.0500000000000000e+01, 5.0750000000000000e+01, 51., 5.1250000000000000e+01, 5.1500000000000000e+01, 5.1750000000000000e+01, 52., 5.2250000000000000e+01, 5.2500000000000000e+01, 5.2750000000000000e+01, 53., 5.3250000000000000e+01, 5.3500000000000000e+01, 5.3750000000000000e+01, 54., 5.4250000000000000e+01, 5.4500000000000000e+01, 5.4750000000000000e+01, 55., 5.5250000000000000e+01, 5.5500000000000000e+01, 5.5750000000000000e+01, 56., 5.6250000000000000e+01, 5.6500000000000000e+01, 5.6750000000000000e+01, 57., 5.7250000000000000e+01, 5.7500000000000000e+01, 5.7750000000000000e+01, 58., 5.8250000000000000e+01, 5.8500000000000000e+01, 5.8750000000000000e+01, 59., 5.9250000000000000e+01, 5.9500000000000000e+01, 5.9750000000000000e+01, 60., 6.0250000000000000e+01, 6.0500000000000000e+01, 6.0750000000000000e+01, 61., 6.1250000000000000e+01, 6.1500000000000000e+01, 6.1750000000000000e+01, 62., 6.2250000000000000e+01, 6.2750000000000000e+01, 6.2750000000000000e+01, 6.2750000000000000e+01, 6.2750000000000000e+01};
static const casadi_real casadi_c3[252] = {-2.8041110351460460e-16, -2.7843886080887727e-16, -5.4232557838519215e-16, -5.3637781812439499e-16, 4.8638870489504022e-16, 1.6751185508157103e-15, 2.1506783059909809e-15, 9.0801213602146596e-16, -1.3580420382529684e-15, -1.8449788549361938e-15, 6.0411940583478213e-16, 2.0167177386727350e-15, -8.9145712504449789e-16, -2.6428663834024188e-15, 1.4001266065893123e-15, 2.5106575955146080e-15, -2.6989681057984551e-15, -1.1359138673881611e-15, 3.2285255715099950e-15, -2.0095504195405792e-15, -9.2898601956451847e-16, 2.8397871453706826e-15, -3.7101394217009230e-15, 2.3320245026486494e-15, 7.9968382290962156e-16, -4.4478881805906747e-15, 5.8443329909493786e-15, -5.0544519872322692e-15, 2.9840086503012380e-15, -1.2888043258627162e-15, -2.6933159105922219e-17, -3.2567430675831410e-17, 1.5226818289385827e-15, -3.3742501027073441e-15, 3.7654315925594583e-15, -6.3154729269547519e-15, 2.3437332027358826e-14, -9.0502058277609052e-14, 3.3405362124807313e-13, -1.2416440351507524e-12, 4.6372556140572441e-12, -1.7312885163088624e-11, 6.4610213029519217e-11, -2.4112152845468453e-10, 8.9988115062681385e-10, -3.3584064734975967e-09, 1.2533739383853919e-08, -4.6776552915032513e-08, 1.7457247152030484e-07, -6.5151333575640848e-07, 2.4314808686028858e-06, -9.0744101389614519e-06, 3.3866159689515636e-05, -1.2639022861669409e-04, 4.7169475477813142e-04, -1.7603887904962920e-03, 6.5698604072061337e-03, -2.4519052838329079e-02, 9.1506350946109594e-02, -3.4150635094610976e-01, -2.2548094716167111e-01, -2.5656986040720625e-01, -2.4823961120950380e-01, -2.5047169475477826e-01, -2.4987360977138340e-01, -2.5003386615968953e-01, -2.4999092558986113e-01, -2.5000243148086843e-01, -2.4999934848666416e-01, -2.5000017457247159e-01, -2.4999995322344795e-01, -2.5000001253374132e-01, -2.4999999664159606e-01, -2.5000000089988228e-01, -2.4999999975887724e-01, -2.5000000006460865e-01, -2.4999999998268799e-01, -2.5000000000463934e-01, -2.4999999999875760e-01, -2.5000000000033140e-01, -2.4999999999991113e-01, -2.5000000000002615e-01, -2.4999999999999126e-01, -2.5000000000000261e-01, -2.4999999999999997e-01, -2.4999999999999961e-01, -2.4999999999999920e-01, -2.5000000000000283e-01, -2.4999999999999517e-01, -2.5000000000000527e-01, -2.4999999999999581e-01, -2.5000000000000150e-01, -2.5000000000000150e-01, -2.4999999999999581e-01, -2.5000000000000527e-01, -2.4999999999999517e-01, -2.5000000000000283e-01, -2.4999999999999920e-01, -2.4999999999999961e-01, -2.4999999999999997e-01, -2.5000000000000261e-01, -2.4999999999999126e-01, -2.5000000000002615e-01, -2.4999999999991113e-01, -2.5000000000033140e-01, -2.4999999999875760e-01, -2.5000000000463934e-01, -2.4999999998268799e-01, -2.5000000006460865e-01, -2.4999999975887724e-01, -2.5000000089988228e-01, -2.4999999664159606e-01, -2.5000001253374132e-01, -2.4999995322344795e-01, -2.5000017457247159e-01, -2.4999934848666416e-01, -2.5000243148086843e-01, -2.4999092558986113e-01, -2.5003386615968953e-01, -2.4987360977138340e-01, -2.5047169475477826e-01, -2.4823961120950380e-01, -2.5656986040720631e-01, -2.2548094716167102e-01, -3.4150635094610970e-01, 9.1506350946109663e-02, -2.4519052838328951e-02, 6.5698604072063809e-03, -1.7603887904960411e-03, 4.7169475477812616e-04, -1.2639022861745463e-04, 3.3866159687713943e-05, -9.0744101411613862e-06, 2.4314808676781465e-06, -6.5151333440327219e-07, 1.7457247336398279e-07, -4.6776553519480624e-08, 1.2533737367052202e-08, -3.3584055820298644e-09, 8.9988379348361596e-10, -2.4112292858239918e-10, 6.4607702371712829e-11, -1.7310186197307980e-11, 4.6383915288821372e-12, -1.2448725595876620e-12, 3.3606317155986261e-13, -8.9573072230409502e-14, 2.0597544732088074e-14, -2.6053334210247870e-15, 1.4334070583859009e-15, -4.1739339146356997e-15, 5.9705700077344604e-15, -5.8769004233766328e-15, 5.0275188264263468e-15, -4.2728129762096397e-15, 4.2728129762789032e-15, -5.0275188262341290e-15, 5.8769004246652315e-15, -5.9705700156105211e-15, 4.1739339339264633e-15, -1.4334070886441740e-15, 2.6053335237732931e-15, -2.0597544803685810e-14, 8.9573072256776566e-14, -3.3606317168132467e-13, 1.2448725607257087e-12, -4.6383915279103866e-12, 1.7310186194983977e-11, -6.4607702371870317e-11, 2.4112292858165580e-10, -8.9988379349093176e-10, 3.3584055820533480e-09, -1.2533737367069651e-08, 4.6776553519463035e-08, -1.7457247336397234e-07, 6.5151333440331867e-07, -2.4314808676782219e-06, 9.0744101411613845e-06, -3.3866159687713137e-05, 1.2639022861745387e-04, -4.7169475477812475e-04, 1.7603887904960351e-03, -6.5698604072063696e-03, 2.4519052838328944e-02, -9.1506350946109649e-02, 3.4150635094610970e-01, 2.2548094716167111e-01, 2.5656986040720625e-01, 2.4823961120950380e-01, 2.5047169475477826e-01, 2.4987360977138359e-01, 2.5003386615969003e-01, 2.4999092558986194e-01, 2.5000243148086948e-01, 2.4999934848666511e-01, 2.5000017457247159e-01, 2.4999995322344681e-01, 2.5000001253373921e-01, 2.4999999664159436e-01, 2.5000000089988161e-01, 2.4999999975887691e-01, 2.5000000006460826e-01, 2.4999999998268813e-01, 2.5000000000463973e-01, 2.4999999999875666e-01, 2.5000000000033012e-01, 2.4999999999991215e-01, 2.5000000000002670e-01, 2.4999999999998909e-01, 2.5000000000000250e-01, 2.5000000000000139e-01, 2.4999999999999908e-01, 2.4999999999999845e-01, 2.5000000000000239e-01, 2.4999999999999659e-01, 2.5000000000000444e-01, 2.4999999999999437e-01, 2.5000000000000439e-01, 2.4999999999999772e-01, 2.4999999999999939e-01, 2.5000000000000211e-01, 2.4999999999999686e-01, 2.5000000000000172e-01, 2.5000000000000006e-01, 2.4999999999999842e-01, 2.5000000000000094e-01, 2.4999999999999972e-01, 2.4999999999999914e-01, 2.5000000000000117e-01, 2.4999999999999939e-01, 2.4999999999999908e-01, 2.5000000000000094e-01, 2.4999999999999956e-01, 2.4999999999999864e-01, 2.4999999999999997e-01, 2.4999999999999939e-01, 2.4999999999999875e-01, 2.5000000000000133e-01, 2.5000000000000278e-01, 2.5000000000000117e-01, 2.4999999999999931e-01, 2.4999999999999944e-01, 2.5000000000000006e-01, 2.5000000000000000e-01, 2.4999999999999958e-01, 2.4999999999999956e-01, 2.4999999999999978e-01, 2.4999999999999997e-01, 2.4999999999999983e-01, 2.4999999999999994e-01, 2.5000000000000011e-01, 2.5000000000000044e-01};

/* jac_kapparef_s:(x,out_f[1x1,0nz])->(jac_f_x) */
static int casadi_f1(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real w0, w1;
  /* #0: @0 = input[0][0] */
  w0 = arg[0] ? arg[0][0] : 0;
  /* #1: @1 = BSpline(@0) */
  casadi_clear((&w1), 1);
  CASADI_PREFIX(nd_boor_eval)((&w1),1,casadi_c0,casadi_s0,casadi_s1,casadi_s2,casadi_c1,1,(&w0),casadi_s1, iw, w);
  /* #2: output[0][0] = @1 */
  if (res[0]) res[0][0] = w1;
  return 0;
}

/* kapparef_s:(x)->(f) */
static int casadi_f2(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real w0, w1;
  /* #0: @0 = input[0][0] */
  w0 = arg[0] ? arg[0][0] : 0;
  /* #1: @1 = BSpline(@0) */
  casadi_clear((&w1), 1);
  CASADI_PREFIX(nd_boor_eval)((&w1),1,casadi_c2,casadi_s3,casadi_s4,casadi_s2,casadi_c3,1,(&w0),casadi_s1, iw, w);
  /* #2: output[0][0] = @1 */
  if (res[0]) res[0][0] = w1;
  return 0;
}

/* tricycle_frenet_expl_vde_adj:(i0[5],i1[5],i2[2],i3[16])->(o0[7]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real **res1=res+1;
  const casadi_real **arg1=arg+4;
  casadi_real w0, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, *w14=w+25, w15, w16, w17, w18, w19, w20;
  /* #0: @0 = input[0][0] */
  w0 = arg[0] ? arg[0][0] : 0;
  /* #1: @1 = 00 */
  /* #2: @2 = jac_kapparef_s(@0, @1) */
  arg1[0]=(&w0);
  arg1[1]=0;
  res1[0]=(&w2);
  if (casadi_f1(arg1, res1, iw, w, 0)) return 1;
  /* #3: @3 = input[0][1] */
  w3 = arg[0] ? arg[0][1] : 0;
  /* #4: @4 = input[0][3] */
  w4 = arg[0] ? arg[0][3] : 0;
  /* #5: @5 = input[0][4] */
  w5 = arg[0] ? arg[0][4] : 0;
  /* #6: @6 = cos(@5) */
  w6 = cos( w5 );
  /* #7: @7 = (@4*@6) */
  w7  = (w4*w6);
  /* #8: @8 = input[0][2] */
  w8 = arg[0] ? arg[0][2] : 0;
  /* #9: @9 = cos(@8) */
  w9 = cos( w8 );
  /* #10: @10 = (@7*@9) */
  w10  = (w7*w9);
  /* #11: @11 = 1 */
  w11 = 1.;
  /* #12: @12 = kapparef_s(@0) */
  arg1[0]=(&w0);
  res1[0]=(&w12);
  if (casadi_f2(arg1, res1, iw, w, 0)) return 1;
  /* #13: @13 = (@12*@3) */
  w13  = (w12*w3);
  /* #14: @11 = (@11-@13) */
  w11 -= w13;
  /* #15: @10 = (@10/@11) */
  w10 /= w11;
  /* #16: @13 = (@10/@11) */
  w13  = (w10/w11);
  /* #17: @14 = input[1][0] */
  casadi_copy(arg[1], 5, w14);
  /* #18: {@15, @16, @17, @18, @19} = vertsplit(@14) */
  w15 = w14[0];
  w16 = w14[1];
  w17 = w14[2];
  w18 = w14[3];
  w19 = w14[4];
  /* #19: @20 = kapparef_s(@0) */
  arg1[0]=(&w0);
  res1[0]=(&w20);
  if (casadi_f2(arg1, res1, iw, w, 0)) return 1;
  /* #20: @20 = (@20*@17) */
  w20 *= w17;
  /* #21: @15 = (@15-@20) */
  w15 -= w20;
  /* #22: @13 = (@13*@15) */
  w13 *= w15;
  /* #23: @3 = (@3*@13) */
  w3 *= w13;
  /* #24: @2 = (@2*@3) */
  w2 *= w3;
  /* #25: @1 = 00 */
  /* #26: @3 = jac_kapparef_s(@0, @1) */
  arg1[0]=(&w0);
  arg1[1]=0;
  res1[0]=(&w3);
  if (casadi_f1(arg1, res1, iw, w, 0)) return 1;
  /* #27: @10 = (@10*@17) */
  w10 *= w17;
  /* #28: @3 = (@3*@10) */
  w3 *= w10;
  /* #29: @2 = (@2-@3) */
  w2 -= w3;
  /* #30: output[0][0] = @2 */
  if (res[0]) res[0][0] = w2;
  /* #31: @12 = (@12*@13) */
  w12 *= w13;
  /* #32: output[0][1] = @12 */
  if (res[0]) res[0][1] = w12;
  /* #33: @12 = cos(@8) */
  w12 = cos( w8 );
  /* #34: @13 = (@7*@16) */
  w13  = (w7*w16);
  /* #35: @12 = (@12*@13) */
  w12 *= w13;
  /* #36: @13 = sin(@8) */
  w13 = sin( w8 );
  /* #37: @15 = (@15/@11) */
  w15 /= w11;
  /* #38: @7 = (@7*@15) */
  w7 *= w15;
  /* #39: @13 = (@13*@7) */
  w13 *= w7;
  /* #40: @12 = (@12-@13) */
  w12 -= w13;
  /* #41: output[0][2] = @12 */
  if (res[0]) res[0][2] = w12;
  /* #42: @12 = sin(@5) */
  w12 = sin( w5 );
  /* #43: @13 = 0.970874 */
  w13 = 9.7087378640776700e-01;
  /* #44: @13 = (@13*@17) */
  w13 *= w17;
  /* #45: @12 = (@12*@13) */
  w12 *= w13;
  /* #46: @8 = sin(@8) */
  w8 = sin( w8 );
  /* #47: @8 = (@8*@16) */
  w8 *= w16;
  /* #48: @9 = (@9*@15) */
  w9 *= w15;
  /* #49: @8 = (@8+@9) */
  w8 += w9;
  /* #50: @6 = (@6*@8) */
  w6 *= w8;
  /* #51: @12 = (@12+@6) */
  w12 += w6;
  /* #52: output[0][3] = @12 */
  if (res[0]) res[0][3] = w12;
  /* #53: @12 = cos(@5) */
  w12 = cos( w5 );
  /* #54: @13 = (@4*@13) */
  w13  = (w4*w13);
  /* #55: @12 = (@12*@13) */
  w12 *= w13;
  /* #56: @5 = sin(@5) */
  w5 = sin( w5 );
  /* #57: @4 = (@4*@8) */
  w4 *= w8;
  /* #58: @5 = (@5*@4) */
  w5 *= w4;
  /* #59: @12 = (@12-@5) */
  w12 -= w5;
  /* #60: output[0][4] = @12 */
  if (res[0]) res[0][4] = w12;
  /* #61: output[0][5] = @18 */
  if (res[0]) res[0][5] = w18;
  /* #62: output[0][6] = @19 */
  if (res[0]) res[0][6] = w19;
  return 0;
}

CASADI_SYMBOL_EXPORT int tricycle_frenet_expl_vde_adj(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int tricycle_frenet_expl_vde_adj_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int tricycle_frenet_expl_vde_adj_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tricycle_frenet_expl_vde_adj_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int tricycle_frenet_expl_vde_adj_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tricycle_frenet_expl_vde_adj_release(int mem) {
}

CASADI_SYMBOL_EXPORT void tricycle_frenet_expl_vde_adj_incref(void) {
}

CASADI_SYMBOL_EXPORT void tricycle_frenet_expl_vde_adj_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int tricycle_frenet_expl_vde_adj_n_in(void) { return 4;}

CASADI_SYMBOL_EXPORT casadi_int tricycle_frenet_expl_vde_adj_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real tricycle_frenet_expl_vde_adj_default_in(casadi_int i) {
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tricycle_frenet_expl_vde_adj_name_in(casadi_int i) {
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tricycle_frenet_expl_vde_adj_name_out(casadi_int i) {
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tricycle_frenet_expl_vde_adj_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s5;
    case 1: return casadi_s5;
    case 2: return casadi_s6;
    case 3: return casadi_s7;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tricycle_frenet_expl_vde_adj_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s8;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int tricycle_frenet_expl_vde_adj_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 7;
  if (sz_res) *sz_res = 6;
  if (sz_iw) *sz_iw = 8;
  if (sz_w) *sz_w = 36;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
