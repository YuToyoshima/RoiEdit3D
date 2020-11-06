// elipsoid_em_22.cpp for Linux 

/*
 Comments written in Japanese are removed for compiling in Linux.
 Please see elipsoid_em_9.cpp for windows version.
 
 em_15: test for regularization of \mu 
 em_16: another test for regularization of \mu
        lambda_mu was specified by lambda_mu_coeff, lambda_mu_diag, and lambda_mu_nondiag.
        diagonal element of lambda_mu equals lambda_mu_diag.
        nondiagonal element of lambda_mu equals
           lambda_mu_nondiag / (||\mu_i-\mu_j||^2)^lambda_mu_coeff.
        Caution! this version is not stable. Do not use.
 em_17: Calculate Hessian explicitly based on elipsoid_mex_26 and others.
        Porting from elipsoid_em_15 and elipsoid_mex_26.
        Introducing step-size control based on Levenberg-Marquardt method.
        Fitting a mixture of 2 ellipsoids are dramatically improved by the step-size control.
        note: Mean consuming time was 3.2e9 clocks / iteration.
        note: Explicit calculation of hessian was not time consuming (~4.5e7 clocks).
        note: dsysv in LM step for sigma maybe time consuming (3e8~4e8 clocks).
              -> Use sparse matrix?
        note: calc_ofv in LM step (called twice) is time consuming (2.4e8 clocks).
              -> Trust-region method may be fast because that requires single call.
                 But solving trust-region subproblem may take time.
              -> If LM parameter (lmp) becomes too small, skip the call for small lmp?
        note: calc_ofv for update pi_k is slower than calc_H (2.4e8 vs 1.3e8 clocks).
              -> Use calc_H?
 em_18: Introducing L2 regularization for mu and sigma
 em_19: Introducing L2 regularization for pi_k
        For improving stability, use dgelsd instead of dsysv
 em_20: For improving performance use dgelsy instead of dgelsd
 em_21: Changing L2 regularization term for mu.
        old: inner product of moving vectors 
             sum_i(sum_j(lambda_mu_ij* (mu_i-mu_i^init)*(mu_j-mu_j^init) ))
        new: square of difference of moving vectors
             sum_i(sum_j(lambda_mu_ij* ((mu_i-mu_i^init)-(mu_j-mu_j^init))^2 ))
        With this change, summation of lambda_mu is calculated in init_params().
        Also, diagonal elements of lambda_mu become useless and set to 0.
        This version may include bug in initialization; first and second trial return different values.
        -> initialization of jpvt for dgelsy seemed to improve the symptom.
 em_22: Introducing fixvol option
        correct bug in hess_mu calculation in calc_W
 
 todo: 
    Brush up the code 
    switch dgelsy and related function depending on the rank of the matrix (or using try-catch)
*/


/* include header files */
// include header file for matlab
#include "mex.h"
#include "matrix.h"

// include header file for intel mkl
#include "mkl.h"

// include header file for container vector
#include <vector>

// include header for OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// detect memory leak using msvc
//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>
//#include <crtdbg.h>
//#ifdef _DEBUG
//#define   new                   new(_NORMAL_BLOCK, __FILE__, __LINE__)
//#define   malloc(s)             _malloc_dbg(s, _NORMAL_BLOCK, __FILE__, __LINE__)
//#define   calloc(c, s)          _calloc_dbg(c, s, _NORMAL_BLOCK, __FILE__, __LINE__)
//#define   realloc(p, s)         _realloc_dbg(p, s, _NORMAL_BLOCK, __FILE__, __LINE__)
//#define   _recalloc(p, c, s)    _recalloc_dbg(p, c, s, _NORMAL_BLOCK, __FILE__, __LINE__)
//#define   _expand(p, s)         _expand_dbg(p, s, _NORMAL_BLOCK, __FILE__, __LINE__)
//#define   _aligned_malloc(s,m)  _aligned_malloc_dbg(s, m, __FILE__, __LINE__)
//#endif



// include header file for c standard
#ifdef _MSC_VER
    #include <malloc.h>
#else
    #include <stdlib.h>
    static inline void *_aligned_malloc(size_t size, size_t alignment)
    {
        void *p;
        int ret = posix_memalign(&p, alignment, size);
        return (ret == 0) ? p : 0;
    }
    #define _aligned_free(a)    free(a)
#endif

#include <math.h>
#include <stdio.h>
#ifdef _MSC_VER
    #include <intrin.h>
#else
    #include <immintrin.h>
#endif


/* define macro and SIMD settings */
/* parameter settings */
#define NUM_PARAMS_INPUT 10 // 10 pararameters per an ellipsoid for input
#define NUM_PARAMS_PROC 24  // 24 pararameters per an ellipsoid for calc ellipsoid
#define NUM_PARAMS_CALCG 10 // 10 parameters per an ellipsoid for calc_G
#define NUM_PARAMS_SEQ 35 // 35 parameters per an ellipsoid for sequential calculation of V_k
#define NUM_PARAMS_PI_K 1 // 1 parameters per an ellipsoid for pi_k
#define NUM_PARAMS_MU 3 // 3 parameters per an ellipsoid for mu
#define NUM_PARAMS_SIGMA 6 // 6 parameters per an ellipsoid for sigma
#define M_PI 3.14159265358979323846 // math_pi


/* SIMD Vector settings */
//#define ENABLE_SINGLE_PRECISION // uncomment if you use single float for calculation.

#if defined(ENABLE_SINGLE_PRECISION)
    #define FLT_MAX_PRECISION 3.402823466e+38F // max limit of single float
    #define FLT_EXP_MIN_INPUT -87.0 // minimum limit of input for exponential function, single float

    #define SIMD_ELEMENT_TYPE float // single_float
    #define SIMD_VECTOR_LENGTH 8  // AVX is 256-bit SIMD (single_float*8) 
    #define __MM __m256  // AVX is 256-bit SIMD (single_float*8)
    #define __MMI __m256i // AVX is 256-bit SIMD (32bit integer*8)

    #define _MM_MASKLOAD(a,b) _mm256_maskload_ps((a), (_mm256_castps_si256((b))) )
    #define _MM_MOVEMASK(a)   _mm256_movemask_ps((a))
    #define _MM_SET1(a)       _mm256_set1_ps((a))
    #define _MM_CMP(a,b,c)    _mm256_cmp_ps((a), (b), (c))
    #define _MM_AND(a,b)      _mm256_and_ps((a), (b))
    #define _MM_OR(a,b)       _mm256_or_ps((a), (b))
    #define _MM_MAX(a,b)      _mm256_max_ps((a), (b))
    #define _MM_MIN(a,b)      _mm256_min_ps((a), (b))
    #define _MM_ADD(a,b)      _mm256_add_ps((a), (b))
    #define _MM_SUB(a,b)      _mm256_sub_ps((a), (b))
    #define _MM_MUL(a,b)      _mm256_mul_ps((a), (b))
    #define _MM_DIV(a,b)      _mm256_div_ps((a), (b))
    #define _MM_SQR(a)        _mm256_mul_ps((a), (a))
    #define _MM_SQRT(a)       _mm256_sqrt_ps((a))
    #define _MM_EXP(a)        _mm256_exp_ps((a))   // intel compiler intrinsics
    #define _MM_CEIL(a)       _mm256_ceil_ps((a)) // intel compiler intrinsics
    #define _MM_FLOOR(a)      _mm256_floor_ps((a)) // intel compiler intrinsics
    #define _MM_CVTF2I(a)     _mm256_cvtps_epi32((a))
    #define _MM_CVTI2F(a)     _mm256_cvtepi32_ps((a))
    #define _MM_BLENDV(a,b,c) _mm256_blendv_ps((a), (b), (c))
    #define _MM_PERM(a,b)     _mm256_permute_ps((a), (b))    
    #define _MM_PERM2(a,b,c)  _mm256_permute2f128_ps((a), (b), (c))
    #define _MM_TESTZ(a,b)    _mm256_testz_ps( (a), (b) )
    
    typedef union {
    uint32_T uint32;
    float f;
    } aliasing;
    aliasing aliasing_mask_true  = { 0xFFFFFFFF };
    aliasing aliasing_mask_false = { 0x00000000 };
    float mask_true  = aliasing_mask_true.f;
    float mask_false = aliasing_mask_false.f;
    __MM mask_true_mm = _MM_SET1(mask_true);
    __MM mask_false_mm = _MM_SET1(mask_false);
    #define _MM_NOT(a)   _mm256_andnot_ps((a), (mask_true_mm))
    #define _MM_TESTZ1(a) _mm256_testz_ps( (a), (mask_true_mm))
    __MM addx_mm = _mm256_setr_ps(0,1,2,3,4,5,6,7);


#else
    #define FLT_MAX_PRECISION 1.7976931348623158e+308 // max limit of double float
    #define FLT_EXP_MIN_INPUT -708.0 // minimum limit of input for exponential function, double float

    #define SIMD_ELEMENT_TYPE double // double_float
    #define SIMD_VECTOR_LENGTH 4  // AVX is 256-bit SIMD (double_float*4) 
    #define __MM __m256d  // AVX is 256-bit SIMD (double_float*4)
    #define __MMI __m128i // SSE2 is 128-bit SIMD (32bit integer*4)

    #define _MM_MASKLOAD(a,b) _mm256_maskload_pd((a), (_mm256_castpd_si256((b))) )
    #define _MM_MOVEMASK(a)   _mm256_movemask_pd((a))
    #define _MM_SET1(a)       _mm256_set1_pd((a))
    #define _MM_CMP(a,b,c)    _mm256_cmp_pd((a), (b), (c))
    #define _MM_AND(a,b)      _mm256_and_pd((a), (b))
    #define _MM_OR(a,b)       _mm256_or_pd((a), (b))
    #define _MM_MAX(a,b)      _mm256_max_pd((a), (b))
    #define _MM_MIN(a,b)      _mm256_min_pd((a), (b))
    #define _MM_ADD(a,b)      _mm256_add_pd((a), (b))
    #define _MM_SUB(a,b)      _mm256_sub_pd((a), (b))
    #define _MM_MUL(a,b)      _mm256_mul_pd((a), (b))
    #define _MM_DIV(a,b)      _mm256_div_pd((a), (b))
    #define _MM_SQR(a)        _mm256_mul_pd((a), (a))
    #define _MM_SQRT(a)       _mm256_sqrt_pd((a))
    #define _MM_EXP(a)        _mm256_exp_pd((a))   // intel compiler intrinsics
    #define _MM_CEIL(a)       _mm256_ceil_pd((a)) // intel compiler intrinsics
    #define _MM_FLOOR(a)      _mm256_floor_pd((a)) // intel compiler intrinsics
    #define _MM_CVTF2I(a)     _mm256_cvtpd_epi32((a))
    #define _MM_CVTI2F(a)     _mm256_cvtepi32_pd((a))
    #define _MM_BLENDV(a,b,c) _mm256_blendv_pd((a), (b), (c))
    #define _MM_PERM(a,b)     _mm256_permute_pd((a), (b))    
    #define _MM_PERM2(a,b,c)  _mm256_permute2f128_pd((a), (b), (c))
    #define _MM_TESTZ(a,b)    _mm256_testz_pd( (a), (b) )

    typedef union {
    uint64_T uint64;
    double f;
    } aliasing;
    aliasing aliasing_mask_true  = { 0xFFFFFFFFFFFFFFFF };
    aliasing aliasing_mask_false = { 0x0000000000000000 };
    double mask_true  = aliasing_mask_true.f;
    double mask_false = aliasing_mask_false.f;
    __MM mask_true_mm  = _MM_SET1(mask_true);
    __MM mask_false_mm = _MM_SET1(mask_false);
    #define _MM_NOT(a)   _mm256_andnot_pd((a), (mask_true_mm))
    #define _MM_TESTZ1(a) _mm256_testz_pd( (a), (mask_true_mm))
    __MM addx_mm = _mm256_setr_pd(0,1,2,3);

#endif

#define MEMORY_ALIGNMENT 64 // alignment (bit) for memory allocation
typedef SIMD_ELEMENT_TYPE (*fvec)[SIMD_VECTOR_LENGTH]; // SIMD Vector
typedef               int (*ivec)[SIMD_VECTOR_LENGTH]; // SIMD Vector

#ifndef _MSC_VER
// macro of max and min for gcc
#define __max(x, y) (((x) > (y)) ? (x) : (y))
#define __min(x, y) (((x) < (y)) ? (x) : (y))
#define _isnan(x) isnan(x)
#define _finite(x) isfinite(x)
#endif


// #define MY_ASSERT 0 // use clock when defined

// #define ENABLE_TWOSUM // use twosum (more precise version of summation)


// declaration of static variables
SIMD_ELEMENT_TYPE *mask = NULL; // buffer for mask
SIMD_ELEMENT_TYPE *im_orig = NULL; // buffer for original image
SIMD_ELEMENT_TYPE *im_synth = NULL; // buffer for synthesized image
SIMD_ELEMENT_TYPE *idx_Zstart = NULL; // buffer for start point of image calculation
SIMD_ELEMENT_TYPE *idx_Zend = NULL; // buffer for end point of image calculation
SIMD_ELEMENT_TYPE *mem_int0 = NULL; // buffer for start point of image calculation
SIMD_ELEMENT_TYPE *mem_fz0 = NULL; // buffer for start point of image calculation
SIMD_ELEMENT_TYPE *mem_int = NULL; // buffer for image calculation
SIMD_ELEMENT_TYPE *mem_fz = NULL; // buffer for image calculation
SIMD_ELEMENT_TYPE *new_mu_x = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_y = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_z = NULL; // 
SIMD_ELEMENT_TYPE *new_S_11 = NULL; // 
SIMD_ELEMENT_TYPE *new_S_12 = NULL; // 
SIMD_ELEMENT_TYPE *new_S_13 = NULL; // 
SIMD_ELEMENT_TYPE *new_S_22 = NULL; // 
SIMD_ELEMENT_TYPE *new_S_23 = NULL; // 
SIMD_ELEMENT_TYPE *new_S_33 = NULL; // 
SIMD_ELEMENT_TYPE *sumrzg = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_x_err = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_y_err = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_z_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_11_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_12_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_13_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_22_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_23_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_33_err = NULL; // 
SIMD_ELEMENT_TYPE *sumrzg_err = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_x_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_y_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_z_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_11_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_12_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_13_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_22_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_23_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_33_tmp = NULL; // 
SIMD_ELEMENT_TYPE *sumrzg_tmp = NULL; // 
SIMD_ELEMENT_TYPE *mem_sxm_1 = NULL; // buffer for gradient calculation
SIMD_ELEMENT_TYPE *mem_sxm_2 = NULL; // buffer for gradient calculation
SIMD_ELEMENT_TYPE *mem_sxm_3 = NULL; // buffer for gradient calculation
SIMD_ELEMENT_TYPE *gky = NULL; // buffer for calculation of g_k * y
SIMD_ELEMENT_TYPE *gky_err = NULL; // buffer for calculation of g_k * y
SIMD_ELEMENT_TYPE *seq = NULL; // buffer for sequential calculation of V_k
SIMD_ELEMENT_TYPE *seq_err = NULL; // buffer for sequential calculation of V_k
SIMD_ELEMENT_TYPE *params = NULL; // buffer for parameters 
SIMD_ELEMENT_TYPE *params_old = NULL; // buffer for parameters, for temporal isolation
SIMD_ELEMENT_TYPE *params_diff = NULL; // buffer for parameters
SIMD_ELEMENT_TYPE *params_init = NULL; // buffer for parameters
SIMD_ELEMENT_TYPE *params_proc = NULL; // buffer for parameters
SIMD_ELEMENT_TYPE *volume = NULL; // volume of ellipsoid, without pi
SIMD_ELEMENT_TYPE *detS_old = NULL; // relevant to volume of ellipsoid, for fixvol
SIMD_ELEMENT_TYPE *params_calcG = NULL; // buffer for variables in calc_G

SIMD_ELEMENT_TYPE thrint, thrdist, tol;

// declaration of static variables for SIMD
__MM *mask_mm     = NULL;
__MM *im_orig_mm  = NULL;
__MM *im_synth_mm  = NULL;
__MM *idx_Zstart_mm = NULL;
__MM *idx_Zend_mm = NULL;
__MM *mem_int0_mm = NULL;
__MM *mem_fz0_mm = NULL;
__MM *mem_int_mm = NULL;
__MM *mem_fz_mm = NULL;
__MM *new_mu_x_mm = NULL;
__MM *new_mu_y_mm = NULL;
__MM *new_mu_z_mm = NULL;
__MM *new_S_11_mm = NULL;
__MM *new_S_12_mm = NULL;
__MM *new_S_13_mm = NULL;
__MM *new_S_22_mm = NULL;
__MM *new_S_23_mm = NULL;
__MM *new_S_33_mm = NULL;
__MM *sumrzg_mm = NULL;
__MM *new_mu_x_err_mm = NULL;
__MM *new_mu_y_err_mm = NULL;
__MM *new_mu_z_err_mm = NULL;
__MM *new_S_11_err_mm = NULL;
__MM *new_S_12_err_mm = NULL;
__MM *new_S_13_err_mm = NULL;
__MM *new_S_22_err_mm = NULL;
__MM *new_S_23_err_mm = NULL;
__MM *new_S_33_err_mm = NULL;
__MM *sumrzg_err_mm = NULL;
__MM *new_mu_x_tmp_mm = NULL;
__MM *new_mu_y_tmp_mm = NULL;
__MM *new_mu_z_tmp_mm = NULL;
__MM *new_S_11_tmp_mm = NULL;
__MM *new_S_12_tmp_mm = NULL;
__MM *new_S_13_tmp_mm = NULL;
__MM *new_S_22_tmp_mm = NULL;
__MM *new_S_23_tmp_mm = NULL;
__MM *new_S_33_tmp_mm = NULL;
__MM *sumrzg_tmp_mm = NULL;
__MM *mem_sxm_1_mm = NULL;
__MM *mem_sxm_2_mm = NULL;
__MM *mem_sxm_3_mm = NULL;
__MM *gky_mm = NULL; 
__MM *gky_err_mm = NULL; 
__MM *seq_mm = NULL; 
__MM *seq_err_mm = NULL; 
__MM *params_proc_mm = NULL;
__MM *rss_current_mm = NULL;
__MM *rss_current_err_mm = NULL;

__MM mtwos_mm = _MM_SET1(-2);
__MM mones_mm = _MM_SET1(-1);
__MM mhalf_mm = _MM_SET1(-0.5);
__MM zeros_mm = _MM_SET1(0);
__MM ones_mm  = _MM_SET1(1);
__MM twos_mm  = _MM_SET1(2);
__MM numz_mm, numzm1_mm;

__int64* verbose_etime = NULL;
double* verbose_score = NULL;
double* grad_pi_k = NULL; // buffer for gradient of pi_k (g_k*y)
double* grad_mu = NULL; // buffer for gradient of mu
double* grad_sigma = NULL; // buffer for gradient of sigma
double* lambda_mu = NULL; // buffer for lambda_mu
double* lambda_mu_sum = NULL; // buffer for summation of lambda_mu
double* hess_pi_k = NULL; // buffer for calculation of pi_k (g_k * g_l)
double* hess_mu = NULL; // buffer for calculation of mu
double* hess_sigma = NULL; // buffer for calculation of sigma

double rss_orig_inv, score, lambda_mu_coeff, /*lambda_mu_diag,*/ lambda_mu_nondiag, lambda_sigma, lambda_pi, penalty;
double lmp_mu; // Levenberg-Marquardt parameter for update mu
double lmp_sigma; // Levenberg-Marquardt parameter for update sigma
double lmp_orig; // Levenberg-Marquardt parameter for original;
double lmd; // Levenberg-Marquardt damping parameter 
    

int *idx_validelip = NULL;
int *idx_minZstart = NULL;
int *idx_maxZend   = NULL;

int maxLengthZ, fixpik, fixmu, fixvol, /*m1, m2,*/ verbose;
// int idx_update; // specify what should be updatated. -1:calc_ofv, 0:imsynth, 1:pi_k, 2:mu, 3:sigma

unsigned int *idx_VXstart = NULL;
unsigned int *idx_Ystart = NULL;
unsigned int *idx_VXend = NULL;
unsigned int *idx_Yend = NULL;

unsigned int numx, numy, numz, padlength, numXWithPad, numVectorX, maxiter, finiter;
unsigned int numelipsoids = 0;
unsigned int lm_maxiter; // Levenberg-Marquardt loop max iteration



unsigned int s2i[3][3][3][3] = // table for subscripts to index
{ 
  {
    { {  0,  1,  2},
      {  1,  3,  4},
      {  2,  4,  5} },
    { {  1,  3,  4},
      {  3,  6,  7},
      {  4,  7,  8} },
    { {  2,  4,  5},
      {  4,  7,  8},
      {  5,  8,  9} } 
  }, 
  { 
    { {  1,  3,  4},
      {  3,  6,  7},
      {  4,  7,  8} },
    { {  3,  6,  7},
      {  6, 10, 11},
      {  7, 11, 12} },
    { {  4,  7,  8},
      {  7, 11, 12},
      {  8, 12, 13} }
  },
  {
    { {  2,  4,  5},
      {  4,  7,  8},
      {  5,  8,  9} },
    { {  4,  7,  8},
      {  7, 11, 12},
      {  8, 12, 13} },
    { {  5,  8,  9},
      {  8, 12, 13},
      {  9, 13, 14} }
  }
}; 




// The table indicates all of the voxels at a specified xy index of SIMD vector are lower than thrint.
// Will be used in mstep.
std::vector< std::vector<unsigned int> > table_mask_xy(0, std::vector<unsigned int>(0)); 

// The table contains indices of valid elipsoids at a specified xy index of SIMD vector. 
// Will be used in mstep.
std::vector< std::vector<unsigned int> > table_elipsoid_xy(0, std::vector<unsigned int>(0)); 

// The table contains indices of valid elipsoids at a specified z index of SIMD vector. 
// Will be used in mstep.
std::vector< std::vector<unsigned int> > table_elipsoid_z(0, std::vector<unsigned int>(0)); 


const mwSize *dims;
mwSize numdims;
mxArray *prhs_mldivide[2] = {NULL,NULL};
mxArray *plhs_mldivide[1] = {NULL};
mxArray *plhs_rcond[1] = {NULL};




/* function prototype */
void setverbose(double *pverbose);
void setresult(double *presult);
void init_params(const mxArray *prhs0);
void setparam();
static void closefun(void);
void free_im();
void free_params();
void initialize(int nrhs, const mxArray *prhs[]);
inline void twosum(__MM *sum_mm, __MM *err_mm, __MM *add_mm);
inline void twosub(__MM *sum_mm, __MM *err_mm, __MM *sub_mm);
inline void twosum(SIMD_ELEMENT_TYPE *sumnum, SIMD_ELEMENT_TYPE *errnum, SIMD_ELEMENT_TYPE *addnum);
inline void twosub(SIMD_ELEMENT_TYPE *sumnum, SIMD_ELEMENT_TYPE *errnum, SIMD_ELEMENT_TYPE *subnum);
inline int hmax(__MM sbj_mm);
inline int hmin(__MM sbj_mm);
void calc_penalty();
void calc_ofv(int idx_update);
void calc_W(int idx_update);
void update_params(int idx_update);
void optimize();
void setsynth2(double *psynth);


inline void twosum(__MM *sum_mm, __MM *err_mm, __MM *add_mm){
    __MM predictor_mm = _MM_ADD(*sum_mm,*add_mm);
    #ifdef ENABLE_TWOSUM
    __MM corrector_mm = _MM_SUB(predictor_mm,*sum_mm);
    *err_mm = _MM_ADD(*err_mm,_MM_ADD(_MM_SUB(*sum_mm,_MM_SUB(predictor_mm,corrector_mm)),
                                      _MM_SUB(*add_mm,corrector_mm)));
    #endif
    *sum_mm = predictor_mm;
}


inline void twosub(__MM *sum_mm, __MM *err_mm, __MM *sub_mm){
    __MM predictor_mm = _MM_SUB(*sum_mm,*sub_mm);
    #ifdef ENABLE_TWOSUM
    __MM corrector_mm = _MM_SUB(predictor_mm,*sum_mm);
    *err_mm = _MM_ADD(*err_mm,_MM_SUB(_MM_SUB(*sum_mm,_MM_SUB(predictor_mm,corrector_mm)),
                                      _MM_ADD(*sub_mm,corrector_mm)));
    #endif
    *sum_mm = predictor_mm;
}


inline void twosum(SIMD_ELEMENT_TYPE *sumnum, SIMD_ELEMENT_TYPE *errnum, SIMD_ELEMENT_TYPE *addnum){
    volatile SIMD_ELEMENT_TYPE predictor = *sumnum + *addnum;
    #ifdef ENABLE_TWOSUM
    volatile SIMD_ELEMENT_TYPE corrector = predictor - *sumnum;
    volatile SIMD_ELEMENT_TYPE tmpnum = predictor - corrector;
    *errnum += *sumnum - tmpnum + (*addnum-corrector);
    #endif
    *sumnum = predictor;
}


inline void twosub(SIMD_ELEMENT_TYPE *sumnum, SIMD_ELEMENT_TYPE *errnum, SIMD_ELEMENT_TYPE *subnum){
    volatile SIMD_ELEMENT_TYPE predictor = *sumnum - *subnum;
    #ifdef ENABLE_TWOSUM
    volatile SIMD_ELEMENT_TYPE corrector = predictor - *sumnum;
    volatile SIMD_ELEMENT_TYPE tmpnum = predictor - corrector;
    *errnum += *sumnum - tmpnum - (*subnum+corrector);
    #endif
    *sumnum = predictor;
}


inline int hmax(__MM sbj_mm) {

    __MMI sbj_mmi =_MM_CVTF2I(sbj_mm);

    int *sbj = (int *)&sbj_mmi;

    #if defined(ENABLE_SINGLE_PRECISION) // horizontal min and max for single float mm256
        return __max(__max(__max(sbj[0],sbj[1]),
                           __max(sbj[2],sbj[3])),
                     __max(__max(sbj[4],sbj[5]),
                           __max(sbj[6],sbj[7])));

    #else // horizontal min and max for double float mm256d
        return __max(__max(sbj[0],sbj[1]),
                     __max(sbj[2],sbj[3]));
    #endif

}


inline int hmin(__MM sbj_mm) {

    __MMI sbj_mmi =_MM_CVTF2I(sbj_mm);

    int *sbj = (int *)&sbj_mmi;

    #if defined(ENABLE_SINGLE_PRECISION) // horizontal min and min for single float mm256
        return __min(__min(__min(sbj[0],sbj[1]),
                           __min(sbj[2],sbj[3])),
                     __min(__min(sbj[4],sbj[5]),
                           __min(sbj[6],sbj[7])));

    #else // horizontal min and min for double float mm256d
        return __min(__min(sbj[0],sbj[1]),
                     __min(sbj[2],sbj[3]));
    #endif

}


void calc_penalty() {

    penalty = 0.0;
    double lambda1, pd, xd1, yd1, zd1, sd11, sd12, sd13, sd22, sd23, sd33;
    for (unsigned int celipsoid1=0; celipsoid1<numelipsoids; celipsoid1++) {
        
        // for penalty of pi_k
        pd = params_diff[NUM_PARAMS_INPUT*celipsoid1 + 0];
        penalty += lambda_pi*pd*pd;

        // for penalty of mu
        //lambda1 = lambda_mu[celipsoid1*numelipsoids+celipsoid1];
        xd1 = params_diff[NUM_PARAMS_INPUT*celipsoid1 + 1];
        yd1 = params_diff[NUM_PARAMS_INPUT*celipsoid1 + 2];
        zd1 = params_diff[NUM_PARAMS_INPUT*celipsoid1 + 3];
        //penalty += lambda1 * (xd1*xd1+yd1*yd1+zd1*zd1);

        double lambda12, xd2, yd2, zd2, xd12, yd12, zd12;
        for (unsigned int celipsoid2=celipsoid1+1; celipsoid2<numelipsoids; celipsoid2++) {
            lambda12 = lambda_mu[celipsoid1*numelipsoids+celipsoid2];
            xd2 = params_diff[NUM_PARAMS_INPUT*celipsoid2 + 1];
            yd2 = params_diff[NUM_PARAMS_INPUT*celipsoid2 + 2];
            zd2 = params_diff[NUM_PARAMS_INPUT*celipsoid2 + 3];
            xd12 = xd1-xd2;
            yd12 = yd1-yd2;
            zd12 = zd1-zd2;
            penalty += 2 * lambda12 * (xd12*xd12+yd12*yd12+zd12*zd12);
        }

        // for penalty of sigma
        sd11 = params_diff[NUM_PARAMS_INPUT*celipsoid1+4];
        sd12 = params_diff[NUM_PARAMS_INPUT*celipsoid1+5];
        sd13 = params_diff[NUM_PARAMS_INPUT*celipsoid1+6];
        sd22 = params_diff[NUM_PARAMS_INPUT*celipsoid1+7];
        sd23 = params_diff[NUM_PARAMS_INPUT*celipsoid1+8];
        sd33 = params_diff[NUM_PARAMS_INPUT*celipsoid1+9];

        penalty += lambda_sigma * (sd11*sd11+sd22*sd22+sd33*sd33+2*(sd12*sd12+sd13*sd13+sd23*sd23));

    }

}


void setresult(double *presult) {

    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
        presult[NUM_PARAMS_INPUT*celipsoid + 0] = params[NUM_PARAMS_INPUT*celipsoid + 0]; // pi
        presult[NUM_PARAMS_INPUT*celipsoid + 1] = params[NUM_PARAMS_INPUT*celipsoid + 1] + 1.0; // xc
        presult[NUM_PARAMS_INPUT*celipsoid + 2] = params[NUM_PARAMS_INPUT*celipsoid + 2] + 1.0; // yc
        presult[NUM_PARAMS_INPUT*celipsoid + 3] = params[NUM_PARAMS_INPUT*celipsoid + 3] + 1.0; // zc
        presult[NUM_PARAMS_INPUT*celipsoid + 4] = params[NUM_PARAMS_INPUT*celipsoid + 4]; // S_11
        presult[NUM_PARAMS_INPUT*celipsoid + 5] = params[NUM_PARAMS_INPUT*celipsoid + 5]; // S_12
        presult[NUM_PARAMS_INPUT*celipsoid + 6] = params[NUM_PARAMS_INPUT*celipsoid + 6]; // S_13
        presult[NUM_PARAMS_INPUT*celipsoid + 7] = params[NUM_PARAMS_INPUT*celipsoid + 7]; // S_22
        presult[NUM_PARAMS_INPUT*celipsoid + 8] = params[NUM_PARAMS_INPUT*celipsoid + 8]; // S_23
        presult[NUM_PARAMS_INPUT*celipsoid + 9] = params[NUM_PARAMS_INPUT*celipsoid + 9]; // S_33
    }

}


void setverbose(double *pverbose) {

    for (unsigned int citer=0; citer<finiter+1; citer++) {
        pverbose[citer] = (double)verbose_etime[citer];
    }
    for (unsigned int citer=0; citer<finiter+1; citer++) {
        pverbose[citer + finiter+1] = verbose_score[citer];
    }

}


void setparam() {    

    SIMD_ELEMENT_TYPE c = __min(2*thrdist*thrdist,-FLT_EXP_MIN_INPUT);

    for (unsigned int cvxy=0; cvxy<table_elipsoid_xy.size(); cvxy++) {
        table_elipsoid_xy[cvxy].clear();
    }

    for (unsigned int counter=0; counter<NUM_PARAMS_INPUT*numelipsoids; counter++) {
        params_diff[counter] = params[counter] - params_init[counter];
    }

    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {

        SIMD_ELEMENT_TYPE pi   = params[NUM_PARAMS_INPUT*celipsoid + 0];
        SIMD_ELEMENT_TYPE xc   = params[NUM_PARAMS_INPUT*celipsoid + 1];
        SIMD_ELEMENT_TYPE yc   = params[NUM_PARAMS_INPUT*celipsoid + 2];
        SIMD_ELEMENT_TYPE zc   = params[NUM_PARAMS_INPUT*celipsoid + 3];
        SIMD_ELEMENT_TYPE S_11 = params[NUM_PARAMS_INPUT*celipsoid + 4];
        SIMD_ELEMENT_TYPE S_12 = params[NUM_PARAMS_INPUT*celipsoid + 5];
        SIMD_ELEMENT_TYPE S_13 = params[NUM_PARAMS_INPUT*celipsoid + 6];
        SIMD_ELEMENT_TYPE S_22 = params[NUM_PARAMS_INPUT*celipsoid + 7];
        SIMD_ELEMENT_TYPE S_23 = params[NUM_PARAMS_INPUT*celipsoid + 8];
        SIMD_ELEMENT_TYPE S_33 = params[NUM_PARAMS_INPUT*celipsoid + 9];


        SIMD_ELEMENT_TYPE detS =    S_11*S_22*S_33 
                                + 2*S_12*S_13*S_23
                                -   S_11*S_23*S_23 
                                -   S_12*S_12*S_33
                                -   S_13*S_22*S_13;

        SIMD_ELEMENT_TYPE detS_inv = 1/detS;

        SIMD_ELEMENT_TYPE Sinv_11 = detS_inv*(S_22*S_33 - S_23*S_23);
        SIMD_ELEMENT_TYPE Sinv_22 = detS_inv*(S_11*S_33 - S_13*S_13);
        SIMD_ELEMENT_TYPE Sinv_33 = detS_inv*(S_11*S_22 - S_12*S_12);
        SIMD_ELEMENT_TYPE Sinv_12 = detS_inv*(S_13*S_23 - S_12*S_33);
        SIMD_ELEMENT_TYPE Sinv_13 = detS_inv*(S_12*S_23 - S_13*S_22);
        SIMD_ELEMENT_TYPE Sinv_23 = detS_inv*(S_12*S_13 - S_11*S_23);


        SIMD_ELEMENT_TYPE B_1 = xc*Sinv_11 + yc*Sinv_12 + zc*Sinv_13;
        SIMD_ELEMENT_TYPE B_2 = xc*Sinv_12 + yc*Sinv_22 + zc*Sinv_23;
        SIMD_ELEMENT_TYPE B_3 = xc*Sinv_13 + yc*Sinv_23 + zc*Sinv_33;
        SIMD_ELEMENT_TYPE C_1 = B_1*xc + B_2*yc + B_3*zc;
        

        SIMD_ELEMENT_TYPE exp2bz = exp(-Sinv_33);

        SIMD_ELEMENT_TYPE Sinv_33_inv = 1/Sinv_33;
        SIMD_ELEMENT_TYPE Sp_11 = Sinv_11 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE Sp_12 = Sinv_12 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE Sp_13 = Sinv_13 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE Sp_22 = Sinv_22 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE Sp_23 = Sinv_23 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE Sp_33 = Sinv_33 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE cp    = c       * Sinv_33_inv;

        SIMD_ELEMENT_TYPE cx2 = Sp_13*Sp_13 - Sp_11;
        SIMD_ELEMENT_TYPE cy2 = Sp_23*Sp_23 - Sp_22;
        SIMD_ELEMENT_TYPE cxy = Sp_13*Sp_23 - Sp_12;
        SIMD_ELEMENT_TYPE comden_inv = cp / ( cxy*cxy - cx2*cy2 );

        SIMD_ELEMENT_TYPE xsqrt = sqrt( cy2 * comden_inv );
        SIMD_ELEMENT_TYPE ysqrt = sqrt( cx2 * comden_inv );
        
        unsigned int vxstart = (unsigned int) (__min(__max( ceil(xc - xsqrt),0),numx-1)/SIMD_VECTOR_LENGTH);
        unsigned int vxend   = (unsigned int) (__min(__max(floor(xc + xsqrt),0),numx-1)/SIMD_VECTOR_LENGTH);
        unsigned int ystart  = (unsigned int) (__min(__max( ceil(yc - ysqrt),0),numy-1));        
        unsigned int yend    = (unsigned int) (__min(__max(floor(yc + ysqrt),0),numy-1));


        for (unsigned int cy=ystart; cy<=yend; cy++) {
            for (unsigned int cvx=vxstart; cvx<=vxend; cvx++) {
                table_elipsoid_xy[cy*numVectorX+cvx].push_back(celipsoid);
            }
        }


        volume[celipsoid] = sqrt(8*M_PI*M_PI*M_PI)*sqrt(detS);

        idx_VXstart[celipsoid] = vxstart;
        idx_VXend  [celipsoid] = vxend;
        idx_Ystart [celipsoid] = ystart;
        idx_Yend   [celipsoid] = yend;

        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  0] = _MM_SET1(pi);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  1] = _MM_SET1(xc);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  2] = _MM_SET1(yc);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  3] = _MM_SET1(zc);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  4] = _MM_SET1(Sinv_11);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  5] = _MM_SET1(Sinv_12);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  6] = _MM_SET1(Sinv_13);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  7] = _MM_SET1(Sinv_22); 
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  8] = _MM_SET1(Sinv_23);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  9] = _MM_SET1(Sinv_33);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 10] = _MM_SET1(-Sp_13); // cx1
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 11] = _MM_SET1(-Sp_23); // cy1
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 12] = _MM_SET1(cx2); 
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 13] = _MM_SET1(cy2);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 14] = _MM_SET1(2*cxy);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 15] = _MM_SET1(cp);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 16] = _MM_SET1(-0.5*Sinv_33); // -bz
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 17] = _MM_SET1(exp2bz);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 18] = _MM_SET1(S_11);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 19] = _MM_SET1(S_12);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 20] = _MM_SET1(S_13);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 21] = _MM_SET1(S_22);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 22] = _MM_SET1(S_23);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 23] = _MM_SET1(S_33);

        params_calcG[NUM_PARAMS_CALCG*celipsoid + 0] = Sinv_11; // A_11
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 1] = Sinv_12; // A_12
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 2] = Sinv_13; // A_13
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 3] = Sinv_22; // A_22
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 4] = Sinv_23; // A_23
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 5] = Sinv_33; // A_33
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 6] = B_1; // B_1
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 7] = B_2; // B_2
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 8] = B_3; // B_3
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 9] = C_1; // C
        

    }

}


void init_params(const mxArray *prhs0){

	/* get the number of parameters and parameter sets */
	unsigned int m = (unsigned int)mxGetM(prhs0);
	unsigned int n = (unsigned int)mxGetN(prhs0);

    int maxnumthreads = 1;
    #ifdef _OPENMP
        maxnumthreads = omp_get_max_threads();
    #endif

	if ( (m*n) != NUM_PARAMS_INPUT*numelipsoids) {

        numelipsoids = (m*n) / NUM_PARAMS_INPUT;
        
        free_params(); // free memory for params 

        idx_validelip = (int *)_aligned_malloc(sizeof(int)*numelipsoids*maxnumthreads,MEMORY_ALIGNMENT);
        idx_minZstart = (int *)_aligned_malloc(sizeof(int)*numelipsoids,MEMORY_ALIGNMENT);
        idx_maxZend   = (int *)_aligned_malloc(sizeof(int)*numelipsoids,MEMORY_ALIGNMENT);


        idx_VXstart = (unsigned int *)_aligned_malloc(sizeof(unsigned int)*numelipsoids,MEMORY_ALIGNMENT);
        idx_Ystart  = (unsigned int *)_aligned_malloc(sizeof(unsigned int)*numelipsoids,MEMORY_ALIGNMENT);
        idx_VXend   = (unsigned int *)_aligned_malloc(sizeof(unsigned int)*numelipsoids,MEMORY_ALIGNMENT);
        idx_Yend    = (unsigned int *)_aligned_malloc(sizeof(unsigned int)*numelipsoids,MEMORY_ALIGNMENT);
        
        size_t tmpsize = sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*numelipsoids*maxnumthreads;
        idx_Zstart   = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        idx_Zend     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_int0     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_fz0      = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_int      = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_fz       = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_x     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_y     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_z     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_11     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_12     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_13     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_22     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_23     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_33     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        sumrzg       = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_x_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_y_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_z_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_11_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_12_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_13_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_22_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_23_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_33_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        sumrzg_err   = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_x_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_y_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_z_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_11_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_12_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_13_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_22_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_23_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_33_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        sumrzg_tmp   = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_sxm_1    = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_sxm_2    = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_sxm_3    = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        gky          = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        gky_err      = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        seq          = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize*NUM_PARAMS_SEQ,MEMORY_ALIGNMENT);
        seq_err      = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize*NUM_PARAMS_SEQ,MEMORY_ALIGNMENT);

        params_proc = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*numelipsoids*NUM_PARAMS_PROC,MEMORY_ALIGNMENT);
        params_init = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numelipsoids*NUM_PARAMS_INPUT,MEMORY_ALIGNMENT);
        params_diff = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numelipsoids*NUM_PARAMS_INPUT,MEMORY_ALIGNMENT);
        params_old  = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numelipsoids*NUM_PARAMS_INPUT,MEMORY_ALIGNMENT);
        params      = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numelipsoids*NUM_PARAMS_INPUT,MEMORY_ALIGNMENT);
        volume      = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numelipsoids,MEMORY_ALIGNMENT);
        detS_old    = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numelipsoids,MEMORY_ALIGNMENT);
        params_calcG = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numelipsoids*NUM_PARAMS_CALCG,MEMORY_ALIGNMENT);

        grad_pi_k  = (double *)_aligned_malloc(sizeof(double)*numelipsoids*NUM_PARAMS_PI_K,MEMORY_ALIGNMENT);
        grad_mu    = (double *)_aligned_malloc(sizeof(double)*numelipsoids*NUM_PARAMS_MU,MEMORY_ALIGNMENT);
        grad_sigma = (double *)_aligned_malloc(sizeof(double)*numelipsoids*NUM_PARAMS_SIGMA,MEMORY_ALIGNMENT);        
        lambda_mu  = (double *)_aligned_malloc(sizeof(double)*numelipsoids*numelipsoids,MEMORY_ALIGNMENT);
        lambda_mu_sum  = (double *)_aligned_malloc(sizeof(double)*numelipsoids,MEMORY_ALIGNMENT);
        hess_pi_k  = (double *)_aligned_malloc(sizeof(double)*numelipsoids*numelipsoids*NUM_PARAMS_PI_K*NUM_PARAMS_PI_K,MEMORY_ALIGNMENT);
        hess_mu    = (double *)_aligned_malloc(sizeof(double)*numelipsoids*numelipsoids*NUM_PARAMS_MU*NUM_PARAMS_MU,MEMORY_ALIGNMENT);
        hess_sigma = (double *)_aligned_malloc(sizeof(double)*numelipsoids*numelipsoids*NUM_PARAMS_SIGMA*NUM_PARAMS_SIGMA,MEMORY_ALIGNMENT);
        

        /* convert to SIMD array */
        idx_Zstart_mm = (__MM *) idx_Zstart;
        idx_Zend_mm = (__MM *) idx_Zend;
        mem_int0_mm = (__MM *) mem_int0;
        mem_fz0_mm  = (__MM *) mem_fz0;
        mem_int_mm = (__MM *) mem_int;
        mem_fz_mm  = (__MM *) mem_fz;
        new_mu_x_mm = (__MM *) new_mu_x;
        new_mu_y_mm = (__MM *) new_mu_y;
        new_mu_z_mm = (__MM *) new_mu_z;
        new_S_11_mm = (__MM *) new_S_11;
        new_S_12_mm = (__MM *) new_S_12;
        new_S_13_mm = (__MM *) new_S_13;
        new_S_22_mm = (__MM *) new_S_22;
        new_S_23_mm = (__MM *) new_S_23;
        new_S_33_mm = (__MM *) new_S_33;
        sumrzg_mm = (__MM *) sumrzg;
        new_mu_x_err_mm = (__MM *) new_mu_x_err;
        new_mu_y_err_mm = (__MM *) new_mu_y_err;
        new_mu_z_err_mm = (__MM *) new_mu_z_err;
        new_S_11_err_mm = (__MM *) new_S_11_err;
        new_S_12_err_mm = (__MM *) new_S_12_err;
        new_S_13_err_mm = (__MM *) new_S_13_err;
        new_S_22_err_mm = (__MM *) new_S_22_err;
        new_S_23_err_mm = (__MM *) new_S_23_err;
        new_S_33_err_mm = (__MM *) new_S_33_err;
        sumrzg_err_mm = (__MM *) sumrzg_err;
        new_mu_x_tmp_mm = (__MM *) new_mu_x_tmp;
        new_mu_y_tmp_mm = (__MM *) new_mu_y_tmp;
        new_mu_z_tmp_mm = (__MM *) new_mu_z_tmp;
        new_S_11_tmp_mm = (__MM *) new_S_11_tmp;
        new_S_12_tmp_mm = (__MM *) new_S_12_tmp;
        new_S_13_tmp_mm = (__MM *) new_S_13_tmp;
        new_S_22_tmp_mm = (__MM *) new_S_22_tmp;
        new_S_23_tmp_mm = (__MM *) new_S_23_tmp;
        new_S_33_tmp_mm = (__MM *) new_S_33_tmp;
        sumrzg_tmp_mm = (__MM *) sumrzg_tmp;
        mem_sxm_1_mm = (__MM *) mem_sxm_1;
        mem_sxm_2_mm = (__MM *) mem_sxm_2;
        mem_sxm_3_mm = (__MM *) mem_sxm_3;
        gky_mm = (__MM *) gky;
        gky_err_mm = (__MM *) gky_err;
        seq_mm = (__MM *) seq;
        seq_err_mm = (__MM *) seq_err;
        params_proc_mm = (__MM *) params_proc;
        

	}

	double *pin = mxGetPr(prhs0);

    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {

        params_init[NUM_PARAMS_INPUT*celipsoid + 0] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 0]; // pi
        params_init[NUM_PARAMS_INPUT*celipsoid + 1] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 1] - 1; // xc, 0-start
        params_init[NUM_PARAMS_INPUT*celipsoid + 2] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 2] - 1; // yc, 0-start
        params_init[NUM_PARAMS_INPUT*celipsoid + 3] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 3] - 1; // zc, 0-start
        params_init[NUM_PARAMS_INPUT*celipsoid + 4] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 4]; // S_11
        params_init[NUM_PARAMS_INPUT*celipsoid + 5] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 5]; // S_12
        params_init[NUM_PARAMS_INPUT*celipsoid + 6] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 6]; // S_13
        params_init[NUM_PARAMS_INPUT*celipsoid + 7] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 7]; // S_22
        params_init[NUM_PARAMS_INPUT*celipsoid + 8] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 8]; // S_23
        params_init[NUM_PARAMS_INPUT*celipsoid + 9] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 9]; // S_33

        params[NUM_PARAMS_INPUT*celipsoid + 0] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 0]; // pi
        params[NUM_PARAMS_INPUT*celipsoid + 1] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 1] - 1; // xc, 0-start
        params[NUM_PARAMS_INPUT*celipsoid + 2] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 2] - 1; // yc, 0-start
        params[NUM_PARAMS_INPUT*celipsoid + 3] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 3] - 1; // zc, 0-start
        params[NUM_PARAMS_INPUT*celipsoid + 4] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 4]; // S_11
        params[NUM_PARAMS_INPUT*celipsoid + 5] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 5]; // S_12
        params[NUM_PARAMS_INPUT*celipsoid + 6] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 6]; // S_13
        params[NUM_PARAMS_INPUT*celipsoid + 7] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 7]; // S_22
        params[NUM_PARAMS_INPUT*celipsoid + 8] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 8]; // S_23
        params[NUM_PARAMS_INPUT*celipsoid + 9] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 9]; // S_33

        // obtain detS_old
        SIMD_ELEMENT_TYPE S_11 = params[NUM_PARAMS_INPUT*celipsoid + 4];
        SIMD_ELEMENT_TYPE S_12 = params[NUM_PARAMS_INPUT*celipsoid + 5];
        SIMD_ELEMENT_TYPE S_13 = params[NUM_PARAMS_INPUT*celipsoid + 6];
        SIMD_ELEMENT_TYPE S_22 = params[NUM_PARAMS_INPUT*celipsoid + 7];
        SIMD_ELEMENT_TYPE S_23 = params[NUM_PARAMS_INPUT*celipsoid + 8];
        SIMD_ELEMENT_TYPE S_33 = params[NUM_PARAMS_INPUT*celipsoid + 9];


        SIMD_ELEMENT_TYPE detS =    S_11*S_22*S_33 
                                + 2*S_12*S_13*S_23
                                -   S_11*S_23*S_23 
                                -   S_12*S_12*S_33
                                -   S_13*S_22*S_13;
        detS_old[celipsoid] = detS;

    }

    // obtain lambda_mu
    for (unsigned int celipsoid1=0; celipsoid1<numelipsoids; celipsoid1++) {
        //lambda_mu[celipsoid1*numelipsoids+celipsoid1] = lambda_mu_diag;
        lambda_mu[celipsoid1*numelipsoids+celipsoid1] = 0;

        SIMD_ELEMENT_TYPE mu1_x = params[NUM_PARAMS_INPUT*celipsoid1 + 1]; // mu_i_x
        SIMD_ELEMENT_TYPE mu1_y = params[NUM_PARAMS_INPUT*celipsoid1 + 2]; // mu_i_y
        SIMD_ELEMENT_TYPE mu1_z = params[NUM_PARAMS_INPUT*celipsoid1 + 3]; // mu_i_z
        
        for (unsigned int celipsoid2=celipsoid1+1; celipsoid2<numelipsoids; celipsoid2++) {
            SIMD_ELEMENT_TYPE mu2_x = params[NUM_PARAMS_INPUT*celipsoid2 + 1]; // mu_j_x
            SIMD_ELEMENT_TYPE mu2_y = params[NUM_PARAMS_INPUT*celipsoid2 + 2]; // mu_j_y
            SIMD_ELEMENT_TYPE mu2_z = params[NUM_PARAMS_INPUT*celipsoid2 + 3]; // mu_j_z
            
            SIMD_ELEMENT_TYPE mudist_x = (mu1_x-mu2_x)*(mu1_x-mu2_x);
            SIMD_ELEMENT_TYPE mudist_y = (mu1_y-mu2_y)*(mu1_y-mu2_y);
            SIMD_ELEMENT_TYPE mudist_z = (mu1_z-mu2_z)*(mu1_z-mu2_z);
            
            // SIMD_ELEMENT_TYPE distance = pow(sqrt(mudist_x+mudist_y+mudist_z),-lambda_mu_coeff);
            SIMD_ELEMENT_TYPE distance = exp(-(mudist_x+mudist_y+mudist_z)/(lambda_mu_coeff*lambda_mu_coeff));
            
            lambda_mu[numelipsoids*celipsoid1+celipsoid2] = lambda_mu_nondiag*distance;
            lambda_mu[numelipsoids*celipsoid2+celipsoid1] = lambda_mu_nondiag*distance;
        }
    }

    for (unsigned int celipsoid1=0; celipsoid1<numelipsoids; celipsoid1++) {
        lambda_mu_sum[celipsoid1] = 0.0;
        for (unsigned int celipsoid2=celipsoid1+1; celipsoid2<numelipsoids; celipsoid2++) {
        //lambda_mu[celipsoid1*numelipsoids+celipsoid1] = lambda_mu_diag;
        lambda_mu_sum[celipsoid1] += lambda_mu[numelipsoids*celipsoid1+celipsoid2];
        }
    }



    // Initialize Levenberg-Marquardt parameter
    lmp_mu = lmp_orig;
    lmp_sigma = lmp_orig;    

}


void free_im(){

    /* free memory for global variables */
    
    _aligned_free(mask);
    _aligned_free(im_orig);
    _aligned_free(im_synth);
    _aligned_free(verbose_etime);
    _aligned_free(verbose_score);
    //_aligned_free(im_init_int);
    //_aligned_free(im_init_fz);
    //_aligned_free(im_end_int);
    //_aligned_free(flag_invalid_allZ);
    _aligned_free(rss_current_mm);
    _aligned_free(rss_current_err_mm);
    
	/* assign NULL for pointers of global variables */

    mask = NULL;
    im_orig = NULL;
    im_synth = NULL;
    verbose_etime = NULL;
    verbose_score = NULL;
    //im_init_int = NULL;
    //im_init_fz = NULL;
    //im_end_int = NULL;
    //flag_invalid_allZ = NULL;
    
    mask_mm     = NULL;
    im_orig_mm  = NULL;
    im_synth_mm  = NULL;
    //imInitInt_mm = NULL;
    //imInitFz_mm  = NULL;
    //imEndInt_mm  = NULL;
    rss_current_mm = NULL;
    rss_current_err_mm = NULL;
   
    ///* detect memory leaks using msvc */
    //#ifdef _CRTDBG_MAP_ALLOC
    //_CrtDumpMemoryLeaks();
    //#endif

}


void free_params() {

	/* free memory for global variables */
    _aligned_free(idx_validelip);
    _aligned_free(idx_minZstart);
    _aligned_free(idx_maxZend);

    _aligned_free(idx_VXstart);
    _aligned_free(idx_Ystart);
    _aligned_free(idx_VXend);
    _aligned_free(idx_Yend);
    
    _aligned_free(idx_Zstart);
    _aligned_free(idx_Zend);
    _aligned_free(mem_int0);
    _aligned_free(mem_fz0);
    _aligned_free(mem_int);
    _aligned_free(mem_fz);
    _aligned_free(new_mu_x);
    _aligned_free(new_mu_y);
    _aligned_free(new_mu_z);
    _aligned_free(new_S_11);
    _aligned_free(new_S_12);
    _aligned_free(new_S_13);
    _aligned_free(new_S_22);
    _aligned_free(new_S_23);
    _aligned_free(new_S_33);
    _aligned_free(sumrzg);
    _aligned_free(new_mu_x_err);
    _aligned_free(new_mu_y_err);
    _aligned_free(new_mu_z_err);
    _aligned_free(new_S_11_err);
    _aligned_free(new_S_12_err);
    _aligned_free(new_S_13_err);
    _aligned_free(new_S_22_err);
    _aligned_free(new_S_23_err);
    _aligned_free(new_S_33_err);
    _aligned_free(sumrzg_err);
    _aligned_free(new_mu_x_tmp);
    _aligned_free(new_mu_y_tmp);
    _aligned_free(new_mu_z_tmp);
    _aligned_free(new_S_11_tmp);
    _aligned_free(new_S_12_tmp);
    _aligned_free(new_S_13_tmp);
    _aligned_free(new_S_22_tmp);
    _aligned_free(new_S_23_tmp);
    _aligned_free(new_S_33_tmp);
    _aligned_free(sumrzg_tmp);
    _aligned_free(mem_sxm_1);
    _aligned_free(mem_sxm_2);
    _aligned_free(mem_sxm_3);
    _aligned_free(gky);
    _aligned_free(gky_err);
    _aligned_free(seq);
    _aligned_free(seq_err);
    _aligned_free(params_proc);
    _aligned_free(params_init);
    _aligned_free(params_diff);
    _aligned_free(params_old);
    _aligned_free(params);
    _aligned_free(volume);
    _aligned_free(detS_old);
    _aligned_free(params_calcG);
    _aligned_free(grad_pi_k);
    _aligned_free(grad_mu);
    _aligned_free(grad_sigma);
    _aligned_free(lambda_mu);
    _aligned_free(lambda_mu_sum);
    _aligned_free(hess_pi_k);
    _aligned_free(hess_mu);
    _aligned_free(hess_sigma);


	/* assign NULL for pointers of global variables */
    idx_validelip = NULL;
    idx_minZstart = NULL;
    idx_maxZend = NULL;

    idx_VXstart = NULL;
    idx_Ystart = NULL;
    idx_VXend = NULL;
    idx_Yend = NULL;
    
    idx_Zstart = NULL;
    idx_Zend = NULL;
    mem_int0 = NULL;
    mem_fz0 = NULL;
    mem_int = NULL;
    mem_fz = NULL;
    new_mu_x = NULL;
    new_mu_y = NULL;
    new_mu_z = NULL;
    new_S_11 = NULL;
    new_S_12 = NULL;
    new_S_13 = NULL;
    new_S_22 = NULL;
    new_S_23 = NULL;
    new_S_33 = NULL;
    sumrzg = NULL;
    new_mu_x_err = NULL;
    new_mu_y_err = NULL;
    new_mu_z_err = NULL;
    new_S_11_err = NULL;
    new_S_12_err = NULL;
    new_S_13_err = NULL;
    new_S_22_err = NULL;
    new_S_23_err = NULL;
    new_S_33_err = NULL;
    sumrzg_err = NULL;
    new_mu_x_tmp = NULL;
    new_mu_y_tmp = NULL;
    new_mu_z_tmp = NULL;
    new_S_11_tmp = NULL;
    new_S_12_tmp = NULL;
    new_S_13_tmp = NULL;
    new_S_22_tmp = NULL;
    new_S_23_tmp = NULL;
    new_S_33_tmp = NULL;
    sumrzg_tmp = NULL;
    mem_sxm_1 = NULL;
    mem_sxm_2 = NULL;
    mem_sxm_3 = NULL;
    gky = NULL;
    gky_err = NULL;
    seq = NULL;
    seq_err = NULL;
    params_proc = NULL;
    params_init = NULL;
    params_diff = NULL;
    params_old = NULL;
    params = NULL;
    volume = NULL;
    detS_old = NULL;
    params_calcG = NULL;
    grad_pi_k = NULL;
    grad_mu = NULL;
    grad_sigma = NULL;
    lambda_mu = NULL;
    lambda_mu_sum = NULL;
    hess_pi_k = NULL;
    hess_mu = NULL;
    hess_sigma = NULL;

    /* assign NULL for pointers of SIMD array */
    idx_Zstart_mm = NULL;
    idx_Zend_mm = NULL;
    mem_int0_mm = NULL;
    mem_fz0_mm = NULL;
    mem_int_mm = NULL;
    mem_fz_mm = NULL;
    new_mu_x_mm = NULL;
    new_mu_y_mm = NULL;
    new_mu_z_mm = NULL;
    new_S_11_mm = NULL;
    new_S_12_mm = NULL;
    new_S_13_mm = NULL;
    new_S_22_mm = NULL;
    new_S_23_mm = NULL;
    new_S_33_mm = NULL;
    sumrzg_mm = NULL;
    new_mu_x_err_mm = NULL;
    new_mu_y_err_mm = NULL;
    new_mu_z_err_mm = NULL;
    new_S_11_err_mm = NULL;
    new_S_12_err_mm = NULL;
    new_S_13_err_mm = NULL;
    new_S_22_err_mm = NULL;
    new_S_23_err_mm = NULL;
    new_S_33_err_mm = NULL;
    sumrzg_err_mm = NULL;
    new_mu_x_tmp_mm = NULL;
    new_mu_y_tmp_mm = NULL;
    new_mu_z_tmp_mm = NULL;
    new_S_11_tmp_mm = NULL;
    new_S_12_tmp_mm = NULL;
    new_S_13_tmp_mm = NULL;
    new_S_22_tmp_mm = NULL;
    new_S_23_tmp_mm = NULL;
    new_S_33_tmp_mm = NULL;
    sumrzg_tmp_mm = NULL;
    mem_sxm_1_mm = NULL;
    mem_sxm_2_mm = NULL;
    mem_sxm_3_mm = NULL;
    gky_mm = NULL;
    gky_err_mm = NULL;
    seq_mm = NULL;
    seq_err_mm = NULL;
    params_proc_mm = NULL;
    
    ///* detect memory leaks using msvc */
    //#ifdef _CRTDBG_MAP_ALLOC
    //_CrtDumpMemoryLeaks();
    //#endif

}


static void closefun(void) {
    free_im();
    free_params();

    ///* detect memory leaks using msvc */
    //#ifdef _CRTDBG_MAP_ALLOC
    //_CrtDumpMemoryLeaks();
    //#endif
}


void initialize(int nrhs, const mxArray *prhs[]) {
	/* declaration of variables */
	
    /* register memory freeing function that is called at finishing of the mex function */    
	mexAtExit(closefun);
	
    /* free memory if already allocated */
	free_im();		

    /* processing option settings */

    // set default values
    thrint = -FLT_MAX_PRECISION;
    thrdist = sqrt(-FLT_EXP_MIN_INPUT);
    tol = 5e-5;
    maxiter = 1000;
    fixpik = 0;
    fixmu = 0;
    fixvol = 0;
    //m1 = 0;
    //m2 = 0;
    verbose = 0;
    lambda_pi = 0.0;
    lambda_mu_coeff = 1.0;
    //lambda_mu_diag    = 0.0;
    lambda_mu_nondiag = 0.0;
    lambda_sigma = 0.0;
    lmp_orig = 1;
    lmd = 2;
    lm_maxiter = 100;

    // set the specified values
    if (nrhs==3 && mxIsStruct(prhs[2])) {
        mxArray *tmparr;
        if (mxGetFieldNumber(prhs[2],"thrint")>=0) {
            tmparr = mxGetField(prhs[2],0,"thrint");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                thrint = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"thrdist")>=0) {
            tmparr = mxGetField(prhs[2],0,"thrdist");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                thrdist = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"tol")>=0) {
            tmparr = mxGetField(prhs[2],0,"tol");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                tol = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"maxiter")>=0) {
            tmparr = mxGetField(prhs[2],0,"maxiter");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                maxiter = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"fixpik")>=0) {
            tmparr = mxGetField(prhs[2],0,"fixpik");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                fixpik = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"fixmu")>=0) {
            tmparr = mxGetField(prhs[2],0,"fixmu");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                fixmu = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"fixvol")>=0) {
            tmparr = mxGetField(prhs[2],0,"fixvol");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                fixvol = mxGetScalar(tmparr);
            }
        }
        //if (mxGetFieldNumber(prhs[2],"m1")>=0) {
        //    tmparr = mxGetField(prhs[2],0,"m1");
        //    if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
        //        m1 = mxGetScalar(tmparr);
        //    }
        //}
        //if (mxGetFieldNumber(prhs[2],"m2")>=0) {
        //    tmparr = mxGetField(prhs[2],0,"m2");
        //    if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
        //        m2 = mxGetScalar(tmparr);
        //    }
        //}
        if (mxGetFieldNumber(prhs[2],"verbose")>=0) {
            tmparr = mxGetField(prhs[2],0,"verbose");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                verbose = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"lambda_pi")>=0) {
            tmparr = mxGetField(prhs[2],0,"lambda_pi");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                lambda_pi = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"lambda_mu_coeff")>=0) {
            tmparr = mxGetField(prhs[2],0,"lambda_mu_coeff");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                lambda_mu_coeff = mxGetScalar(tmparr);
            }
        }
        //if (mxGetFieldNumber(prhs[2],"lambda_mu_diag")>=0) {
        //    tmparr = mxGetField(prhs[2],0,"lambda_mu_diag");
        //    if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
        //        lambda_mu_diag = mxGetScalar(tmparr);
        //    }
        //}
        if (mxGetFieldNumber(prhs[2],"lambda_mu_nondiag")>=0) {
            tmparr = mxGetField(prhs[2],0,"lambda_mu_nondiag");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                lambda_mu_nondiag = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"lambda_sigma")>=0) {
            tmparr = mxGetField(prhs[2],0,"lambda_sigma");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                lambda_sigma = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"lm_init")>=0) {
            tmparr = mxGetField(prhs[2],0,"lm_init");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                lmp_orig = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"lm_damp")>=0) {
            tmparr = mxGetField(prhs[2],0,"lm_damp");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                lmd = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"lm_maxiter")>=0) {
            tmparr = mxGetField(prhs[2],0,"lm_maxiter");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                lm_maxiter = mxGetScalar(tmparr);
            }
        }
    }
	
	/* get the size of the image: [numx,numy,numz] = size(im)*/
	dims = mxGetDimensions(prhs[1]);
	numdims = mxGetNumberOfDimensions(prhs[1]);
	numx = (unsigned int)dims[0]; // size(im,1)
	numy = (unsigned int)dims[1]; // size(im,2)
	if ((unsigned int) numdims >= 3) {
		numz = (unsigned int)dims[2]; // size(im,3)
	} else {
		numz = 1;  // size(im,3)
	}
    numz_mm = _MM_SET1(numz);
    numzm1_mm = _MM_SET1(numz-1);
	
    padlength = (SIMD_VECTOR_LENGTH - (numx % SIMD_VECTOR_LENGTH))%SIMD_VECTOR_LENGTH;
    numXWithPad    = numx + padlength;   
    numVectorX = numXWithPad / SIMD_VECTOR_LENGTH;
    
    unsigned int numElementXYZ = numXWithPad * numy * numz;

    int maxnumthreads = 1;
    #ifdef _OPENMP
        maxnumthreads = omp_get_max_threads();
    #endif

    /* memory allocation */
    mask        = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numElementXYZ,MEMORY_ALIGNMENT);
    im_orig     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numElementXYZ,MEMORY_ALIGNMENT);
    im_synth    = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numElementXYZ,MEMORY_ALIGNMENT);
    if (verbose==1) {
        verbose_etime = (__int64 *)_aligned_malloc(sizeof(__int64)*maxiter,MEMORY_ALIGNMENT);
        verbose_score = ( double *)_aligned_malloc(sizeof( double)*maxiter,MEMORY_ALIGNMENT);
    }
    //im_init_int = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*(numz+2),MEMORY_ALIGNMENT);
    //im_init_fz  = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*(numz+2),MEMORY_ALIGNMENT);
    //im_end_int  = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*(numz+2),MEMORY_ALIGNMENT);
    //flag_invalid_allZ     =       (int *)_aligned_malloc(sizeof(int)*numVectorX*numy,MEMORY_ALIGNMENT);
    rss_current_mm     = (__MM *) _aligned_malloc(sizeof(__MM)*maxnumthreads,MEMORY_ALIGNMENT);
    rss_current_err_mm = (__MM *) _aligned_malloc(sizeof(__MM)*maxnumthreads,MEMORY_ALIGNMENT);


    /* convert to SIMD-vector array */
    mask_mm     = (__MM *)mask;
    im_orig_mm  = (__MM *)im_orig;
    im_synth_mm  = (__MM *)im_synth;
    //imInitInt_mm = (__MM *)im_init_int;
    //imInitFz_mm  = (__MM *)im_init_fz;
    //imEndInt_mm  = (__MM *)im_end_int;
    

    ///* initialize im_init */
    //for (unsigned int c=0; c<numz+2; c++) {
    //    imInitInt_mm[c] = zeros_mm;
    //    imInitFz_mm[c]  = zeros_mm;
    //    imEndInt_mm[c]  = mask_true_mm;
    //}

    /* initialize table_mask_xy */
    table_mask_xy.resize(numy);
    for (unsigned int c=0; c<table_mask_xy.size(); c++) {
        table_mask_xy[c].clear();
        table_mask_xy[c].reserve(numVectorX/4);
    }
    
    /* initialize table_elipsoid_xy */
    table_elipsoid_xy.resize(numVectorX*numy);
    for (unsigned int c=0; c<table_elipsoid_xy.size(); c++) {
        table_elipsoid_xy[c].resize(1);
    }

    /* initialize table_elipsoid_z */
    table_elipsoid_z.resize(numz);
    for (unsigned int c=0; c<table_elipsoid_z.size(); c++) {
        table_elipsoid_z[c].resize(1);
    }

    /* get pointer of real part of input */
	double *pim;
    pim = mxGetPr(prhs[1]); // im


    /* copy original image data to im_orig */        
    unsigned int cvxyz, cin;
    __MM tmpsqr_mm, mask_padding_mm;
    __MM thrint_mm = _MM_SET1( (SIMD_ELEMENT_TYPE) thrint );
    __MM rss_orig_mm = zeros_mm;
    __MM rss_orig_err_mm = zeros_mm;       
    
    // make mask for padding
    SIMD_ELEMENT_TYPE *mask_padding = (SIMD_ELEMENT_TYPE *)&mask_padding_mm;
    for (unsigned int ce=0; ce<SIMD_VECTOR_LENGTH-padlength; ce++) {
        mask_padding[ce] = mask_true;
    }
    for (unsigned int ce=SIMD_VECTOR_LENGTH-padlength; ce<SIMD_VECTOR_LENGTH; ce++) {
        mask_padding[ce] = mask_false;
    }

    
    for (unsigned int cy=0; cy<numy; cy++) {        
        
        // for elements without padding        
        for (unsigned int cvx=0; cvx<numVectorX-1; cvx++) {
            __MM tmpmask;
            //__MM tmpmask_sum = mask_false_mm;
            for (unsigned int cz=0; cz<numz; cz++) {
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                cin   = numx*numy*cz + numx*cy + SIMD_VECTOR_LENGTH*cvx;
                for (unsigned int ce=0; ce<SIMD_VECTOR_LENGTH; ce++) {
                    im_orig[ cvxyz*SIMD_VECTOR_LENGTH + ce ] = (SIMD_ELEMENT_TYPE) pim[ cin + ce ];                          
                }
                tmpmask = _MM_CMP(im_orig_mm[cvxyz],thrint_mm,_CMP_GE_OS);
                im_orig_mm[cvxyz] = _MM_AND(im_orig_mm[cvxyz],tmpmask);
                //mask_mm[cvxyz] = tmpmask;
                mask_mm[cvxyz] = mask_true_mm;
                //tmpmask_sum = _MM_OR(tmpmask_sum,tmpmask);
                tmpsqr_mm = _MM_SQR(_MM_AND(tmpmask,im_orig_mm[cvxyz]));
                twosum(&rss_orig_mm,&rss_orig_err_mm,&tmpsqr_mm);
            }
            //flag_invalid_allZ[cy*numVectorX+cvx] = _MM_TESTZ1(tmpmask_sum);
            //if (!_MM_TESTZ1(tmpmask_sum)) {
                table_mask_xy[cy].push_back(cvx);
            //}
            
        }

        // for elements with padding
        unsigned int cvx = numVectorX - 1; {
            __MM tmpmask;
            //__MM tmpmask_sum = mask_false_mm;
            for (unsigned int cz=0; cz<numz; cz++) {
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                cin   = numx*numy*cz + numx*cy + SIMD_VECTOR_LENGTH*cvx;
                for (unsigned int ce=0; ce<SIMD_VECTOR_LENGTH-padlength; ce++) {
                    im_orig[ cvxyz*SIMD_VECTOR_LENGTH + ce ] = (SIMD_ELEMENT_TYPE) pim[ cin + ce ];                          
                }
                for (unsigned int ce=SIMD_VECTOR_LENGTH-padlength; ce<SIMD_VECTOR_LENGTH; ce++) {
                    im_orig[ cvxyz*SIMD_VECTOR_LENGTH + ce ] = 0; // padding with 0
                }            
                tmpmask = _MM_AND(mask_padding_mm,_MM_CMP(im_orig_mm[cvxyz],thrint_mm,_CMP_GE_OS));
                im_orig_mm[cvxyz] = _MM_AND(im_orig_mm[cvxyz],tmpmask);
                //mask_mm[cvxyz] = tmpmask;
                mask_mm[cvxyz] = mask_padding_mm;
                //tmpmask_sum = _MM_OR(tmpmask_sum,tmpmask);
                tmpsqr_mm = _MM_SQR(_MM_AND(tmpmask,im_orig_mm[cvxyz]));
                twosum(&rss_orig_mm,&rss_orig_err_mm,&tmpsqr_mm);        
            }
            //flag_invalid_allZ[cy*numVectorX+cvx] = _MM_TESTZ1(tmpmask_sum);
            //if (!_MM_TESTZ1(tmpmask_sum)) {
                table_mask_xy[cy].push_back(cvx);
            //}
        }
    }

    SIMD_ELEMENT_TYPE *rss_orig     = (SIMD_ELEMENT_TYPE *)&rss_orig_mm;
    SIMD_ELEMENT_TYPE *rss_orig_err = (SIMD_ELEMENT_TYPE *)&rss_orig_err_mm;
    SIMD_ELEMENT_TYPE tmpsum = 0.0;
    SIMD_ELEMENT_TYPE tmperr = 0.0;
    for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
        twosum(&tmpsum,&tmperr,&rss_orig[counter]);
        twosum(&tmpsum,&tmperr,&rss_orig_err[counter]);
    }
    tmpsum += tmperr;
    rss_orig_inv = 1.0/((double)tmpsum);    

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    #ifdef MY_ASSERT
    /* declare variables */
    unsigned __int64 t_start, t_end;

    t_start = __rdtsc();
    #endif

    // /* Check for memory leak */
    //_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
    //_CrtSetBreakAlloc(117670);
    
	/* check numbers of input and output arguments */
	if (nrhs > 3 || nrhs < 1) {
		mexErrMsgTxt("Number of input argurments should be 1 to 3.");
	} else if (nlhs > 4) {
		mexErrMsgTxt("Too many output arguments!");
	}

	/* If number of input is 2 or more, initialize before main calculation.
	 If number of input is 1, go to main calculation without initialization. */
    if(nrhs > 1) { 
        initialize(nrhs,prhs); /* initialize host-side; memory allocation and store image */
    }

    //////#ifdef _OPENMP
    //////// set openmp stack size through KMP_STACKSIZE environment variable.
    ////////if (nlhs==3) {
    //////    int hoge = kmp_get_stacksize_s();
    //////    kmp_set_stacksize_s(numXWithPad*numy*numz*sizeof(SIMD_ELEMENT_TYPE)+16*1024*1024);
    //////    hoge = kmp_get_stacksize_s();
    ////////}
    //////#endif


    #ifdef MY_ASSERT
    t_end = __rdtsc();
    mexPrintf("initialize: %I64d[clocks]\n", t_end-t_start);
    t_start = __rdtsc();
    #endif


    /* store parameters */
    init_params(prhs[0]);

    #ifdef MY_ASSERT
    t_end = __rdtsc();
    mexPrintf("init_params: %I64d[clocks]\n", t_end-t_start);
    t_start = __rdtsc();
    #endif

    /* EM-like optimization */
    //em_optim();
    optimize();

    #ifdef MY_ASSERT
    t_end = __rdtsc();
    mexPrintf("em_optim: %I64d[clocks]\n", t_end-t_start);
    t_start = __rdtsc();
    #endif


    switch (nlhs) {
    case 4 :
        plhs[3] = mxCreateDoubleMatrix(finiter+1,2,mxREAL);
        setverbose(mxGetPr(plhs[3]));
    case 3 :
        plhs[2] = mxCreateNumericArray(numdims,dims,mxDOUBLE_CLASS,mxREAL);
        //setsynth(mxGetPr(plhs[2]));
        setsynth2(mxGetPr(plhs[2]));
    case 2 :
        plhs[1] = mxCreateDoubleScalar(score);
    case 1 :        
        plhs[0] = mxCreateDoubleMatrix(NUM_PARAMS_INPUT,numelipsoids,mxREAL);
        setresult(mxGetPr(plhs[0]));
    }

    #ifdef MY_ASSERT
    t_end = __rdtsc();
    mexPrintf("setresult: %I64d[clocks]\n", t_end-t_start);
    #endif
}



void calc_ofv(int idx_update){
    /* calculate RSS (residual sum of squares) */
    // int idx_update; // specify what should be updatated.  -1:calc_ofv, 0:imsynth, 1:pi_k, 2:mu, 3:sigma

    #ifdef MY_ASSERT    
        unsigned __int64 t_start, t_init=0, t_loop=0, t_sumup=0;
    #endif
   
    // setup OpenMP
    unsigned int numthreads = 1;
    unsigned int idx_thread = 0;
    #ifdef _OPENMP
        #pragma omp parallel firstprivate(table_elipsoid_z,idx_thread)
        { 
            numthreads = omp_get_num_threads();
            idx_thread = omp_get_thread_num();
    #endif
    unsigned int omp_offset = numelipsoids*idx_thread;

    
    // initialize result matrix
    rss_current_mm[idx_thread] = zeros_mm;
    rss_current_err_mm[idx_thread] = zeros_mm;
    switch (idx_update) {
    case -1: // for update ofv
        // nothing to do here
        break;
    case 0: // for update imsynth
        for (unsigned int cvxyz=0; cvxyz<numVectorX*numy*numz; cvxyz++) {
            im_synth_mm[cvxyz] = zeros_mm;
        }
        break;
    case 1: // for update pi_k
        #ifdef _OPENMP
            #pragma omp for
        #endif
        for (unsigned int celipsoid=0; celipsoid<numelipsoids*numthreads; celipsoid++) {            
            gky_mm[celipsoid] = zeros_mm;
            gky_err_mm[celipsoid] = zeros_mm;
        } // end of omp_for
        break;
    case 2: // for update mu
    case 3: // for update sigma
        #ifdef _OPENMP
            #pragma omp for
        #endif
        for (unsigned int celipsoid=0; celipsoid<numelipsoids*numthreads; celipsoid++) {
            for (unsigned int cseq=0; cseq<NUM_PARAMS_SEQ; cseq++) {
                seq_mm[celipsoid*NUM_PARAMS_SEQ+cseq] = zeros_mm;
                seq_err_mm[celipsoid*NUM_PARAMS_SEQ+cseq] = zeros_mm;
            }
        } // end of omp_for
        break;
    }


    // loop for xy vectors
    // once the parameters of ellipsoids in the xy position are calculated, 
    // then loop for z will be done with parallel handling of the ellipsoids.
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (unsigned int cy = 0; cy<numy; cy++) {
        __MM posY_mm = _MM_SET1(cy);

        for (unsigned int cx=0; cx<table_mask_xy[cy].size(); cx++) {
            unsigned int cvx = table_mask_xy[cy][cx];

            __MM posX_mm = _MM_ADD(_MM_SET1(cvx*SIMD_VECTOR_LENGTH),addx_mm);

            unsigned int cvalid = omp_offset;
            unsigned int tableidx = cy*numVectorX+cvx;
            __MM min_zstart_mm = numz_mm;
            __MM max_zend_mm = zeros_mm;

            for (unsigned int cz=0; cz<table_elipsoid_z.size(); cz++) {
                table_elipsoid_z[cz].clear();
            }

            #ifdef MY_ASSERT
                t_start = __rdtsc();
            #endif

            // loop for ellipsoid
            for (unsigned int c=0; c<table_elipsoid_xy[tableidx].size(); c++) {

                // load parameters of the ellipsoid
                unsigned int offset = table_elipsoid_xy[tableidx][c]*NUM_PARAMS_PROC;
                __MM xc_mm      = params_proc_mm[offset +  1];
                __MM yc_mm      = params_proc_mm[offset +  2];
                __MM zc_mm      = params_proc_mm[offset +  3];
                __MM Sinv_11_mm = params_proc_mm[offset +  4];
                __MM Sinv_12_mm = params_proc_mm[offset +  5];
                __MM Sinv_13_mm = params_proc_mm[offset +  6];
                __MM Sinv_22_mm = params_proc_mm[offset +  7];
                __MM Sinv_23_mm = params_proc_mm[offset +  8];
                __MM Sinv_33_mm = params_proc_mm[offset +  9];
                __MM cx1_mm     = params_proc_mm[offset + 10];
                __MM cy1_mm     = params_proc_mm[offset + 11];
                __MM cx2_mm     = params_proc_mm[offset + 12];
                __MM cy2_mm     = params_proc_mm[offset + 13];
                __MM cxy_mm     = params_proc_mm[offset + 14];
                __MM cp_mm      = params_proc_mm[offset + 15];
                __MM bz_mm      = params_proc_mm[offset + 16];


                // calc interemediate parameters of the ellipsoid at the current xy position
                __MM xd_mm  = _MM_SUB(posX_mm, xc_mm);
                __MM yd_mm  = _MM_SUB(posY_mm, yc_mm);
                __MM xd2_mm = _MM_SQR(xd_mm);
                __MM yd2_mm = _MM_SQR(yd_mm);
                __MM xyd_mm = _MM_MUL(xd_mm,yd_mm);
        
                __MM det_mm = _MM_ADD(_MM_ADD(_MM_ADD(cp_mm,_MM_MUL(cx2_mm,xd2_mm)),
                                                            _MM_MUL(cy2_mm,yd2_mm)),
                                                            _MM_MUL(cxy_mm,xyd_mm));
                __MM flagPosDet = _MM_CMP(det_mm,zeros_mm,_CMP_GE_OS);
                if (_MM_TESTZ1(flagPosDet)) { continue; }

                __MM sqrtdet_mm = _MM_SQRT(_MM_AND(det_mm,flagPosDet));
                __MM cx1xcy1y_mm = _MM_ADD(_MM_ADD(_MM_MUL(cx1_mm,xd_mm),_MM_MUL(cy1_mm,yd_mm)),zc_mm);
                __MM tmpIdxMinZ_mm =  _MM_CEIL(_MM_SUB(cx1xcy1y_mm,sqrtdet_mm)); 
                __MM tmpIdxMaxZ_mm = _MM_FLOOR(_MM_ADD(cx1xcy1y_mm,sqrtdet_mm)); 
            
                __MM zstart_mm = _MM_MAX(tmpIdxMinZ_mm, zeros_mm); 
                __MM   zend_mm = _MM_MIN(tmpIdxMaxZ_mm,numzm1_mm); 

                __MM flagInRange = _MM_AND(flagPosDet,_MM_AND(
                        _MM_CMP(zstart_mm,numzm1_mm,_CMP_LE_OS), 
                        _MM_CMP(  zend_mm, zeros_mm,_CMP_GE_OS)));
                zstart_mm = _MM_BLENDV( numz_mm,zstart_mm,flagInRange);
                zend_mm = _MM_BLENDV(zeros_mm,_MM_ADD(zend_mm,ones_mm),flagInRange);

                __MM zd_mm = _MM_SUB(zstart_mm,zc_mm);

                __MM xs_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_11_mm),
                                             _MM_MUL(yd_mm,Sinv_12_mm)),
                                             _MM_MUL(zd_mm,Sinv_13_mm));
                __MM ys_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_12_mm),
                                             _MM_MUL(yd_mm,Sinv_22_mm)),
                                             _MM_MUL(zd_mm,Sinv_23_mm));
                __MM zs_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_13_mm),
                                             _MM_MUL(yd_mm,Sinv_23_mm)),
                                             _MM_MUL(zd_mm,Sinv_33_mm));
            
                __MM intensity0_mm = _MM_EXP(_MM_MUL(mhalf_mm,_MM_ADD(_MM_ADD(_MM_MUL(xd_mm,xs_mm),
                                                                              _MM_MUL(yd_mm,ys_mm)),
                                                                              _MM_MUL(zd_mm,zs_mm))));
                __MM az_mm = zs_mm;
                __MM fz0_mm = _MM_EXP(_MM_SUB(bz_mm,az_mm));


                int tmpzmin = hmin(zstart_mm);
                int tmpzmax = hmax(  zend_mm);
                for (int cz=tmpzmin; cz<tmpzmax; cz++) {
                    table_elipsoid_z[cz].push_back(cvalid);
                }               
                
                min_zstart_mm = _MM_MIN(min_zstart_mm,zstart_mm);
                max_zend_mm = _MM_MAX( max_zend_mm, zend_mm);
                idx_validelip[cvalid] = table_elipsoid_xy[tableidx][c];
                mem_int0_mm[cvalid] = intensity0_mm;
                mem_fz0_mm[cvalid] = fz0_mm;
                mem_int_mm[cvalid] = zeros_mm;
                mem_fz_mm[cvalid] = zeros_mm;
                idx_Zstart_mm[cvalid] = zstart_mm;
                idx_Zend_mm[cvalid] = zend_mm;


                switch (idx_update) { 
                case -1: // for update ofv
                case 0: // for update imsynth
                    // nothing to do here
                    break;
                case 1: // for update pi_k
                    // nothing to do here
                    break;
                case 2: // for update mu
                case 3: // for update sigma
                    mem_sxm_1_mm[cvalid] = _MM_ADD(_MM_MUL(xd_mm,Sinv_11_mm),
                        _MM_MUL(yd_mm,Sinv_12_mm));
                    mem_sxm_2_mm[cvalid] = _MM_ADD(_MM_MUL(xd_mm,Sinv_12_mm),
                        _MM_MUL(yd_mm,Sinv_22_mm));
                    mem_sxm_3_mm[cvalid] = _MM_ADD(_MM_MUL(xd_mm,Sinv_13_mm),
                        _MM_MUL(yd_mm,Sinv_23_mm));
                    break;
                }
                
                cvalid++;           

            }
            unsigned int numvalid = cvalid;


            int minZStart = hmin(min_zstart_mm);
            int   maxZEnd = hmax(  max_zend_mm);


            #ifdef MY_ASSERT
                t_init += __rdtsc() - t_start;
                t_start = __rdtsc();
            #endif


             /* loop for z */
            unsigned int cvxyz;
            for (int cz = minZStart; cz<maxZEnd; cz++ ){
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;

                if (table_elipsoid_z[cz].size()==0) { continue; }

                // update intensity of each ellipsoid at the voxel
                __MM posZ_mm = _MM_SET1(cz);
                for (int ct=0; ct<table_elipsoid_z[cz].size(); ct++) {
                    cvalid = table_elipsoid_z[cz][ct];
                    __MM flagstart_mm = _MM_CMP(idx_Zstart_mm[cvalid],posZ_mm,_CMP_EQ_OS);
                    __MM flagend_mm   = _MM_CMP(  idx_Zend_mm[cvalid],posZ_mm,_CMP_NEQ_OS);
                    __MM exp2bz_mm = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC + 17];
                
                    mem_int_mm[cvalid] = _MM_AND(_MM_ADD(_MM_MUL(mem_int_mm[cvalid],mem_fz_mm[cvalid]),
                                                         _MM_AND(mem_int0_mm[cvalid],flagstart_mm)),flagend_mm);
                    mem_fz_mm[cvalid]  =         _MM_ADD(_MM_MUL(mem_fz_mm[cvalid],exp2bz_mm),
                                                         _MM_AND(mem_fz0_mm[cvalid],flagstart_mm));
                }

                // when the voxel was masked, skip the following calculations
                if (_MM_TESTZ1(mask_mm[cvxyz])) { continue; }
                
                // calc integrals 
                __MM intensity_mm = zeros_mm;
                for (int ct=0; ct<table_elipsoid_z[cz].size(); ct++) {
                    cvalid = table_elipsoid_z[cz][ct];
                    __MM pi_k_mm   = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC +  0];
                    __MM tmpint_mm = _MM_MUL(mem_int_mm[cvalid],pi_k_mm);
                    intensity_mm = _MM_ADD(intensity_mm,tmpint_mm);

                    // for update pi_k
                    if (idx_update==1) {
                        unsigned int tmpidx = idx_validelip[cvalid] + omp_offset;
                        __MM tmpgky_mm = _MM_MUL(mem_int_mm[cvalid],im_orig_mm[cvxyz]);
                        twosum(&gky_mm[tmpidx],&gky_err_mm[tmpidx],&tmpgky_mm);
                    }

                }
                
                // update residuals
                __MM residual_mm = _MM_SUB(im_orig_mm[cvxyz], intensity_mm);
                __MM rss_old_mm  = _MM_SQR(im_orig_mm[cvxyz]);
                __MM rss_tmp_mm  = _MM_SQR(residual_mm);
                twosub(&rss_current_mm[idx_thread],&rss_current_err_mm[idx_thread],&rss_old_mm);
                twosum(&rss_current_mm[idx_thread],&rss_current_err_mm[idx_thread],&rss_tmp_mm);

                // for update ofv
                if (idx_update==-1) {
                    // nothing to do here
                    continue;
                }

                // for update imsynth
                if (idx_update==0) { 
                    im_synth_mm[cvxyz] = _MM_AND(intensity_mm,mask_mm[cvxyz]);
                    continue; 
                }

                // for update pi_k
                if (idx_update==1) {
                    // nothing to do here
                    continue; 
                } 

                // for update mu and sigma
                for (int ct=0; ct<table_elipsoid_z[cz].size(); ct++) {
                    int cvalid = table_elipsoid_z[cz][ct];

                    __MM zc_mm = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC + 3];
                    __MM Sinv_13_mm = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC + 6];
                    __MM Sinv_23_mm = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC + 8];
                    __MM Sinv_33_mm = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC + 9];
                    __MM zd_mm = _MM_SUB(posZ_mm,zc_mm);
                    __MM sxm_1_mm = _MM_ADD(mem_sxm_1_mm[cvalid],_MM_MUL(zd_mm,Sinv_13_mm));
                    __MM sxm_2_mm = _MM_ADD(mem_sxm_2_mm[cvalid],_MM_MUL(zd_mm,Sinv_23_mm));
                    __MM sxm_3_mm = _MM_ADD(mem_sxm_3_mm[cvalid],_MM_MUL(zd_mm,Sinv_33_mm));
                    __MM sxm_mm[3] = {sxm_1_mm,sxm_2_mm,sxm_3_mm};

                    unsigned int tmpidx = (idx_validelip[cvalid] + omp_offset)*NUM_PARAMS_SEQ;
                    unsigned int d1=1, d2=4, d3=10, d4=20;
                    __MM v0_mm, v1_mm, v2_mm, v3_mm, v4_mm;

                    // for updating mu, v0, v1, and v2 are required.
                    // for updating sigma, v2 and v4 are required.
                    

                    v0_mm = _MM_MUL(mem_int_mm[cvalid],residual_mm); // g_k*r
                    twosum(&seq_mm[tmpidx+0],&seq_err_mm[tmpidx+0],&v0_mm);

                    if (idx_update==2) { // for update mu
                        v1_mm = _MM_MUL(v0_mm,sxm_mm[0]);
                        twosum(&seq_mm[ tmpidx+d1+s2i[0][0][0][0]],
                            &seq_err_mm[tmpidx+d1+s2i[0][0][0][0]],&v1_mm);
                        v1_mm = _MM_MUL(v0_mm,sxm_mm[1]);
                        twosum(&seq_mm[ tmpidx+d1+s2i[1][0][0][0]],
                            &seq_err_mm[tmpidx+d1+s2i[1][0][0][0]],&v1_mm);
                        v1_mm = _MM_MUL(v0_mm,sxm_mm[2]);
                        twosum(&seq_mm[ tmpidx+d1+s2i[2][0][0][0]],
                            &seq_err_mm[tmpidx+d1+s2i[2][0][0][0]],&v1_mm);
                    }

                    // for updating mu and sigma
                    v2_mm = _MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[0]);
                    twosum(&seq_mm[ tmpidx+d2+s2i[0][0][0][0]],
                        &seq_err_mm[tmpidx+d2+s2i[0][0][0][0]],&v2_mm);
                    v2_mm = _MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[1]);
                    twosum(&seq_mm[ tmpidx+d2+s2i[0][1][0][0]],
                        &seq_err_mm[tmpidx+d2+s2i[0][1][0][0]],&v2_mm);
                    v2_mm = _MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d2+s2i[0][2][0][0]],
                        &seq_err_mm[tmpidx+d2+s2i[0][2][0][0]],&v2_mm);
                    v2_mm = _MM_MUL(_MM_MUL(v0_mm,sxm_mm[1]),sxm_mm[1]);
                    twosum(&seq_mm[ tmpidx+d2+s2i[1][1][0][0]],
                        &seq_err_mm[tmpidx+d2+s2i[1][1][0][0]],&v2_mm);
                    v2_mm = _MM_MUL(_MM_MUL(v0_mm,sxm_mm[1]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d2+s2i[1][2][0][0]],
                        &seq_err_mm[tmpidx+d2+s2i[1][2][0][0]],&v2_mm);
                    v2_mm = _MM_MUL(_MM_MUL(v0_mm,sxm_mm[2]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d2+s2i[2][2][0][0]],
                        &seq_err_mm[tmpidx+d2+s2i[2][2][0][0]],&v2_mm);

                    if (idx_update==2) { continue; }

                    // for update sigma
                    /* v3 was not used for alternative optimization method
                    v3_mm = _MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[0]),sxm_mm[0]);
                    twosum(&seq_mm[ tmpidx+d3+s2i[0][0][0][0]],
                    &seq_err_mm[tmpidx+d3+s2i[0][0][0][0]],&v3_mm);
                    v3_mm = _MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[0]),sxm_mm[1]);
                    twosum(&seq_mm[ tmpidx+d3+s2i[0][0][1][0]],
                    &seq_err_mm[tmpidx+d3+s2i[0][0][1][0]],&v3_mm);
                    v3_mm = _MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[0]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d3+s2i[0][0][2][0]],
                    &seq_err_mm[tmpidx+d3+s2i[0][0][2][0]],&v3_mm);
                    v3_mm = _MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[1]),sxm_mm[1]);
                    twosum(&seq_mm[ tmpidx+d3+s2i[0][1][1][0]],
                    &seq_err_mm[tmpidx+d3+s2i[0][1][1][0]],&v3_mm);
                    v3_mm = _MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[1]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d3+s2i[0][1][2][0]],
                    &seq_err_mm[tmpidx+d3+s2i[0][1][2][0]],&v3_mm);
                    v3_mm = _MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[2]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d3+s2i[0][2][2][0]],
                    &seq_err_mm[tmpidx+d3+s2i[0][2][2][0]],&v3_mm);
                    v3_mm = _MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[1]),sxm_mm[1]),sxm_mm[1]);
                    twosum(&seq_mm[ tmpidx+d3+s2i[1][1][1][0]],
                    &seq_err_mm[tmpidx+d3+s2i[1][1][1][0]],&v3_mm);
                    v3_mm = _MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[1]),sxm_mm[1]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d3+s2i[1][1][2][0]],
                    &seq_err_mm[tmpidx+d3+s2i[1][1][2][0]],&v3_mm);
                    v3_mm = _MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[1]),sxm_mm[2]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d3+s2i[1][2][2][0]],
                    &seq_err_mm[tmpidx+d3+s2i[1][2][2][0]],&v3_mm);
                    v3_mm = _MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[2]),sxm_mm[2]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d3+s2i[2][2][2][0]],
                    &seq_err_mm[tmpidx+d3+s2i[2][2][2][0]],&v3_mm);
                    */

                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[0]),sxm_mm[0]),sxm_mm[0]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[0][0][0][0]],
                        &seq_err_mm[tmpidx+d4+s2i[0][0][0][0]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[0]),sxm_mm[0]),sxm_mm[1]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[0][0][0][1]],
                        &seq_err_mm[tmpidx+d4+s2i[0][0][0][1]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[0]),sxm_mm[0]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[0][0][0][2]],
                        &seq_err_mm[tmpidx+d4+s2i[0][0][0][2]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[0]),sxm_mm[1]),sxm_mm[1]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[0][0][1][1]],
                        &seq_err_mm[tmpidx+d4+s2i[0][0][1][1]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[0]),sxm_mm[1]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[0][0][1][2]],
                        &seq_err_mm[tmpidx+d4+s2i[0][0][1][2]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[0]),sxm_mm[2]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[0][0][2][2]],
                        &seq_err_mm[tmpidx+d4+s2i[0][0][2][2]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[1]),sxm_mm[1]),sxm_mm[1]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[0][1][1][1]],
                        &seq_err_mm[tmpidx+d4+s2i[0][1][1][1]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[1]),sxm_mm[1]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[0][1][1][2]],
                        &seq_err_mm[tmpidx+d4+s2i[0][1][1][2]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[1]),sxm_mm[2]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[0][1][2][2]],
                        &seq_err_mm[tmpidx+d4+s2i[0][1][2][2]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[0]),sxm_mm[2]),sxm_mm[2]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[0][2][2][2]],
                        &seq_err_mm[tmpidx+d4+s2i[0][2][2][2]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[1]),sxm_mm[1]),sxm_mm[1]),sxm_mm[1]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[1][1][1][1]],
                        &seq_err_mm[tmpidx+d4+s2i[1][1][1][1]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[1]),sxm_mm[1]),sxm_mm[1]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[1][1][1][2]],
                        &seq_err_mm[tmpidx+d4+s2i[1][1][1][2]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[1]),sxm_mm[1]),sxm_mm[2]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[1][1][2][2]],
                        &seq_err_mm[tmpidx+d4+s2i[1][1][2][2]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[1]),sxm_mm[2]),sxm_mm[2]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[1][2][2][2]],
                        &seq_err_mm[tmpidx+d4+s2i[1][2][2][2]],&v4_mm);
                    v4_mm = _MM_MUL(_MM_MUL(_MM_MUL(_MM_MUL(v0_mm,sxm_mm[2]),sxm_mm[2]),sxm_mm[2]),sxm_mm[2]);
                    twosum(&seq_mm[ tmpidx+d4+s2i[2][2][2][2]],
                        &seq_err_mm[tmpidx+d4+s2i[2][2][2][2]],&v4_mm);

                } // end of for update mu and sigma

            } // end of loop for z

            #ifdef MY_ASSERT
                t_loop += __rdtsc() - t_start;
            #endif

        }
    } // end of omp_for


    
    #ifdef _OPENMP

    // for updating  pi_k
    // aggregate the local copies into global array
    if (idx_update==1) {
        #pragma omp for
        for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
            for (unsigned int p=1; p<numthreads; p++) {
                twosum(&gky_mm[celipsoid],&gky_err_mm[celipsoid],&gky_mm[numelipsoids*p+celipsoid]);
                twosum(&gky_mm[celipsoid],&gky_err_mm[celipsoid],&gky_err_mm[numelipsoids*p+celipsoid]);
            }
        } // end of omp_for
    }

    // for updating  mu and sigma
    // aggregate the local copies into global array
    if (idx_update>=2) {
        #pragma omp for
        for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
            for (unsigned int p=1; p<numthreads; p++) {            

                unsigned int offset1 = NUM_PARAMS_SEQ*celipsoid;
                unsigned int offset2 = NUM_PARAMS_SEQ*(numelipsoids*p+celipsoid);
                for (unsigned int cseq=0; cseq<NUM_PARAMS_SEQ; cseq++) {
                    twosum(&seq_mm[offset1+cseq],&seq_err_mm[offset1+cseq],&seq_mm[offset2+cseq]);
                    twosum(&seq_mm[offset1+cseq],&seq_err_mm[offset1+cseq],&seq_err_mm[offset2+cseq]);
                }

            }
        } // end of omp_for
    }


    } // end of parallel

    // for updateing pi, mu, sigma
    // aggregate the local copies of rss_current_mm into global array
    if (idx_update>=1) {
        for (unsigned int p=1; p<numthreads; p++) {
            twosum(&rss_current_mm[0],&rss_current_err_mm[0],&rss_current_mm[p]);
            twosum(&rss_current_mm[0],&rss_current_err_mm[0],&rss_current_err_mm[p]);
        }
    }
    #endif
    

    #ifdef MY_ASSERT
        t_start = __rdtsc();
    #endif


    /* sum up g_k*y for SIMD */
    if (idx_update==1) {
        for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
            unsigned int tmpidx = SIMD_VECTOR_LENGTH*celipsoid;
            for (unsigned int counter=1; counter<SIMD_VECTOR_LENGTH; counter++) {
                twosum(&gky[tmpidx],&gky_err[tmpidx],&gky[tmpidx+counter]);
                twosum(&gky[tmpidx],&gky_err[tmpidx],&gky_err[tmpidx+counter]);
            }
            twosum(&gky[tmpidx],&gky_err[tmpidx],&gky_err[tmpidx]);
        }
    }


    /* sum up V_k for SIMD */
    if (idx_update>2) {
        for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
            for (unsigned int cseq=0; cseq<NUM_PARAMS_SEQ; cseq++) {
                unsigned int tmpidx = SIMD_VECTOR_LENGTH*(NUM_PARAMS_SEQ*celipsoid+cseq);
                for (unsigned int counter=1; counter<SIMD_VECTOR_LENGTH; counter++) {                
                    twosum(&seq[tmpidx],&seq_err[tmpidx],&seq[tmpidx+counter]);
                    twosum(&seq[tmpidx],&seq_err[tmpidx],&seq_err[tmpidx+counter]);
                }
                twosum(&seq[tmpidx],&seq_err[tmpidx],&seq_err[tmpidx]);
            }     
        }
    }

    // sum up score for SIMD then update score
    if (idx_update>=1 || idx_update==-1) {
        SIMD_ELEMENT_TYPE tmpsum = 0.0;
        SIMD_ELEMENT_TYPE tmpsum_err = 0.0;
        SIMD_ELEMENT_TYPE *rss_current = (SIMD_ELEMENT_TYPE *)&rss_current_mm[0];
        SIMD_ELEMENT_TYPE *rss_current_err = (SIMD_ELEMENT_TYPE *)&rss_current_err_mm[0];
        for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
            twosum(&tmpsum,&tmpsum_err,&rss_current[counter]);
            twosum(&tmpsum,&tmpsum_err,&rss_current_err[counter]);
        }
        tmpsum += tmpsum_err;
        // score = 1 + ((double)tmpsum) * rss_orig_inv;        
        calc_penalty();
        score = 1 + (((double)tmpsum) + penalty) * rss_orig_inv;

        
    }

    #ifdef MY_ASSERT
        t_sumup = __rdtsc() - t_start;
    #endif

    #ifdef MY_ASSERT    
        mexPrintf("calc_ofv.t_init: %I64d[clocks]\n", t_init);
        mexPrintf("calc_ofv.t_loop: %I64d[clocks]\n", t_loop);
        mexPrintf("calc_ofv.t_sumup: %I64d[clocks]\n", t_sumup);
    #endif

}


void calc_W(int idx_update) {
/* calc W (integral of gk*gl) for update pi_k, mu, and sigma */
    // int idx_update; // specify what should be updatated. -1:calc_ofv, 0:imsynth, 1:pi_k, 2:mu, 3:sigma

    // for update im_synth and calc_ofv, nothing to do in this function
    if (idx_update<=0) { return; }

    // initialize dsyevd for diagonalizing A
    lapack_int lwork = -1;
    lapack_int liwork = -1;
    if (idx_update>=2) { // for update mu and sigma    
        lapack_int n = 3;
        lapack_int tmpinfo;
        double *work = new double[1];
        lapack_int *iwork = new lapack_int[1];
        double *A = new double[n*n];
        double *w = new double[n];
        A[0]=10; A[1]= 0; A[2]= 0;
        A[3]= 0; A[4]=10; A[5]= 0;
        A[6]= 0; A[7]= 0; A[8]=10;
        dsyevd("V","U",&n,A,&n,w,work,&lwork,iwork,&liwork,&tmpinfo);
        lwork = __max(1,((lapack_int)(work[0])));
        liwork = __max(1,iwork[0]);
        delete[] A;
        delete[] w;
        delete[] work;
        delete[] iwork;
    }


    // initialize return matrix
    switch (idx_update) {
    case 1: // update for pi_k
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for (unsigned int k=0; k<numelipsoids*NUM_PARAMS_PI_K; k++) {   
            grad_pi_k[k] = 0.0;
            for (unsigned int l=0; l<numelipsoids*NUM_PARAMS_PI_K; l++) {
                hess_pi_k[numelipsoids*NUM_PARAMS_PI_K*k+l] = 0.0;                
            }
        }
        break;
    case 2: // update for mu
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for (unsigned int k=0; k<numelipsoids*NUM_PARAMS_MU; k++) {   
            grad_mu[k] = 0.0;
            for (unsigned int l=0; l<numelipsoids*NUM_PARAMS_MU; l++) {
                hess_mu[numelipsoids*NUM_PARAMS_MU*k+l] = 0.0;
            }
        }
        break;
    case 3: // update for sigma
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for (unsigned int k=0; k<numelipsoids*NUM_PARAMS_SIGMA; k++) {   
            grad_sigma[k] = 0.0;
            for (unsigned int l=0; l<numelipsoids*NUM_PARAMS_SIGMA; l++) {
                hess_sigma[numelipsoids*NUM_PARAMS_SIGMA*k+l] = 0.0;               
            }
        }
        break;
    }

    // loop for ellipsoid
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (unsigned int k=0; k<numelipsoids; k++) {        
        SIMD_ELEMENT_TYPE A_k_11   = params_calcG[k*NUM_PARAMS_CALCG + 0];
        SIMD_ELEMENT_TYPE A_k_12   = params_calcG[k*NUM_PARAMS_CALCG + 1];
        SIMD_ELEMENT_TYPE A_k_13   = params_calcG[k*NUM_PARAMS_CALCG + 2];
        SIMD_ELEMENT_TYPE A_k_22   = params_calcG[k*NUM_PARAMS_CALCG + 3];
        SIMD_ELEMENT_TYPE A_k_23   = params_calcG[k*NUM_PARAMS_CALCG + 4];
        SIMD_ELEMENT_TYPE A_k_33   = params_calcG[k*NUM_PARAMS_CALCG + 5];
        SIMD_ELEMENT_TYPE B_k_1    = params_calcG[k*NUM_PARAMS_CALCG + 6];
        SIMD_ELEMENT_TYPE B_k_2    = params_calcG[k*NUM_PARAMS_CALCG + 7];
        SIMD_ELEMENT_TYPE B_k_3    = params_calcG[k*NUM_PARAMS_CALCG + 8];
        SIMD_ELEMENT_TYPE C_k      = params_calcG[k*NUM_PARAMS_CALCG + 9];


        // for update mu and sigma
        lapack_int n=3, info, *iwork;
        double *work, *A, *w;        
        SIMD_ELEMENT_TYPE pi_k, xc_k, yc_k, zc_k;
        double V_phi, V_i[3], V_ij[3][3], /*V_ijp[3][3][3],*/ V_ijpq[3][3][3][3];

        if (idx_update>=2) {
            // setup dsyevd for diagonalizing A
            work = new double[lwork];
            iwork = new lapack_int[liwork];
            A = new double[n*n];
            w = new double[n];

            // load params
            pi_k   = params[NUM_PARAMS_INPUT*k + 0];
            xc_k   = params[NUM_PARAMS_INPUT*k + 1];
            yc_k   = params[NUM_PARAMS_INPUT*k + 2];
            zc_k   = params[NUM_PARAMS_INPUT*k + 3];

            // set V_k
            unsigned int offset = SIMD_VECTOR_LENGTH*NUM_PARAMS_SEQ*k;
            V_phi = seq[offset + 0*SIMD_VECTOR_LENGTH];            
            for (unsigned int i=0; i<3; i++) {
                V_i[i] = seq[offset + (s2i[i][0][0][0]+1)*SIMD_VECTOR_LENGTH];
                for (unsigned int j=0; j<3; j++) {                
                    V_ij[i][j] = seq[offset + (s2i[i][j][0][0]+4)*SIMD_VECTOR_LENGTH];
                    for (unsigned int p=0; p<3; p++) {
                        //V_ijp[i][j][p] = seq[offset + (s2i[i][j][p][0]+10)*SIMD_VECTOR_LENGTH];
                        if (idx_update==3) // for update sigma
                        for (unsigned int q=0; q<3; q++) {
                            V_ijpq[i][j][p][q] = seq[offset + (s2i[i][j][p][q]+20)*SIMD_VECTOR_LENGTH];
                        }
                    }
                }
            } // end of loop for set V_k
        } // end of if


        // update gradient
        {
            SIMD_ELEMENT_TYPE lambda, xd1, yd1, zd1, xd2, yd2, zd2;
            SIMD_ELEMENT_TYPE pen_x = 0.0;
            SIMD_ELEMENT_TYPE pen_y = 0.0;
            SIMD_ELEMENT_TYPE pen_z = 0.0;
            switch (idx_update) {
            case 1: // for update pi_k
                grad_pi_k[k] = gky[SIMD_VECTOR_LENGTH*k] + lambda_pi*params_init[NUM_PARAMS_INPUT*k];
                break;
            case 2: // for update mu; because H\Delta\theta = -G, the sign of gradient is inverted. 
                xd1 = params_diff[NUM_PARAMS_INPUT*k + 1];
                yd1 = params_diff[NUM_PARAMS_INPUT*k + 2];
                zd1 = params_diff[NUM_PARAMS_INPUT*k + 3];
                for (unsigned int celipsoid2=0; celipsoid2<numelipsoids; celipsoid2++) {
                    lambda = lambda_mu[k*numelipsoids+celipsoid2];
                    xd2 = params_diff[NUM_PARAMS_INPUT*celipsoid2 + 1];
                    yd2 = params_diff[NUM_PARAMS_INPUT*celipsoid2 + 2];
                    zd2 = params_diff[NUM_PARAMS_INPUT*celipsoid2 + 3];
                    pen_x += lambda * (xd1-xd2);
                    pen_y += lambda * (yd1-yd2);
                    pen_z += lambda * (zd1-zd2);
                }
                grad_mu[NUM_PARAMS_MU*k+0] = 2*(pi_k*V_i[0]-2*pen_x)*rss_orig_inv;
                grad_mu[NUM_PARAMS_MU*k+1] = 2*(pi_k*V_i[1]-2*pen_y)*rss_orig_inv;
                grad_mu[NUM_PARAMS_MU*k+2] = 2*(pi_k*V_i[2]-2*pen_z)*rss_orig_inv;
                break;
            case 3: // for update sigma; because H\Delta\theta = -G, the sign of gradient is inverted.
                unsigned int idx1 = NUM_PARAMS_SIGMA*k;
                unsigned int idx2 = NUM_PARAMS_INPUT*k;
                grad_sigma[idx1+0] =   (pi_k*V_ij[0][0] - 2*lambda_sigma*params_diff[idx2+4])*rss_orig_inv;
                grad_sigma[idx1+1] = 2*(pi_k*V_ij[0][1] - 2*lambda_sigma*params_diff[idx2+5])*rss_orig_inv;
                grad_sigma[idx1+2] = 2*(pi_k*V_ij[0][2] - 2*lambda_sigma*params_diff[idx2+6])*rss_orig_inv;
                grad_sigma[idx1+3] =   (pi_k*V_ij[1][1] - 2*lambda_sigma*params_diff[idx2+7])*rss_orig_inv;
                grad_sigma[idx1+4] = 2*(pi_k*V_ij[1][2] - 2*lambda_sigma*params_diff[idx2+8])*rss_orig_inv;
                grad_sigma[idx1+5] =   (pi_k*V_ij[2][2] - 2*lambda_sigma*params_diff[idx2+9])*rss_orig_inv;
                break;
            }
        }

        // loop for 2nd ellipsoid (l)
        for (unsigned int l=k; l<numelipsoids; l++) {
            SIMD_ELEMENT_TYPE A_l_11   = params_calcG[l*NUM_PARAMS_CALCG + 0];
            SIMD_ELEMENT_TYPE A_l_12   = params_calcG[l*NUM_PARAMS_CALCG + 1];
            SIMD_ELEMENT_TYPE A_l_13   = params_calcG[l*NUM_PARAMS_CALCG + 2];
            SIMD_ELEMENT_TYPE A_l_22   = params_calcG[l*NUM_PARAMS_CALCG + 3];
            SIMD_ELEMENT_TYPE A_l_23   = params_calcG[l*NUM_PARAMS_CALCG + 4];
            SIMD_ELEMENT_TYPE A_l_33   = params_calcG[l*NUM_PARAMS_CALCG + 5];
            SIMD_ELEMENT_TYPE B_l_1    = params_calcG[l*NUM_PARAMS_CALCG + 6];
            SIMD_ELEMENT_TYPE B_l_2    = params_calcG[l*NUM_PARAMS_CALCG + 7];
            SIMD_ELEMENT_TYPE B_l_3    = params_calcG[l*NUM_PARAMS_CALCG + 8];
            SIMD_ELEMENT_TYPE C_l      = params_calcG[l*NUM_PARAMS_CALCG + 9];

            // for update mu and sigma
            SIMD_ELEMENT_TYPE pi_l, xc_l, yc_l, zc_l;
            if (idx_update>=2) {
                pi_l   = params[NUM_PARAMS_INPUT*l + 0];
                xc_l   = params[NUM_PARAMS_INPUT*l + 1];
                yc_l   = params[NUM_PARAMS_INPUT*l + 2];
                zc_l   = params[NUM_PARAMS_INPUT*l + 3];
            }
            
            // for update mu, regularization term should be added before pruning
            if (idx_update==2) {
                unsigned int offset1 = NUM_PARAMS_MU*NUM_PARAMS_MU*numelipsoids*k+NUM_PARAMS_MU*l;
                unsigned int offset2 = NUM_PARAMS_MU*NUM_PARAMS_MU*numelipsoids*l+NUM_PARAMS_MU*k;
                if (k==l) {
                    for (unsigned int i=0; i<3; i++) { // mu_{k,i}
                        //hess_mu[offset1+NUM_PARAMS_MU*numelipsoids*i+i] +=  4*lambda_mu_sum[k];                        
                        hess_mu[offset1+NUM_PARAMS_MU*numelipsoids*i+i] +=  4*lambda_mu_sum[k]*rss_orig_inv;                        
                    }
                } else {
                    for (unsigned int i=0; i<3; i++) { // mu_{k,i}
                        //hess_mu[offset1+NUM_PARAMS_MU*numelipsoids*i+i] += -4*lambda_mu[k*numelipsoids+l];
                        //hess_mu[offset2+NUM_PARAMS_MU*numelipsoids*i+i] += -4*lambda_mu[k*numelipsoids+l];
                        hess_mu[offset1+NUM_PARAMS_MU*numelipsoids*i+i] += -4*lambda_mu[k*numelipsoids+l]*rss_orig_inv;
                        hess_mu[offset2+NUM_PARAMS_MU*numelipsoids*i+i] += -4*lambda_mu[k*numelipsoids+l]*rss_orig_inv;
                    }
                }
            }

            if ( idx_VXstart[k] > idx_VXend[l]
              || idx_VXstart[l] > idx_VXend[k]
              || idx_Ystart[k] > idx_Yend[l]
              || idx_Ystart[l] > idx_Yend[k] ) { continue; }


            SIMD_ELEMENT_TYPE A_11 = A_k_11 + A_l_11;
            SIMD_ELEMENT_TYPE A_22 = A_k_22 + A_l_22;
            SIMD_ELEMENT_TYPE A_33 = A_k_33 + A_l_33;
            SIMD_ELEMENT_TYPE A_12 = A_k_12 + A_l_12;
            SIMD_ELEMENT_TYPE A_13 = A_k_13 + A_l_13;
            SIMD_ELEMENT_TYPE A_23 = A_k_23 + A_l_23;

            SIMD_ELEMENT_TYPE B_1 = B_k_1 + B_l_1;
            SIMD_ELEMENT_TYPE B_2 = B_k_2 + B_l_2;
            SIMD_ELEMENT_TYPE B_3 = B_k_3 + B_l_3;

            SIMD_ELEMENT_TYPE C = C_k + C_l;


            SIMD_ELEMENT_TYPE detA =    A_11*A_22*A_33 
                + 2*A_12*A_13*A_23
                -   A_11*A_23*A_23 
                -   A_12*A_12*A_33
                -   A_13*A_22*A_13;

            SIMD_ELEMENT_TYPE detA_inv = 1/detA;

            SIMD_ELEMENT_TYPE Ainv_11 = (A_22*A_33 - A_23*A_23)*detA_inv;
            SIMD_ELEMENT_TYPE Ainv_22 = (A_11*A_33 - A_13*A_13)*detA_inv;
            SIMD_ELEMENT_TYPE Ainv_33 = (A_11*A_22 - A_12*A_12)*detA_inv;
            SIMD_ELEMENT_TYPE Ainv_12 = (A_13*A_23 - A_12*A_33)*detA_inv;
            SIMD_ELEMENT_TYPE Ainv_13 = (A_12*A_23 - A_13*A_22)*detA_inv;
            SIMD_ELEMENT_TYPE Ainv_23 = (A_12*A_13 - A_11*A_23)*detA_inv;

            SIMD_ELEMENT_TYPE AB[3] = { Ainv_11*B_1 + Ainv_12*B_2 + Ainv_13*B_3,
                                        Ainv_12*B_1 + Ainv_22*B_2 + Ainv_23*B_3,
                                        Ainv_13*B_1 + Ainv_23*B_2 + Ainv_33*B_3};

            SIMD_ELEMENT_TYPE BAB =  B_1*AB[0] + B_2*AB[1] + B_3*AB[2];

            SIMD_ELEMENT_TYPE gg = sqrt(8*M_PI*M_PI*M_PI*detA_inv)*exp(0.5*(BAB-C));


            if (idx_update==1) { // for update pi_k            
                hess_pi_k[numelipsoids*k+l] = (double) gg + (k==l)*lambda_pi;
                hess_pi_k[numelipsoids*l+k] = (double) gg + (k==l)*lambda_pi;
                continue;
            }


            // diagonalize A using dsyevd
            // eigenvector should be normalized (returned matrix A should be orthonormal)
            // [A[0],A[3],A[6]] from returned matrix A is an eigenvector            
            A[0]=A_11; A[1]=A_12; A[2] = A_13;
            A[3]=A_12; A[4]=A_22; A[5] = A_23;
            A[6]=A_13; A[7]=A_23; A[8] = A_33;
            dsyevd("V","U",&n,A,&n,w,work,&lwork,iwork,&liwork,&info);
            
            double U[3][3] = {{A[0], A[3], A[6]},
                              {A[1], A[4], A[7]},
                              {A[2], A[5], A[8]}};
            double D[3] = {1/w[0],1/w[1],1/w[2]}; // inverse of eigenvalue


            // calc W_kl

            // S_k^-1
            double Sk[3][3] = {{A_k_11,A_k_12,A_k_13},
                               {A_k_12,A_k_22,A_k_23},
                               {A_k_13,A_k_23,A_k_33}};
            double Sl[3][3] = {{A_l_11,A_l_12,A_l_13},
                               {A_l_12,A_l_22,A_l_23},
                               {A_l_13,A_l_23,A_l_33}};

            // Rk = Sk*U
            double Rk[3][3] = {{Sk[0][0]*U[0][0] + Sk[0][1]*U[1][0] + Sk[0][2]*U[2][0],   // Rk_11
                                Sk[0][0]*U[0][1] + Sk[0][1]*U[1][1] + Sk[0][2]*U[2][1],   // Rk_12
                                Sk[0][0]*U[0][2] + Sk[0][1]*U[1][2] + Sk[0][2]*U[2][2]},  // Rk_13
                               {Sk[1][0]*U[0][0] + Sk[1][1]*U[1][0] + Sk[1][2]*U[2][0],   // Rk_21
                                Sk[1][0]*U[0][1] + Sk[1][1]*U[1][1] + Sk[1][2]*U[2][1],   // Rk_22
                                Sk[1][0]*U[0][2] + Sk[1][1]*U[1][2] + Sk[1][2]*U[2][2]},  // Rk_23
                               {Sk[2][0]*U[0][0] + Sk[2][1]*U[1][0] + Sk[2][2]*U[2][0],   // Rk_31
                                Sk[2][0]*U[0][1] + Sk[2][1]*U[1][1] + Sk[2][2]*U[2][1],   // Rk_32
                                Sk[2][0]*U[0][2] + Sk[2][1]*U[1][2] + Sk[2][2]*U[2][2]}}; // Rk_33
                              
            // Rl = Sl*U
            double Rl[3][3] = {{Sl[0][0]*U[0][0] + Sl[0][1]*U[1][0] + Sl[0][2]*U[2][0],   // Rl_11
                                Sl[0][0]*U[0][1] + Sl[0][1]*U[1][1] + Sl[0][2]*U[2][1],   // Rl_12
                                Sl[0][0]*U[0][2] + Sl[0][1]*U[1][2] + Sl[0][2]*U[2][2]},  // Rl_13
                               {Sl[1][0]*U[0][0] + Sl[1][1]*U[1][0] + Sl[1][2]*U[2][0],   // Rl_21
                                Sl[1][0]*U[0][1] + Sl[1][1]*U[1][1] + Sl[1][2]*U[2][1],   // Rl_22
                                Sl[1][0]*U[0][2] + Sl[1][1]*U[1][2] + Sl[1][2]*U[2][2]},  // Rl_23
                               {Sl[2][0]*U[0][0] + Sl[2][1]*U[1][0] + Sl[2][2]*U[2][0],   // Rl_31
                                Sl[2][0]*U[0][1] + Sl[2][1]*U[1][1] + Sl[2][2]*U[2][1],   // Rl_32
                                Sl[2][0]*U[0][2] + Sl[2][1]*U[1][2] + Sl[2][2]*U[2][2]}}; // Rl_33


            //// if k=l, A^-1*B = mu_k, then tk=vk=0 and Rvk=Rvl=0.
            double Rvk[3] = {0,0,0};
            double Rvl[3] = {0,0,0};
            if (k!=l) {
                // tk = A^-1*B-mu_k 
                double tk[3] = {AB[0]-xc_k,
                                AB[1]-yc_k,
                                AB[2]-zc_k};

                // tl = A^-1*B-mu_l
                double tl[3] = {AB[0]-xc_l,
                                AB[1]-yc_l,
                                AB[2]-zc_l};

                // Rv_k = Rk*vk = Sk*tk
                Rvk[0] = Sk[0][0]*tk[0] + Sk[0][1]*tk[1] + Sk[0][2]*tk[2];
                Rvk[1] = Sk[1][0]*tk[0] + Sk[1][1]*tk[1] + Sk[1][2]*tk[2];
                Rvk[2] = Sk[2][0]*tk[0] + Sk[2][1]*tk[1] + Sk[2][2]*tk[2];

                // Rv_l = Rl*vl = Sl*tl
                Rvl[0] = Sl[0][0]*tl[0] + Sl[0][1]*tl[1] + Sl[0][2]*tl[2];
                Rvl[1] = Sl[1][0]*tl[0] + Sl[1][1]*tl[1] + Sl[1][2]*tl[2];
                Rvl[2] = Sl[2][0]*tl[0] + Sl[2][1]*tl[1] + Sl[2][2]*tl[2];
            } 

            
            // W_kl_phi; [1], symmetric
            double W_phi = gg;

            // W_kl_li; [i], assymetric
            /* not used in the alternative optimization method
            double W_ki[3] = {Rvk[0]*gg,
                              Rvk[1]*gg,
                              Rvk[2]*gg};
            double W_li[3] = {Rvl[0]*gg,
                              Rvl[1]*gg,
                              Rvl[2]*gg};
            */

            // W_kl_kilj; [i,j], symmetric
            double W_kilj[3][3];
            for (unsigned int i=0; i<3; i++) {
                for (unsigned int j=0; j<3; j++) {
                    double p0 = Rvk[i]*Rvl[j];
                    double p2d = 0.0;
                    for (unsigned int p=0; p<3; p++) {
                        p2d += Rk[i][p]*Rl[j][p]*D[p];
                    }
                    W_kilj[i][j] = (p0 + p2d)*gg;
                }
            }
            /* not used in the alternative optimization method
            double W_kikj[3][3];
            for (unsigned int i=0; i<3; i++) {
                for (unsigned int j=0; j<3; j++) {
                    double p0 = Rvk[i]*Rvk[j];
                    double p2d = 0.0;
                    for (unsigned int p=0; p<3; p++) {
                        p2d += Rk[i][p]*Rk[j][p]*D[p];
                    }
                    W_kikj[i][j] = (p0 + p2d)*gg;
                }
            }
            double W_lilj[3][3];
            for (unsigned int i=0; i<3; i++) {
                for (unsigned int j=0; j<3; j++) {
                    double p0 = Rvl[i]*Rvl[j];
                    double p2d = 0.0;
                    for (unsigned int p=0; p<3; p++) {
                        p2d += Rl[i][p]*Rl[j][p]*D[p];
                    }
                    W_lilj[i][j] = (p0 + p2d)*gg;
                }
            }
            */

            
            // W_kl_khlilj; [h,i,j], assymmetric
            /* not used in the alternative optimization method
            double W_khlilj[3][3][3];
            for (unsigned int h=0; h<3; h++) {
                for (unsigned int i=0; i<3; i++) {
                    for (unsigned int j=0; j<3; j++) {
                        double p0 = Rvk[h]*Rvl[i]*Rvl[j];
                        double p2d = 0.0;
                        for (unsigned int p=0; p<3; p++) {
                            p2d += (  Rk[h][p] * Rl[i][p] * Rvl[j]
                                    + Rk[h][p] * Rvl[i]   * Rl[j][p]
                                    + Rvk[h]   * Rl[i][p] * Rl[j][p])*D[p];
                        }
                        W_khlilj[h][i][j] = (p0 + p2d)*gg;
                    }
                }
            }
            double W_lhkikj[3][3][3];
            for (unsigned int h=0; h<3; h++) {
                for (unsigned int i=0; i<3; i++) {
                    for (unsigned int j=0; j<3; j++) {
                        double p0 = Rvl[h]*Rvk[i]*Rvk[j];
                        double p2d = 0.0;
                        for (unsigned int p=0; p<3; p++) {
                            p2d += (  Rl[h][p] * Rk[i][p] * Rvk[j]
                                    + Rl[h][p] * Rvk[i]   * Rk[j][p]
                                    + Rvl[h]   * Rk[i][p] * Rk[j][p])*D[p];
                        }
                        W_lhkikj[h][i][j] = (p0 + p2d)*gg;
                    }
                }
            }
            */

            // W_kl_kpkqlilj; [p,q,i,j], symmetric
            double W_kpkqlilj[3][3][3][3];
            if (idx_update==3) // for update sigma
            for (unsigned int p=0; p<3; p++) {
                for (unsigned int q=0; q<3; q++) {
                    for (unsigned int i=0; i<3; i++) {
                        for (unsigned int j=0; j<3; j++) {
                            double p0 = Rvk[p]*Rvk[q]*Rvl[i]*Rvl[j];
                            double p2d = 0.0;
                            for (unsigned int h=0; h<3; h++) {
                                p2d += (  Rk[p][h] * Rk[q][h] * Rvl[i]   * Rvl[j]
                                        + Rk[p][h] * Rvk[q]   * Rl[i][h] * Rvl[j]
                                        + Rk[p][h] * Rvk[q]   * Rvl[i]   * Rl[j][h]
                                        + Rvk[p]   * Rk[q][h] * Rl[i][h] * Rvl[j]
                                        + Rvk[p]   * Rk[q][h] * Rvl[i]   * Rl[j][h]
                                        + Rvk[p]   * Rvk[q]   * Rl[i][h] * Rl[j][h])*D[h];
                            }
                            double p22dd = 0.0;
                            for (unsigned int r=0; r<3; r++) {
                                for (unsigned int s=r+1; s<3; s++) {
                                    p22dd += (  Rk[p][r]*Rk[q][r]*Rl[i][s]*Rl[j][s]
                                              + Rk[p][r]*Rk[q][s]*Rl[i][r]*Rl[j][s]
                                              + Rk[p][r]*Rk[q][s]*Rl[i][s]*Rl[j][r]
                                              + Rk[p][s]*Rk[q][r]*Rl[i][r]*Rl[j][s]
                                              + Rk[p][s]*Rk[q][r]*Rl[i][s]*Rl[j][r]
                                              + Rk[p][s]*Rk[q][s]*Rl[i][r]*Rl[j][r])*D[r]*D[s];
                                }
                            }
                            double p4dd = 0.0;
                            for (unsigned int h=0; h<3; h++) {
                                p4dd += Rk[p][h]*Rk[q][h]*Rl[i][h]*Rl[j][h]*D[h]*D[h];
                            }
                            W_kpkqlilj[p][q][i][j] = (3*p4dd+p22dd+p2d+p0)*gg;
                        }
                    }
                }
            }

            // for update mu
            if (idx_update==2) {
                unsigned int offset1 = NUM_PARAMS_MU*NUM_PARAMS_MU*numelipsoids*k+NUM_PARAMS_MU*l;
                unsigned int offset2 = NUM_PARAMS_MU*NUM_PARAMS_MU*numelipsoids*l+NUM_PARAMS_MU*k;
                for (unsigned int i=0; i<3; i++) { // mu_{k,i}
                    for (unsigned int j=0; j<3; j++) { // mu_{l,j}
                        double hess = 2*pi_k*( (k==l)*(Sk[i][j]*V_phi - V_ij[i][j])
                                        + pi_l*W_kilj[i][j]) * rss_orig_inv;
                        hess_mu[offset1+NUM_PARAMS_MU*numelipsoids*i+j] += hess;
                        hess_mu[offset2+NUM_PARAMS_MU*numelipsoids*j+i] += hess;
                    }
                }
            }

            // for update sigma
            if (idx_update==3) {            
                unsigned int offset1 = NUM_PARAMS_SIGMA*NUM_PARAMS_SIGMA*numelipsoids*k+NUM_PARAMS_SIGMA*l;
                unsigned int offset2 = NUM_PARAMS_SIGMA*NUM_PARAMS_SIGMA*numelipsoids*l+NUM_PARAMS_SIGMA*k;
                for (unsigned int p=0; p<3; p++) {
                    for (unsigned int q=p; q<3; q++) {
                        unsigned int idx1 = s2i[p][q][0][0]; // Sigma_{k,pq}
                        for (unsigned int i=0; i<3; i++) {
                            for (unsigned int j=i; j<3; j++) {
                                unsigned int idx2 = s2i[i][j][0][0]; // Sigma_{l,ij}
                                double hess = 0.5*pi_k*(2-(i==j))*(2-(p==q)) * rss_orig_inv
                                             *( (k==l)*( Sk[i][q]*V_ij[j][p] + Sk[j][q]*V_ij[i][p]
                                                       + Sk[i][p]*V_ij[j][q] + Sk[j][p]*V_ij[i][q]
                                                       - V_ijpq[i][j][p][q] )
                                                + pi_l*W_kpkqlilj[p][q][i][j] );
                                if ( (k==l) & (idx1==idx2) ) {
                                    hess += 2*lambda_sigma*(2-(p==q))*rss_orig_inv; // L2 regularization
                                }
                                hess_sigma[offset1+NUM_PARAMS_SIGMA*numelipsoids*idx1+idx2] = hess;
                                hess_sigma[offset2+NUM_PARAMS_SIGMA*numelipsoids*idx2+idx1] = hess;
                            }
                        }
                    }
                }
            }
            
        } // end of loop for l

        // for update mu and sigma
        if (idx_update>=2) {
            // delete matrices
            delete[] work;
            delete[] iwork;
            delete[] w;
            delete[] A;
        }

    } // end of loop for k
}// end of function calc_W


void update_params(int idx_update) {
    // int idx_update; // specify what should be updatated. -1:calc_ofv, 0:imsynth, 1:pi_k, 2:mu, 3:sigma    
    
    #ifdef MY_ASSERT
        unsigned __int64 t_start=0,t_calcofv=0,t_calcw=0,t_LM=0;
        t_start = __rdtsc();
    #endif

    calc_ofv(idx_update); // update V_k and score

    #ifdef MY_ASSERT
        t_calcofv += __rdtsc() - t_start;
        t_start = __rdtsc();
    #endif

    calc_W(idx_update); // update W_kl and set gradient and hessian

    #ifdef MY_ASSERT
        t_calcw += __rdtsc() - t_start;
        t_start = __rdtsc();
    #endif

    unsigned int n, offset;
    double *hess, *grad, *hess_orig, *grad_orig, *delta_oldstep;
    double lmp;
    
    switch (idx_update) {
    case 1: // for update pi_k
        n = NUM_PARAMS_PI_K;
        offset = 0;
        hess = hess_pi_k;
        grad = grad_pi_k;
        break;
    case 2: // for update mu
        n = NUM_PARAMS_MU;
        offset = 1;
        hess = hess_mu;
        grad = grad_mu;
        lmp = lmp_mu;
        break;
    case 3: // for update sigma
        n = NUM_PARAMS_SIGMA;
        offset = 4;
        hess = hess_sigma;
        grad = grad_sigma;
        lmp = lmp_sigma;
        break;
    }

    
    if (idx_update>=2) {
        // isolate grad and hessian because these values are broken by dgelsy
        unsigned int n1 = n*numelipsoids;
        grad_orig = new double[n1];
        for (unsigned int c=0; c<n1; c++) {
            grad_orig[c] = grad[c];
        }
        unsigned int n2 = n1*n1;
        hess_orig = new double[n2];
        for (unsigned int c=0; c<n2; c++) {
            hess_orig[c] = hess[c];
        }

        // Isolate old params
        for (unsigned int c=0; c<NUM_PARAMS_INPUT*numelipsoids; c++) {
            params_old[c] = params[c]; // isolate old param
        }

        // buffer for delta_x
        delta_oldstep = new double[n1];
    }

    // for initializing dgelsy (mldivide)
    // query and allocate the optimal workspace for dgelsy
    lapack_int nlhs = n*numelipsoids;
    lapack_int nrhs = 1;
    lapack_int lda = nlhs;
    lapack_int ldb = nlhs; 
    lapack_int info = 0;
    lapack_int rank = 0;
    lapack_int lwork = -1;
    lapack_int *jpvt = new lapack_int[nlhs];
    //double rcond = -1.0; // rcond<0 means using machine epsillon
    double rcond = 100*DBL_EPSILON;
    double* work = new double[1];        
    ////// dsysv("Upper", &nlhs, &nrhs, hess, &lda, ipiv, grad, &ldb, work, &lwork, &info);
    //// dgelsd(&nlhs,&nlhs,&nrhs,hess,&lda,grad,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);    
    for(int c=0;c<nlhs;c++) { jpvt[c]=0; } // initialize jpvt
    dgelsy(&nlhs,&nlhs,&nrhs,hess,&lda,grad,&ldb,jpvt,&rcond,&rank,work,&lwork,&info);   
    lwork = (lapack_int)work[0];
    delete[] work;
    work = new double[__max(1,lwork)];


    // initialize dsyevd for diagonalizing Sigma
    lapack_int ndim_s, lwork_s, liwork_s, info_s, *iwork_s;
    double *work_s, *Amat_s, *evalue_s;
    if (idx_update==3 && fixvol==1) {
        ndim_s = 3;
        lwork_s = -1;
        liwork_s = -1;
        work_s   = new double[1];
        iwork_s  = new lapack_int[1];
        Amat_s   = new double[ndim_s*ndim_s];    
        evalue_s = new double[ndim_s];    
        Amat_s[0]=10; Amat_s[1]= 0; Amat_s[2]= 0;
        Amat_s[3]= 0; Amat_s[4]=10; Amat_s[5]= 0;
        Amat_s[6]= 0; Amat_s[7]= 0; Amat_s[8]=10;
        dsyevd("V","U",&ndim_s,Amat_s,&ndim_s,evalue_s,work_s,&lwork_s,iwork_s,&liwork_s,&info_s);
        lwork_s = (lapack_int)(work_s[0]);
        liwork_s = iwork_s[0];
        delete[] work_s;
        delete[] iwork_s;
        work_s = new double[__max(1,lwork_s)];
        iwork_s = new lapack_int[__max(1,liwork_s)];
    }


    // for update pi_k     
    if (idx_update==1) {
        ////// use dsysv because hess_pi_k may not be a positive definite matrix
        ////// dposv("Upper", &nlhs, &nrhs, Amat, &lda, Bmat, &ldb, &info);
        ////// dsysv("Upper", &nlhs, &nrhs, hess, &lda, ipiv, grad, &ldb, work, &lwork, &info);
        //// use dgelsd because hess may be a rank-deficit matrix
        ////dgelsd(&nlhs,&nlhs,&nrhs,hess,&lda,grad,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);   
        // use dgelsy instead of dgelsd because of performance reason
        for(int c=0;c<nlhs;c++) { jpvt[c]=0; } // initialize jpvt
        dgelsy(&nlhs,&nlhs,&nrhs,hess,&lda,grad,&ldb,jpvt,&rcond,&rank,work,&lwork,&info);

        delete[] jpvt;
        delete[] work;

        // error occured
        if (info!=0) {
            mexPrintf("Error in dgelsy for update pi_k. Please check illegal ROIs. (idx:%d, info: %d.) \n", idx_update, info);
            return;
        }

        // update pi_k
        for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
            params[NUM_PARAMS_INPUT*celipsoid+offset] = grad[n*celipsoid];
        }
        return;
    } 

    // for update mu and sigma
    if (idx_update>=2) {

        double score_start = score;// isolate old params and score
        double score_oldstep = DBL_MAX;
        double lmp_oldstep;

        // trial for larger step
        lmp = lmp / (lmd*lmd);

        for (unsigned int cloop=0; cloop<lm_maxiter; cloop++) {
            // update lmp
            lmp *= lmd;
            if (!_finite(lmp)) { break; } 

            // set grad and hess with new LM parameter
            unsigned int n1 = n*numelipsoids;
            for (unsigned int c1=0; c1<n1; c1++) {
                grad[c1] = grad_orig[c1];
                for (unsigned int c2=0; c2<n1; c2++) {
                    hess[n1*c1+c2] = hess_orig[n1*c1+c2] * (1.0 + lmp*(c1==c2));
                }
            }

            // solve the equation and obtain delta
            ////// use dsysv because hess_pi_k may not be a positive definite matrix
            ////// dsysv("Upper", &nlhs, &nrhs, hess, &lda, ipiv, grad, &ldb, work, &lwork, &info);
            //// use dgelsd because hess may be a rank-deficit matrix
            //dgelsd(&nlhs,&nlhs,&nrhs,hess,&lda,grad,&ldb,s,&rcond,&rank,work,&lwork,iwork,&info);
            // use dgelsy instead of dgelsd because of performance reason
            for(int c=0;c<nlhs;c++) { jpvt[c]=0; } // initialize jpvt
            dgelsy(&nlhs,&nlhs,&nrhs,hess,&lda,grad,&ldb,jpvt,&rcond,&rank,work,&lwork,&info);

            if (info!=0) { // error occured in dgelsd; continue to next loop
                mexPrintf("Error in dgelsy for update mu & sigma: Try larger LM parameter. (cloop:%d, idx:%d, info: %d, lmp: %f.) \n", cloop, idx_update, info, lmp);
                continue;
            }

            // Calc ofv for new params
            for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
                if (idx_update==3) { // cancel update for sigma that produce negative determinant
                    double old_S_11 = params_old[NUM_PARAMS_INPUT*celipsoid+4];
                    double old_S_12 = params_old[NUM_PARAMS_INPUT*celipsoid+5];
                    double old_S_13 = params_old[NUM_PARAMS_INPUT*celipsoid+6];
                    double old_S_22 = params_old[NUM_PARAMS_INPUT*celipsoid+7];
                    double old_S_23 = params_old[NUM_PARAMS_INPUT*celipsoid+8];
                    double old_S_33 = params_old[NUM_PARAMS_INPUT*celipsoid+9];
                    double tmp_S_11 = old_S_11 + grad[n*celipsoid+0];
                    double tmp_S_12 = old_S_12 + grad[n*celipsoid+1];
                    double tmp_S_13 = old_S_13 + grad[n*celipsoid+2];
                    double tmp_S_22 = old_S_22 + grad[n*celipsoid+3];
                    double tmp_S_23 = old_S_23 + grad[n*celipsoid+4];
                    double tmp_S_33 = old_S_33 + grad[n*celipsoid+5];
                    double detS = (   tmp_S_11*tmp_S_22*tmp_S_33 
                                  + 2*tmp_S_12*tmp_S_13*tmp_S_23
                                  -   tmp_S_11*tmp_S_23*tmp_S_23 
                                  -   tmp_S_12*tmp_S_12*tmp_S_33
                                  -   tmp_S_13*tmp_S_22*tmp_S_13);
                    if (detS<0 || tmp_S_11<0 || tmp_S_22<0 || tmp_S_33<0) {
                        for (unsigned int c2=0; c2<6; c2++) { grad[n*celipsoid+c2] = 0.0; }
                    }

                    if (fixvol==1) { // fix eigenvalue of sigma

                        // diagonalize S using dsyevd
                        // eigenvector should be normalized (returned matrix Amat should be orthonormal)
                        // [Amat[0],Amat[3],Amat[6]] from returned matrix Amat is an eigenvector            
                        Amat_s[0]=old_S_11; Amat_s[1]=old_S_12; Amat_s[2]=old_S_13;
                        Amat_s[3]=old_S_12; Amat_s[4]=old_S_22; Amat_s[5]=old_S_23;
                        Amat_s[6]=old_S_13; Amat_s[7]=old_S_23; Amat_s[8]=old_S_33;
                        dsyevd("V","U",&ndim_s,Amat_s,&ndim_s,evalue_s,work_s,&lwork_s,iwork_s,&liwork_s,&info_s);

                        double old_U[3][3] = {{Amat_s[0], Amat_s[3], Amat_s[6]},
                                              {Amat_s[1], Amat_s[4], Amat_s[7]},
                                              {Amat_s[2], Amat_s[5], Amat_s[8]}};
                        double old_D[3] = {evalue_s[0],evalue_s[1],evalue_s[2]}; //  eigenvalue

                        Amat_s[0]=tmp_S_11; Amat_s[1]=tmp_S_12; Amat_s[2]=tmp_S_13;
                        Amat_s[3]=tmp_S_12; Amat_s[4]=tmp_S_22; Amat_s[5]=tmp_S_23;
                        Amat_s[6]=tmp_S_13; Amat_s[7]=tmp_S_23; Amat_s[8]=tmp_S_33;
                        dsyevd("V","U",&ndim_s,Amat_s,&ndim_s,evalue_s,work_s,&lwork_s,iwork_s,&liwork_s,&info_s);

                        double tmp_U[3][3] = {{Amat_s[0], Amat_s[3], Amat_s[6]},
                                              {Amat_s[1], Amat_s[4], Amat_s[7]},
                                              {Amat_s[2], Amat_s[5], Amat_s[8]}};
                        double tmp_D[3] = {evalue_s[0],evalue_s[1],evalue_s[2]}; //  eigenvalue

                        //tmp_U * old_D * tmpU'
                        double tmp2_S_11 =  tmp_U[0][0]*tmp_U[0][0]*old_D[0]
                                          + tmp_U[0][1]*tmp_U[0][1]*old_D[1] 
                                          + tmp_U[0][2]*tmp_U[0][2]*old_D[2];
                        double tmp2_S_12 =  tmp_U[0][0]*tmp_U[1][0]*old_D[0]
                                          + tmp_U[0][1]*tmp_U[1][1]*old_D[1] 
                                          + tmp_U[0][2]*tmp_U[1][2]*old_D[2];
                        double tmp2_S_13 =  tmp_U[0][0]*tmp_U[2][0]*old_D[0]
                                          + tmp_U[0][1]*tmp_U[2][1]*old_D[1] 
                                          + tmp_U[0][2]*tmp_U[2][2]*old_D[2];
                        double tmp2_S_22 =  tmp_U[1][0]*tmp_U[1][0]*old_D[0]
                                          + tmp_U[1][1]*tmp_U[1][1]*old_D[1] 
                                          + tmp_U[1][2]*tmp_U[1][2]*old_D[2];
                        double tmp2_S_23 =  tmp_U[1][0]*tmp_U[2][0]*old_D[0]
                                          + tmp_U[1][1]*tmp_U[2][1]*old_D[1] 
                                          + tmp_U[1][2]*tmp_U[2][2]*old_D[2];
                        double tmp2_S_33 =  tmp_U[2][0]*tmp_U[2][0]*old_D[0]
                                          + tmp_U[2][1]*tmp_U[2][1]*old_D[1] 
                                          + tmp_U[2][2]*tmp_U[2][2]*old_D[2];

                        grad[n*celipsoid+0] = tmp2_S_11 - old_S_11;
                        grad[n*celipsoid+1] = tmp2_S_12 - old_S_12;
                        grad[n*celipsoid+2] = tmp2_S_13 - old_S_13;                        
                        grad[n*celipsoid+3] = tmp2_S_22 - old_S_22;
                        grad[n*celipsoid+4] = tmp2_S_23 - old_S_23;
                        grad[n*celipsoid+5] = tmp2_S_33 - old_S_33;

                    }
                }
                for (unsigned int c2=0; c2<n; c2++) {
                    unsigned int idx = NUM_PARAMS_INPUT*celipsoid+offset+c2;
                    params[idx] = params_old[idx] + grad[n*celipsoid+c2];
                }
            }
            setparam();
            calc_ofv(-1); // calc score

            // Update parameter based on following rule;
            // If successful step have not been obtained yet, save current state and do next.
            // If at least either of score or score_oldstep is better than score_start,
            // use the better one.
            // If neither, do next.
            if (score_oldstep>=DBL_MAX) { // no successful step yet
                score_oldstep = score;
                lmp_oldstep = lmp;
                for (unsigned int c=0; c<n*numelipsoids; c++) {
                    delta_oldstep[c] = grad[c];
                }
            } else if (score_oldstep<score_start || score<score_start) {
                if (score_oldstep < score) { // old step is better; set params of oldstep
                    lmp = lmp_oldstep;
                    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
                        for (unsigned int c2=0; c2<n; c2++) {
                            unsigned int idx = NUM_PARAMS_INPUT*celipsoid+offset+c2;
                            params[idx] = params_old[idx] + delta_oldstep[n*celipsoid+c2];
                        }
                    }
                }
        
                #ifdef MY_ASSERT
                    mexPrintf("update_params.LM_loop_counter: %d\n", cloop);
                    mexPrintf("update_params.LM_current_lambda: %e\n", lmp);
                #endif
                
                break;
            }
        } // end of for loop

        delete[] jpvt;
        delete[] work;
        delete[] grad_orig;
        delete[] hess_orig;
        delete[] delta_oldstep;
        if(idx_update==3 && fixvol==1) {
            delete[] work_s;
            delete[] iwork_s;
            delete[] evalue_s;
            delete[] Amat_s;
        }

    } // end of if for update mu and sigma

    if (idx_update==2) { lmp_mu = lmp; }
    if (idx_update==3) { lmp_sigma = lmp; }

    #ifdef MY_ASSERT
        t_LM += __rdtsc() - t_start;
        mexPrintf("update_params.calcOFV: %I64d[clocks]\n", t_calcofv);
        mexPrintf("update_params.calcW: %I64d[clocks]\n", t_calcw);
        mexPrintf("update_params.LM: %I64d[clocks]\n", t_LM);
    #endif

}


void optimize(){

    #ifdef MY_ASSERT
        unsigned __int64 t_start=0,t_setparam=0,t_estep=0,t_m1step=0,t_m2step=0;
    #endif

    double oldscore = DBL_MAX;
    finiter = 0;

    __int64 verbose_t = 0;
    if (verbose==1) {
        for (unsigned int citer=0; citer<maxiter; citer++) {
            verbose_etime[citer] = 0;
            verbose_score[citer] = 0.0;
        }
    }

    for (unsigned int citer=0; citer<maxiter; citer++) { 

        mexPrintf("citer = %d,\t", citer);
        mexEvalString("drawnow");
        
        if (verbose==1) {
            verbose_t = __rdtsc();
        }

        if (!fixmu) {

            #ifdef MY_ASSERT
                t_start = __rdtsc();
            #endif

            setparam();
            
            #ifdef MY_ASSERT
                t_setparam += __rdtsc() - t_start;
                t_start = __rdtsc();
            #endif

            if (!fixpik) {      

                update_params(1); // update pi_k
            
                #ifdef MY_ASSERT
                t_estep += __rdtsc() - t_start;
                t_start = __rdtsc();
                #endif

                setparam();

                #ifdef MY_ASSERT
                    t_setparam += __rdtsc() - t_start;
                    t_start = __rdtsc();
                #endif

            }

            update_params(2); // update mu

            #ifdef MY_ASSERT
                t_m1step += __rdtsc() - t_start;
            #endif
        }

        #ifdef MY_ASSERT
            t_start = __rdtsc();
        #endif
        
        setparam();

        #ifdef MY_ASSERT
            t_setparam += __rdtsc() - t_start;
            t_start = __rdtsc();
        #endif
        
        if (!fixpik) {
            update_params(1); // update pi_k

            #ifdef MY_ASSERT
            t_estep += __rdtsc() - t_start;
            t_start = __rdtsc();
            #endif

            setparam();

            #ifdef MY_ASSERT
               t_setparam += __rdtsc() - t_start;
               t_start = __rdtsc();
            #endif

        }

        update_params(3); // update sigma;

        #ifdef MY_ASSERT
            t_m2step += __rdtsc() - t_start;
        #endif

        mexPrintf("score = %15.12f\n", score);


        if (verbose==1) {
            verbose_etime[citer] = __rdtsc() - verbose_t;
            verbose_score[citer] = score;
        }
        finiter = citer;


        if ((oldscore-score)<tol) {break;}
        oldscore = score;
        
    }    

    
    #ifdef MY_ASSERT
        mexPrintf("optimize.setparam: %I64d[clocks]\n", t_setparam);
        mexPrintf("optimize.update_pi_k: %I64d[clocks]\n", t_estep);
        mexPrintf("optimize.update_mu: %I64d[clocks]\n", t_m1step);
        mexPrintf("optimize.update_sigma: %I64d[clocks]\n", t_m2step);
    #endif

}


void setsynth2(double *psynth) {

    /* update synthetic image (im_synth) */    
    setparam();
    calc_ofv(0);

    /* set data */
    unsigned int cvxyz, cout;

    for (unsigned int cy=0; cy<numy; cy++) {

        // for elements without padding
        for (unsigned int cvx=0; cvx<numVectorX-1; cvx++) {
            for (unsigned int cz=0; cz<numz; cz++) {
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                cout = numx*numy*cz + numx*cy + SIMD_VECTOR_LENGTH*cvx;
                for (unsigned int ce=0; ce<SIMD_VECTOR_LENGTH; ce++) {
                    psynth[ cout + ce ] = (double) im_synth[ cvxyz*SIMD_VECTOR_LENGTH + ce ];
                }
            }
        } // end of loop for elements without padding

        // for elements with padding
        unsigned int cvx = numVectorX - 1; {
            for (unsigned int cz=0; cz<numz; cz++) {
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                cout = numx*numy*cz + numx*cy + SIMD_VECTOR_LENGTH*cvx;
                for (unsigned int ce=0; ce<SIMD_VECTOR_LENGTH-padlength; ce++) {
                    psynth[ cout + ce ] = (double) im_synth[ cvxyz*SIMD_VECTOR_LENGTH + ce ];
                }
            }
        } // end of loop for elements with padding
    } // end of loop for y
}
