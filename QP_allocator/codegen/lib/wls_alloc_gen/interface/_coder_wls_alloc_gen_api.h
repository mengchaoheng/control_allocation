/*
 * File: _coder_wls_alloc_gen_api.h
 *
 * MATLAB Coder version            : 24.1
 * C/C++ source code generated on  : 2025-06-16 11:03:22
 */

#ifndef _CODER_WLS_ALLOC_GEN_API_H
#define _CODER_WLS_ALLOC_GEN_API_H

/* Include Files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <string.h>

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void wls_alloc_gen(real32_T B[12], real32_T v[3], real32_T umin[4],
                   real32_T umax[4], real32_T Wv[9], real32_T Wu[16],
                   real32_T ud[4], real32_T gam, real32_T u[4], real32_T W[4],
                   real32_T imax, real32_T m);

void wls_alloc_gen_api(const mxArray *const prhs[12], const mxArray **plhs);

void wls_alloc_gen_atexit(void);

void wls_alloc_gen_initialize(void);

void wls_alloc_gen_terminate(void);

void wls_alloc_gen_xil_shutdown(void);

void wls_alloc_gen_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for _coder_wls_alloc_gen_api.h
 *
 * [EOF]
 */
