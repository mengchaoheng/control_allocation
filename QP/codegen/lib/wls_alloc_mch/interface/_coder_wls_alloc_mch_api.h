/*
 * File: _coder_wls_alloc_mch_api.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 08-Aug-2020 14:31:14
 */

#ifndef _CODER_WLS_ALLOC_MCH_API_H
#define _CODER_WLS_ALLOC_MCH_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_wls_alloc_mch_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void wls_alloc_mch(real_T v[3], real_T u[8]);
extern void wls_alloc_mch_api(const mxArray * const prhs[2], int32_T nlhs, const
  mxArray *plhs[1]);
extern void wls_alloc_mch_atexit(void);
extern void wls_alloc_mch_initialize(void);
extern void wls_alloc_mch_terminate(void);
extern void wls_alloc_mch_xil_terminate(void);

#endif

/*
 * File trailer for _coder_wls_alloc_mch_api.h
 *
 * [EOF]
 */
