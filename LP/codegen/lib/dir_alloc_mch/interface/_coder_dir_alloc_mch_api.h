/*
 * File: _coder_dir_alloc_mch_api.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Nov-2019 20:09:11
 */

#ifndef _CODER_DIR_ALLOC_MCH_API_H
#define _CODER_DIR_ALLOC_MCH_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_dir_alloc_mch_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void dir_alloc_mch(real_T v[3], real_T umin[4], real_T umax[4], real_T u
  [4]);
extern void dir_alloc_mch_api(const mxArray * const prhs[3], int32_T nlhs, const
  mxArray *plhs[1]);
extern void dir_alloc_mch_atexit(void);
extern void dir_alloc_mch_initialize(void);
extern void dir_alloc_mch_terminate(void);
extern void dir_alloc_mch_xil_terminate(void);

#endif

/*
 * File trailer for _coder_dir_alloc_mch_api.h
 *
 * [EOF]
 */
