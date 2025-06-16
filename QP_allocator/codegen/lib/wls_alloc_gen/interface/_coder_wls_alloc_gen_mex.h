/*
 * File: _coder_wls_alloc_gen_mex.h
 *
 * MATLAB Coder version            : 24.1
 * C/C++ source code generated on  : 2025-06-16 11:03:22
 */

#ifndef _CODER_WLS_ALLOC_GEN_MEX_H
#define _CODER_WLS_ALLOC_GEN_MEX_H

/* Include Files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
MEXFUNCTION_LINKAGE void mexFunction(int32_T nlhs, mxArray *plhs[],
                                     int32_T nrhs, const mxArray *prhs[]);

emlrtCTX mexFunctionCreateRootTLS(void);

void unsafe_wls_alloc_gen_mexFunction(int32_T nlhs, mxArray *plhs[1],
                                      int32_T nrhs, const mxArray *prhs[12]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for _coder_wls_alloc_gen_mex.h
 *
 * [EOF]
 */
