/*
 * Prerelease License - for engineering feedback and testing purposes
 * only. Not for sale.
 * File: _coder_allocator_dir_LPwrap_4_mex.h
 *
 * MATLAB Coder version            : 24.1
 * C/C++ source code generated on  : 2024-03-24 20:52:32
 */

#ifndef _CODER_ALLOCATOR_DIR_LPWRAP_4_MEX_H
#define _CODER_ALLOCATOR_DIR_LPWRAP_4_MEX_H

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

void unsafe_allocator_dir_LPwrap_4_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                               int32_T nrhs,
                                               const mxArray *prhs[4]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for _coder_allocator_dir_LPwrap_4_mex.h
 *
 * [EOF]
 */
