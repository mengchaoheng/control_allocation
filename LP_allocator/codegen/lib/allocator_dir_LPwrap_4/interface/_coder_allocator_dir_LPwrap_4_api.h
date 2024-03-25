/*
 * Prerelease License - for engineering feedback and testing purposes
 * only. Not for sale.
 * File: _coder_allocator_dir_LPwrap_4_api.h
 *
 * MATLAB Coder version            : 24.1
 * C/C++ source code generated on  : 2024-03-24 20:52:32
 */

#ifndef _CODER_ALLOCATOR_DIR_LPWRAP_4_API_H
#define _CODER_ALLOCATOR_DIR_LPWRAP_4_API_H

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
void allocator_dir_LPwrap_4(real32_T B[12], real32_T v[3], real32_T umin[4],
                            real32_T umax[4], real32_T u[4], real32_T *z,
                            uint16_T *iters);

void allocator_dir_LPwrap_4_api(const mxArray *const prhs[4], int32_T nlhs,
                                const mxArray *plhs[3]);

void allocator_dir_LPwrap_4_atexit(void);

void allocator_dir_LPwrap_4_initialize(void);

void allocator_dir_LPwrap_4_terminate(void);

void allocator_dir_LPwrap_4_xil_shutdown(void);

void allocator_dir_LPwrap_4_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for _coder_allocator_dir_LPwrap_4_api.h
 *
 * [EOF]
 */
