/*
 * File: _coder_qp_ca_mch_api.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 17-Aug-2020 18:49:57
 */

#ifndef _CODER_QP_CA_MCH_API_H
#define _CODER_QP_CA_MCH_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_qp_ca_mch_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void qp_ca_mch(real_T arg[12], real_T B[27], real_T plim[18], real_T
                      rlim[18], real_T T, real_T Wv[9], real_T Wu[81], real_T
                      ud[9], real_T imax, real_T gam, boolean_T only_plim,
                      real_T u[9]);
extern void qp_ca_mch_api(const mxArray * const prhs[11], int32_T nlhs, const
  mxArray *plhs[1]);
extern void qp_ca_mch_atexit(void);
extern void qp_ca_mch_initialize(void);
extern void qp_ca_mch_terminate(void);
extern void qp_ca_mch_xil_terminate(void);

#endif

/*
 * File trailer for _coder_qp_ca_mch_api.h
 *
 * [EOF]
 */
