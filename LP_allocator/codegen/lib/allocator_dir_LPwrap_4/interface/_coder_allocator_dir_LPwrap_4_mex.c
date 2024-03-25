/*
 * Prerelease License - for engineering feedback and testing purposes
 * only. Not for sale.
 * File: _coder_allocator_dir_LPwrap_4_mex.c
 *
 * MATLAB Coder version            : 24.1
 * C/C++ source code generated on  : 2024-03-24 20:52:32
 */

/* Include Files */
#include "_coder_allocator_dir_LPwrap_4_mex.h"
#include "_coder_allocator_dir_LPwrap_4_api.h"

/* Function Definitions */
/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[]
 *                int32_T nrhs
 *                const mxArray *prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&allocator_dir_LPwrap_4_atexit);
  /* Module initialization. */
  allocator_dir_LPwrap_4_initialize();
  /* Dispatch the entry-point. */
  unsafe_allocator_dir_LPwrap_4_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  allocator_dir_LPwrap_4_terminate();
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[3]
 *                int32_T nrhs
 *                const mxArray *prhs[4]
 * Return Type  : void
 */
void unsafe_allocator_dir_LPwrap_4_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                               int32_T nrhs,
                                               const mxArray *prhs[4])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *b_prhs[4];
  const mxArray *outputs[3];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 4, 4,
                        22, "allocator_dir_LPwrap_4");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 22,
                        "allocator_dir_LPwrap_4");
  }
  /* Call the function. */
  b_prhs[0] = prhs[0];
  b_prhs[1] = prhs[1];
  b_prhs[2] = prhs[2];
  b_prhs[3] = prhs[3];
  allocator_dir_LPwrap_4_api(b_prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

/*
 * File trailer for _coder_allocator_dir_LPwrap_4_mex.c
 *
 * [EOF]
 */
