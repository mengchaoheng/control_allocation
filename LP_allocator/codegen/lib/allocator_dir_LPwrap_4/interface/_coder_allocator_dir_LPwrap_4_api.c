/*
 * Prerelease License - for engineering feedback and testing purposes
 * only. Not for sale.
 * File: _coder_allocator_dir_LPwrap_4_api.c
 *
 * MATLAB Coder version            : 24.1
 * C/C++ source code generated on  : 2024-03-24 20:52:32
 */

/* Include Files */
#include "_coder_allocator_dir_LPwrap_4_api.h"
#include "_coder_allocator_dir_LPwrap_4_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;

emlrtContext emlrtContextGlobal = {
    true,                                                 /* bFirstTime */
    false,                                                /* bInitialized */
    131643U,                                              /* fVersionInfo */
    NULL,                                                 /* fErrorFunction */
    "allocator_dir_LPwrap_4",                             /* fFunctionName */
    NULL,                                                 /* fRTCallStack */
    false,                                                /* bDebugMode */
    {2045744189U, 2170104910U, 2743257031U, 4284093946U}, /* fSigWrd */
    NULL                                                  /* fSigMem */
};

/* Function Declarations */
static real32_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                     const emlrtMsgIdentifier *parentId))[12];

static const mxArray *b_emlrt_marshallOut(const real32_T u);

static real32_T (*c_emlrt_marshallIn(const emlrtStack *sp,
                                     const mxArray *nullptr,
                                     const char_T *identifier))[3];

static const mxArray *c_emlrt_marshallOut(const uint16_T u);

static real32_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                     const emlrtMsgIdentifier *parentId))[3];

static real32_T (*e_emlrt_marshallIn(const emlrtStack *sp,
                                     const mxArray *nullptr,
                                     const char_T *identifier))[4];

static void emlrtExitTimeCleanupDtorFcn(const void *r);

static real32_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[12];

static const mxArray *emlrt_marshallOut(const real32_T u[4]);

static real32_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                     const emlrtMsgIdentifier *parentId))[4];

static real32_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                     const emlrtMsgIdentifier *msgId))[12];

static real32_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                     const emlrtMsgIdentifier *msgId))[3];

static real32_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                     const emlrtMsgIdentifier *msgId))[4];

/* Function Definitions */
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real32_T (*)[12]
 */
static real32_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                     const emlrtMsgIdentifier *parentId))[12]
{
  real32_T(*y)[12];
  y = g_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const real32_T u
 * Return Type  : const mxArray *
 */
static const mxArray *b_emlrt_marshallOut(const real32_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
  *(real32_T *)emlrtMxGetData(m) = u;
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : real32_T (*)[3]
 */
static real32_T (*c_emlrt_marshallIn(const emlrtStack *sp,
                                     const mxArray *nullptr,
                                     const char_T *identifier))[3]
{
  emlrtMsgIdentifier thisId;
  real32_T(*y)[3];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const uint16_T u
 * Return Type  : const mxArray *
 */
static const mxArray *c_emlrt_marshallOut(const uint16_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
  *(uint16_T *)emlrtMxGetData(m) = u;
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real32_T (*)[3]
 */
static real32_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                     const emlrtMsgIdentifier *parentId))[3]
{
  real32_T(*y)[3];
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : real32_T (*)[4]
 */
static real32_T (*e_emlrt_marshallIn(const emlrtStack *sp,
                                     const mxArray *nullptr,
                                     const char_T *identifier))[4]
{
  emlrtMsgIdentifier thisId;
  real32_T(*y)[4];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const void *r
 * Return Type  : void
 */
static void emlrtExitTimeCleanupDtorFcn(const void *r)
{
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : real32_T (*)[12]
 */
static real32_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[12]
{
  emlrtMsgIdentifier thisId;
  real32_T(*y)[12];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const real32_T u[4]
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real32_T u[4])
{
  static const int32_T i = 0;
  static const int32_T i1 = 4;
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxSINGLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, &i1, 1);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real32_T (*)[4]
 */
static real32_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                     const emlrtMsgIdentifier *parentId))[4]
{
  real32_T(*y)[4];
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real32_T (*)[12]
 */
static real32_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                     const emlrtMsgIdentifier *msgId))[12]
{
  static const int32_T dims[2] = {3, 4};
  int32_T iv[2];
  real32_T(*ret)[12];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "single", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real32_T(*)[12])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real32_T (*)[3]
 */
static real32_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                     const emlrtMsgIdentifier *msgId))[3]
{
  static const int32_T dims = 3;
  int32_T i;
  real32_T(*ret)[3];
  boolean_T b = false;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "single", false, 1U,
                            (const void *)&dims, &b, &i);
  ret = (real32_T(*)[3])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real32_T (*)[4]
 */
static real32_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                     const emlrtMsgIdentifier *msgId))[4]
{
  static const int32_T dims = 4;
  int32_T i;
  real32_T(*ret)[4];
  boolean_T b = false;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "single", false, 1U,
                            (const void *)&dims, &b, &i);
  ret = (real32_T(*)[4])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const mxArray * const prhs[4]
 *                int32_T nlhs
 *                const mxArray *plhs[3]
 * Return Type  : void
 */
void allocator_dir_LPwrap_4_api(const mxArray *const prhs[4], int32_T nlhs,
                                const mxArray *plhs[3])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  real32_T(*B)[12];
  real32_T(*u)[4];
  real32_T(*umax)[4];
  real32_T(*umin)[4];
  real32_T(*v)[3];
  real32_T z;
  uint16_T iters;
  st.tls = emlrtRootTLSGlobal;
  u = (real32_T(*)[4])mxMalloc(sizeof(real32_T[4]));
  /* Marshall function inputs */
  B = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "B");
  v = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "v");
  umin = e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "umin");
  umax = e_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "umax");
  /* Invoke the target function */
  allocator_dir_LPwrap_4(*B, *v, *umin, *umax, *u, &z, &iters);
  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*u);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(z);
  }
  if (nlhs > 2) {
    plhs[2] = c_emlrt_marshallOut(iters);
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void allocator_dir_LPwrap_4_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtPushHeapReferenceStackR2021a(
      &st, false, NULL, (void *)&emlrtExitTimeCleanupDtorFcn, NULL, NULL, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  allocator_dir_LPwrap_4_xil_terminate();
  allocator_dir_LPwrap_4_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void allocator_dir_LPwrap_4_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void allocator_dir_LPwrap_4_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_allocator_dir_LPwrap_4_api.c
 *
 * [EOF]
 */
