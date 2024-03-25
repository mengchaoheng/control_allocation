/*
 * Prerelease License - for engineering feedback and testing purposes
 * only. Not for sale.
 * File: _coder_allocator_dir_LPwrap_4_info.c
 *
 * MATLAB Coder version            : 24.1
 * C/C++ source code generated on  : 2024-03-24 20:52:32
 */

/* Include Files */
#include "_coder_allocator_dir_LPwrap_4_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : const mxArray *
 */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[5] = {
      "789c6360f4f465646060e06380801ddc109a17ca1780d24c0ca8005d9e114d1d23aa7206"
      "560616147d30f97e289d9c9f57925a5102e1e425e6a6c275a6e4e766"
      "e625e695845416a43214a516e7e794a5a68065d2327352433273538391397e205eae1b92"
      "149c039202b19d335293b3834b73198a328a112ecc41e6c0c363050e",
      "ffb210191e8e38c243004d3eda35563fb438b5a8583f3739433fa0283f4b1f142445f939"
      "f1893939f9c9892599f979fa3e01fa505e7e517c4a6651bc4f407951"
      "6241bc895e2e9abb0b2874373f0177c3e45d028a9313735253802e71768c7786db9f40a1"
      "fd6c38ed87c89466e695189a21fcbb8042fbac70da872a4f743ca105",
      "0c2c8208850f2c5c09b9179d46a8e700d38597fea7d2d3be8f46c704e9691f0c0c947d15"
      "38cc2336bd89e1b04f004dbed8d2c339b1c4d8c2d3b424dc3f3d2c23"
      "2dac38c9cd03e18e0002f6107207030e3ebdcc5f8043ff60cdb79496abb8c247008d2ece"
      "cc2dc8a9282d284a2d03d678f1ced42a575971da0f910195ab160cd4",
      "8b1f6b9cf6a1ca131d3fe801038d207a9503bd57e85baefe6bef50a3a77d3030dccbd5b0"
      "d0a2a0fc627f8fc2a00caf4ae7722f1fe76c7723e7a15fae0200131f"
      "baf6",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 3024U, &nameCaptureInfo);
  return nameCaptureInfo;
}

/*
 * Arguments    : void
 * Return Type  : mxArray *
 */
mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *propFieldName[9] = {"Version",
                                    "ResolvedFunctions",
                                    "Checksum",
                                    "EntryPoints",
                                    "CoverageInfo",
                                    "IsPolymorphic",
                                    "PropertyList",
                                    "UUID",
                                    "ClassEntryPointIsHandle"};
  const char_T *epFieldName[8] = {
      "Name",     "NumberOfInputs", "NumberOfOutputs", "ConstantInputs",
      "FullPath", "TimeStamp",      "Constructor",     "Visible"};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 1, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 4);
  emlrtSetField(xEntryPoints, 0, "Name",
                emlrtMxCreateString("allocator_dir_LPwrap_4"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(4.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "FullPath",
      emlrtMxCreateString(
          "/Users/mch/Proj/control_allocation/LP/allocator_dir_LPwrap_4.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739335.61800925923));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(
      xResult, 0, "Version",
      emlrtMxCreateString("24.1.0.2498408 (R2024a) Prerelease Update 3"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("x2VdqeFSYzeSNuhEbLpLGH"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/*
 * File trailer for _coder_allocator_dir_LPwrap_4_info.c
 *
 * [EOF]
 */
