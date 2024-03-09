/*
 * File: qp_ca_mch.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 17-Aug-2020 18:49:57
 */

#ifndef QP_CA_MCH_H
#define QP_CA_MCH_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "qp_ca_mch_types.h"

/* Function Declarations */
extern void qp_ca_mch(const double arg[12], const double B[27], const double
                      plim[18], const double rlim[18], double T, const double
                      Wv[9], const double Wu[81], const double ud[9], double
                      imax, double gam, boolean_T only_plim, double u[9]);
extern void qp_ca_mch_initialize(void);
extern void qp_ca_mch_terminate(void);

#endif

/*
 * File trailer for qp_ca_mch.h
 *
 * [EOF]
 */
