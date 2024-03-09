/*
 * File: wls_alloc_m.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 18-Aug-2020 19:20:00
 */

#ifndef WLS_ALLOC_M_H
#define WLS_ALLOC_M_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "wls_alloc_m_types.h"

/* Function Declarations */
extern void wls_alloc_m(const double B[27], const double v[3], double u[9],
  const double p_limits[2], const double v_limits[2], double T, const double Wv
  [9], const double Wu[81], const double ud[9], double imax, double gam,
  boolean_T only_plim);
extern void wls_alloc_m_initialize(void);
extern void wls_alloc_m_terminate(void);

#endif

/*
 * File trailer for wls_alloc_m.h
 *
 * [EOF]
 */
