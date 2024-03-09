/*
 * File: dir_alloc_mch.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Nov-2019 20:09:11
 */

#ifndef DIR_ALLOC_MCH_H
#define DIR_ALLOC_MCH_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "dir_alloc_mch_types.h"

/* Function Declarations */
extern void dir_alloc_mch(const double v[3], const double umin[4], const double
  umax[4], double u[4]);
extern void dir_alloc_mch_initialize(void);
extern void dir_alloc_mch_terminate(void);

#endif

/*
 * File trailer for dir_alloc_mch.h
 *
 * [EOF]
 */
