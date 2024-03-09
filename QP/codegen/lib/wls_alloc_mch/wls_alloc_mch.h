/*
 * File: wls_alloc_mch.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 08-Aug-2020 14:31:14
 */

#ifndef WLS_ALLOC_MCH_H
#define WLS_ALLOC_MCH_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "wls_alloc_mch_types.h"
//----------------------------by mch
#include "userlib.h"
#include "motor.h"
typedef struct
{
	double u[8];//¶æÁ¿£¬µ¥Î»£º»¡¶È
	double p_limits;//¶æ»ú·ùÖµ£¬µ¥Î»£º¶È
}sWls;
//--------------------------------
/* Function Declarations */
extern void wls_alloc_mch(const double v[3], double u[8]);
extern void wls_alloc_mch_initialize(void);
extern void wls_alloc_mch_terminate(void);
//----------------------------by mch
extern sWls wls;
//--------------------------------
#endif

/*
 * File trailer for wls_alloc_mch.h
 *
 * [EOF]
 */
