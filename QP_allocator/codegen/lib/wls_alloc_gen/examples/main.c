/*
 * File: main.c
 *
 * MATLAB Coder version            : 24.1
 * C/C++ source code generated on  : 2025-06-16 11:03:22
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include Files */
#include "main.h"
#include "wls_alloc_gen.h"

/* Function Declarations */
static void argInit_3x1_real32_T(float result[3]);

static void argInit_3x3_real32_T(float result[9]);

static void argInit_3x4_real32_T(float result[12]);

static void argInit_4x1_real32_T(float result[4]);

static void argInit_4x4_real32_T(float result[16]);

static float argInit_real32_T(void);

/* Function Definitions */
/*
 * Arguments    : float result[3]
 * Return Type  : void
 */
static void argInit_3x1_real32_T(float result[3])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real32_T();
  }
}

/*
 * Arguments    : float result[9]
 * Return Type  : void
 */
static void argInit_3x3_real32_T(float result[9])
{
  int i;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 9; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[i] = argInit_real32_T();
  }
}

/*
 * Arguments    : float result[12]
 * Return Type  : void
 */
static void argInit_3x4_real32_T(float result[12])
{
  int i;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 12; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[i] = argInit_real32_T();
  }
}

/*
 * Arguments    : float result[4]
 * Return Type  : void
 */
static void argInit_4x1_real32_T(float result[4])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 4; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real32_T();
  }
}

/*
 * Arguments    : float result[16]
 * Return Type  : void
 */
static void argInit_4x4_real32_T(float result[16])
{
  int i;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 16; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[i] = argInit_real32_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : float
 */
static float argInit_real32_T(void)
{
  return 0.0F;
}

/*
 * Arguments    : int argc
 *                char **argv
 * Return Type  : int
 */
int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;
  /* The initialize function is being called automatically from your entry-point
   * function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
You can call entry-point functions multiple times. */
  main_wls_alloc_gen();
  /* Terminate the application.
You do not need to do this more than one time. */
  wls_alloc_gen_terminate();
  return 0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void main_wls_alloc_gen(void)
{
  static float fv3[16];
  static float fv[12];
  float fv2[9];
  float b_u[4];
  float c_u[4];
  float u[4];
  float fv1[3];
  float gam_tmp;
  /* Initialize function 'wls_alloc_gen' input arguments. */
  /* Initialize function input argument 'B'. */
  /* Initialize function input argument 'v'. */
  /* Initialize function input argument 'umin'. */
  argInit_4x1_real32_T(u);
  /* Initialize function input argument 'umax'. */
  /* Initialize function input argument 'Wv'. */
  /* Initialize function input argument 'Wu'. */
  /* Initialize function input argument 'ud'. */
  gam_tmp = argInit_real32_T();
  /* Initialize function input argument 'u'. */
  /* Initialize function input argument 'W'. */
  /* Call the entry-point 'wls_alloc_gen'. */
  argInit_3x4_real32_T(fv);
  argInit_3x1_real32_T(fv1);
  argInit_3x3_real32_T(fv2);
  argInit_4x4_real32_T(fv3);
  b_u[0] = u[0];
  b_u[1] = u[1];
  b_u[2] = u[2];
  b_u[3] = u[3];
  c_u[0] = u[0];
  c_u[1] = u[1];
  c_u[2] = u[2];
  c_u[3] = u[3];
  wls_alloc_gen(fv, fv1, b_u, b_u, fv2, fv3, b_u, gam_tmp, u, c_u, gam_tmp,
                gam_tmp);
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
