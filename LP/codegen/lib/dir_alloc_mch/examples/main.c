/*
 * File: main.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Nov-2019 20:09:11
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
#include "dir_alloc_mch.h"
#include "main.h"

/* Function Declarations */
static void argInit_3x1_real_T(double result[3]);
static void argInit_4x1_real_T(double result[4]);
static double argInit_real_T(void);
static void main_dir_alloc_mch(void);

/* Function Definitions */

/*
 * Arguments    : double result[3]
 * Return Type  : void
 */
static void argInit_3x1_real_T(double result[3])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[4]
 * Return Type  : void
 */
static void argInit_4x1_real_T(double result[4])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 4; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_dir_alloc_mch(void)
{
  double dv4[3];
  double dv5[4];
  double dv6[4];
  double u[4];

  /* Initialize function 'dir_alloc_mch' input arguments. */
  /* Initialize function input argument 'v'. */
  /* Initialize function input argument 'umin'. */
  /* Initialize function input argument 'umax'. */
  /* Call the entry-point 'dir_alloc_mch'. */
  argInit_3x1_real_T(dv4);
  argInit_4x1_real_T(dv5);
  argInit_4x1_real_T(dv6);
  dir_alloc_mch(dv4, dv5, dv6, u);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_dir_alloc_mch();

  /* Terminate the application.
     You do not need to do this more than one time. */
  dir_alloc_mch_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
