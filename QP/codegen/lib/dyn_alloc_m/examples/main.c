/*
 * File: main.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 18-Aug-2020 19:12:23
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
#include "dyn_alloc_m.h"
#include "main.h"

/* Function Declarations */
static void argInit_1x2_real_T(double result[2]);
static void argInit_3x1_real_T(double result[3]);
static void argInit_3x3_real_T(double result[9]);
static void argInit_3x9_real_T(double result[27]);
static void argInit_9x1_real_T(double result[9]);
static void argInit_9x3_real_T(double result[27]);
static void argInit_9x9_real_T(double result[81]);
static boolean_T argInit_boolean_T(void);
static double argInit_real_T(void);
static void main_dyn_alloc_m(void);

/* Function Definitions */

/*
 * Arguments    : double result[2]
 * Return Type  : void
 */
static void argInit_1x2_real_T(double result[2])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 2; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real_T();
  }
}

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
 * Arguments    : double result[9]
 * Return Type  : void
 */
static void argInit_3x3_real_T(double result[9])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    for (idx1 = 0; idx1 < 3; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 3 * idx1] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : double result[27]
 * Return Type  : void
 */
static void argInit_3x9_real_T(double result[27])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    for (idx1 = 0; idx1 < 9; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 3 * idx1] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : double result[9]
 * Return Type  : void
 */
static void argInit_9x1_real_T(double result[9])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 9; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[27]
 * Return Type  : void
 */
static void argInit_9x3_real_T(double result[27])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 9; idx0++) {
    for (idx1 = 0; idx1 < 3; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 9 * idx1] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : double result[81]
 * Return Type  : void
 */
static void argInit_9x9_real_T(double result[81])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 9; idx0++) {
    for (idx1 = 0; idx1 < 9; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 9 * idx1] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : boolean_T
 */
static boolean_T argInit_boolean_T(void)
{
  return false;
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
static void main_dyn_alloc_m(void)
{
  static double B[27];
  double v[3];
  static double u[9];
  double dv0[2];
  double dv1[2];
  static double dv2[9];
  static double dv3[81];
  static double dv4[81];
  static double dv5[27];

  /* Initialize function 'dyn_alloc_m' input arguments. */
  /* Initialize function input argument 'B'. */
  argInit_3x9_real_T(B);

  /* Initialize function input argument 'v'. */
  argInit_3x1_real_T(v);

  /* Initialize function input argument 'u'. */
  argInit_9x1_real_T(u);

  /* Initialize function input argument 'p_limits'. */
  /* Initialize function input argument 'v_limits'. */
  /* Initialize function input argument 'Wv'. */
  /* Initialize function input argument 'W1'. */
  /* Initialize function input argument 'W2'. */
  /* Initialize function input argument 'S'. */
  /* Call the entry-point 'dyn_alloc_m'. */
  argInit_1x2_real_T(dv0);
  argInit_1x2_real_T(dv1);
  argInit_3x3_real_T(dv2);
  argInit_9x9_real_T(dv3);
  argInit_9x9_real_T(dv4);
  argInit_9x3_real_T(dv5);
  dyn_alloc_m(B, v, u, dv0, dv1, argInit_real_T(), dv2, dv3, dv4, dv5,
              argInit_real_T(), argInit_real_T(), argInit_boolean_T());
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
  main_dyn_alloc_m();

  /* Terminate the application.
     You do not need to do this more than one time. */
  dyn_alloc_m_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
