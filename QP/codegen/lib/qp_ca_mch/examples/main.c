/*
 * File: main.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 17-Aug-2020 18:49:57
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
#include "qp_ca_mch.h"
#include "main.h"

/* Function Declarations */
static void argInit_12x1_real_T(double result[12]);
static void argInit_3x3_real_T(double result[9]);
static void argInit_3x9_real_T(double result[27]);
static void argInit_9x1_real_T(double result[9]);
static void argInit_9x2_real_T(double result[18]);
static void argInit_9x9_real_T(double result[81]);
static boolean_T argInit_boolean_T(void);
static double argInit_real_T(void);
static void main_qp_ca_mch(void);

/* Function Definitions */

/*
 * Arguments    : double result[12]
 * Return Type  : void
 */
static void argInit_12x1_real_T(double result[12])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 12; idx0++) {
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
 * Arguments    : double result[18]
 * Return Type  : void
 */
static void argInit_9x2_real_T(double result[18])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 9; idx0++) {
    for (idx1 = 0; idx1 < 2; idx1++) {
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
static void main_qp_ca_mch(void)
{
  double dv0[12];
  double dv1[27];
  double dv2[18];
  double dv3[18];
  double dv4[9];
  double dv5[81];
  double dv6[9];
  double u[9];

  /* Initialize function 'qp_ca_mch' input arguments. */
  /* Initialize function input argument 'arg'. */
  /* Initialize function input argument 'B'. */
  /* Initialize function input argument 'plim'. */
  /* Initialize function input argument 'rlim'. */
  /* Initialize function input argument 'Wv'. */
  /* Initialize function input argument 'Wu'. */
  /* Initialize function input argument 'ud'. */
  /* Call the entry-point 'qp_ca_mch'. */
  argInit_12x1_real_T(dv0);
  argInit_3x9_real_T(dv1);
  argInit_9x2_real_T(dv2);
  argInit_9x2_real_T(dv3);
  argInit_3x3_real_T(dv4);
  argInit_9x9_real_T(dv5);
  argInit_9x1_real_T(dv6);
  qp_ca_mch(dv0, dv1, dv2, dv3, argInit_real_T(), dv4, dv5, dv6, argInit_real_T(),
            argInit_real_T(), argInit_boolean_T(), u);
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
  main_qp_ca_mch();

  /* Terminate the application.
     You do not need to do this more than one time. */
  qp_ca_mch_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
