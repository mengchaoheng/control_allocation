/*
 * File: qp_ca_mch.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 17-Aug-2020 18:49:57
 */

/* Include Files */
#include <string.h>
#include <math.h>
#include "qp_ca_mch.h"

/* Function Declarations */
static void LSQFromQR(const double A_data[], const int A_size[2], const double
                      tau_data[], const int jpvt_data[], double B[12], int rankA,
                      double Y_data[], int Y_size[1]);
static void qrsolve(const double A_data[], const int A_size[2], const double B
                    [12], double Y_data[], int Y_size[1]);
static int rankFromQR(const double A_data[], const int A_size[2]);
static void rdivide(const double x_data[], const int x_size[1], const double
                    y_data[], double z_data[], int z_size[1]);
static double rt_hypotd(double u0, double u1);
static void wls_alloc_mch(const double B[27], const double v[3], const double
  umin[9], const double umax[9], const double Wv[9], const double Wu[81], const
  double ud[9], double gam, double u[9], double W[9], double imax);
static double xnrm2(int n, const double x_data[], int ix0);
static void xzlarf(int m, int n, int iv0, double tau, double C_data[], int ic0,
                   double work_data[]);

/* Function Definitions */

/*
 * Arguments    : const double A_data[]
 *                const int A_size[2]
 *                const double tau_data[]
 *                const int jpvt_data[]
 *                double B[12]
 *                int rankA
 *                double Y_data[]
 *                int Y_size[1]
 * Return Type  : void
 */
static void LSQFromQR(const double A_data[], const int A_size[2], const double
                      tau_data[], const int jpvt_data[], double B[12], int rankA,
                      double Y_data[], int Y_size[1])
{
  int loop_ub;
  int i;
  double wj;
  Y_size[0] = (signed char)A_size[1];
  loop_ub = (signed char)A_size[1];
  if (0 <= loop_ub - 1) {
    memset(&Y_data[0], 0, (unsigned int)(loop_ub * (int)sizeof(double)));
  }

  for (loop_ub = 0; loop_ub < A_size[1]; loop_ub++) {
    if (tau_data[loop_ub] != 0.0) {
      wj = B[loop_ub];
      for (i = loop_ub + 1; i + 1 < 13; i++) {
        wj += A_data[i + A_size[0] * loop_ub] * B[i];
      }

      wj *= tau_data[loop_ub];
      if (wj != 0.0) {
        B[loop_ub] -= wj;
        for (i = loop_ub + 1; i + 1 < 13; i++) {
          B[i] -= A_data[i + A_size[0] * loop_ub] * wj;
        }
      }
    }
  }

  for (i = 0; i < rankA; i++) {
    Y_data[jpvt_data[i] - 1] = B[i];
  }

  for (loop_ub = rankA - 1; loop_ub + 1 > 0; loop_ub--) {
    Y_data[jpvt_data[loop_ub] - 1] /= A_data[loop_ub + A_size[0] * loop_ub];
    for (i = 0; i < loop_ub; i++) {
      Y_data[jpvt_data[i] - 1] -= Y_data[jpvt_data[loop_ub] - 1] * A_data[i +
        A_size[0] * loop_ub];
    }
  }
}

/*
 * Arguments    : const double A_data[]
 *                const int A_size[2]
 *                const double B[12]
 *                double Y_data[]
 *                int Y_size[1]
 * Return Type  : void
 */
static void qrsolve(const double A_data[], const int A_size[2], const double B
                    [12], double Y_data[], int Y_size[1])
{
  int b_A_size[2];
  int n;
  double b_A_data[108];
  int b_n;
  int jpvt_data[9];
  int yk;
  int k;
  double b_B[12];
  double tau_data[9];
  double work_data[9];
  int i;
  double smax;
  double s;
  int i_i;
  int nmi;
  double absxk;
  double vn1_data[9];
  double vn2_data[9];
  double t;
  int pvt;
  b_A_size[0] = 12;
  b_A_size[1] = A_size[1];
  n = A_size[0] * A_size[1];
  if (0 <= n - 1) {
    memcpy(&b_A_data[0], &A_data[0], (unsigned int)(n * (int)sizeof(double)));
  }

  b_n = A_size[1];
  if (A_size[1] < 1) {
    n = 0;
  } else {
    n = A_size[1];
  }

  if (n > 0) {
    jpvt_data[0] = 1;
    yk = 1;
    for (k = 2; k <= n; k++) {
      yk++;
      jpvt_data[k - 1] = yk;
    }
  }

  if (A_size[1] != 0) {
    n = (signed char)A_size[1];
    if (0 <= n - 1) {
      memset(&work_data[0], 0, (unsigned int)(n * (int)sizeof(double)));
    }

    k = 1;
    for (yk = 0; yk < b_n; yk++) {
      smax = 0.0;
      s = 3.3121686421112381E-170;
      for (n = k; n <= k + 11; n++) {
        absxk = fabs(A_data[n - 1]);
        if (absxk > s) {
          t = s / absxk;
          smax = 1.0 + smax * t * t;
          s = absxk;
        } else {
          t = absxk / s;
          smax += t * t;
        }
      }

      smax = s * sqrt(smax);
      vn1_data[yk] = smax;
      vn2_data[yk] = vn1_data[yk];
      k += 12;
    }

    for (i = 0; i < b_n; i++) {
      i_i = i + i * 12;
      nmi = b_n - i;
      if (nmi < 1) {
        n = 0;
      } else {
        n = 1;
        if (nmi > 1) {
          yk = i;
          smax = fabs(vn1_data[i]);
          for (k = 2; k <= nmi; k++) {
            yk++;
            s = fabs(vn1_data[yk]);
            if (s > smax) {
              n = k;
              smax = s;
            }
          }
        }
      }

      pvt = (i + n) - 1;
      if (pvt + 1 != i + 1) {
        yk = 12 * pvt;
        n = 12 * i;
        for (k = 0; k < 12; k++) {
          smax = b_A_data[yk];
          b_A_data[yk] = b_A_data[n];
          b_A_data[n] = smax;
          yk++;
          n++;
        }

        n = jpvt_data[pvt];
        jpvt_data[pvt] = jpvt_data[i];
        jpvt_data[i] = n;
        vn1_data[pvt] = vn1_data[i];
        vn2_data[pvt] = vn2_data[i];
      }

      s = b_A_data[i_i];
      absxk = 0.0;
      smax = xnrm2(11 - i, b_A_data, i_i + 2);
      if (smax != 0.0) {
        smax = rt_hypotd(b_A_data[i_i], smax);
        if (b_A_data[i_i] >= 0.0) {
          smax = -smax;
        }

        if (fabs(smax) < 1.0020841800044864E-292) {
          yk = 0;
          n = (i_i - i) + 12;
          do {
            yk++;
            for (k = i_i + 1; k < n; k++) {
              b_A_data[k] *= 9.9792015476736E+291;
            }

            smax *= 9.9792015476736E+291;
            s *= 9.9792015476736E+291;
          } while (!(fabs(smax) >= 1.0020841800044864E-292));

          smax = rt_hypotd(s, xnrm2(11 - i, b_A_data, i_i + 2));
          if (s >= 0.0) {
            smax = -smax;
          }

          absxk = (smax - s) / smax;
          s = 1.0 / (s - smax);
          n = (i_i - i) + 12;
          for (k = i_i + 1; k < n; k++) {
            b_A_data[k] *= s;
          }

          for (k = 1; k <= yk; k++) {
            smax *= 1.0020841800044864E-292;
          }

          s = smax;
        } else {
          absxk = (smax - b_A_data[i_i]) / smax;
          s = 1.0 / (b_A_data[i_i] - smax);
          n = (i_i - i) + 12;
          for (k = i_i + 1; k < n; k++) {
            b_A_data[k] *= s;
          }

          s = smax;
        }
      }

      tau_data[i] = absxk;
      b_A_data[i_i] = s;
      if (i + 1 < b_n) {
        s = b_A_data[i_i];
        b_A_data[i_i] = 1.0;
        xzlarf(12 - i, nmi - 1, i_i + 1, tau_data[i], b_A_data, (i + (i + 1) *
                12) + 1, work_data);
        b_A_data[i_i] = s;
      }

      for (yk = i + 1; yk < b_n; yk++) {
        if (vn1_data[yk] != 0.0) {
          smax = fabs(b_A_data[i + 12 * yk]) / vn1_data[yk];
          smax = 1.0 - smax * smax;
          if (smax < 0.0) {
            smax = 0.0;
          }

          s = vn1_data[yk] / vn2_data[yk];
          s = smax * (s * s);
          if (s <= 1.4901161193847656E-8) {
            vn1_data[yk] = xnrm2(11 - i, b_A_data, (i + 12 * yk) + 2);
            vn2_data[yk] = vn1_data[yk];
          } else {
            vn1_data[yk] *= sqrt(smax);
          }
        }
      }
    }
  }

  memcpy(&b_B[0], &B[0], 12U * sizeof(double));
  LSQFromQR(b_A_data, b_A_size, tau_data, jpvt_data, b_B, rankFromQR(b_A_data,
             b_A_size), Y_data, Y_size);
}

/*
 * Arguments    : const double A_data[]
 *                const int A_size[2]
 * Return Type  : int
 */
static int rankFromQR(const double A_data[], const int A_size[2])
{
  int r;
  double tol;
  r = 0;
  if (A_size[1] > 0) {
    tol = 12.0 * fabs(A_data[0]) * 2.2204460492503131E-16;
    while ((r < A_size[1]) && (!(fabs(A_data[r + A_size[0] * r]) <= tol))) {
      r++;
    }
  }

  return r;
}

/*
 * Arguments    : const double x_data[]
 *                const int x_size[1]
 *                const double y_data[]
 *                double z_data[]
 *                int z_size[1]
 * Return Type  : void
 */
static void rdivide(const double x_data[], const int x_size[1], const double
                    y_data[], double z_data[], int z_size[1])
{
  int loop_ub;
  int i0;
  z_size[0] = x_size[0];
  loop_ub = x_size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    z_data[i0] = x_data[i0] / y_data[i0];
  }
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

/*
 * WLS_ALLOC - Control allocation using weighted least squares.
 *
 *   [u,W,iter] = wls_alloc(B,v,umin,umax,[Wv,Wu,ud,gamma,u0,W0,imax],m)
 *
 *  Solves the weighted, bounded least-squares problem
 *
 *    min ||Wu(u-ud)||^2 + gamma ||Wv(Bu-v)||^2
 *
 *    subj. to  umin <= u <= umax
 *
 *  using an active set method.
 *
 *   Inputs:
 *   -------
 *  B     control effectiveness matrix (k x m)
 *  v     commanded virtual control (k x 1)
 *  umin  lower position limits (m x 1)
 *  umax  upper position limits (m x 1)
 *  Wv    virtual control weighting matrix (k x k) [I]
 *  Wu    control weighting matrix (m x m) [I]
 *  ud    desired control (m x 1) [0]
 *  gamma weight (scalar) [1e6]
 *  u0    initial point (m x 1)
 *  W0    initial working set (m x 1) [empty]
 *  imax  max no. of iterations [100]
 *
 *   Outputs:
 *   -------
 *  u     optimal control
 *  W     optimal active set
 *  iter  no. of iterations (= no. of changes in the working set + 1)
 *
 *                             0 if u_i not saturated
 *  Working set syntax: W_i = -1 if u_i = umin_i
 *                            +1 if u_i = umax_i
 *
 *  See also: WLSC_ALLOC, IP_ALLOC, FXP_ALLOC, QP_SIM.
 * Arguments    : const double B[27]
 *                const double v[3]
 *                const double umin[9]
 *                const double umax[9]
 *                const double Wv[9]
 *                const double Wu[81]
 *                const double ud[9]
 *                double gam
 *                double u[9]
 *                double W[9]
 *                double imax
 * Return Type  : void
 */
static void wls_alloc_mch(const double B[27], const double v[3], const double
  umin[9], const double umax[9], const double Wv[9], const double Wu[81], const
  double ud[9], double gam, double u[9], double W[9], double imax)
{
  double gam_sq;
  int i1;
  int aoffset;
  double b_gam_sq[27];
  int partialTrueCount;
  double A[108];
  double c_gam_sq[3];
  double b_Wu[9];
  double d_gam_sq[12];
  double b_A[12];
  int i;
  boolean_T i_free[9];
  double d[12];
  int iter;
  boolean_T exitg1;
  int trueCount;
  int A_free_size[2];
  int tmp_data[9];
  double A_free_data[108];
  double p_free_data[9];
  int p_free_size[1];
  double u_opt[9];
  int b_tmp_data[9];
  boolean_T y;
  boolean_T x_data[9];
  boolean_T b_u_opt[9];
  double p[9];
  boolean_T exitg2;
  boolean_T bv0[9];
  int umin_size[1];
  int c_tmp_data[9];
  double p_data[9];
  double d_tmp_data[9];
  int tmp_size[1];
  int umax_size[1];
  int e_tmp_data[9];
  int idx;

  /*  Number of variables   */
  gam_sq = sqrt(gam);
  for (i1 = 0; i1 < 3; i1++) {
    for (aoffset = 0; aoffset < 9; aoffset++) {
      b_gam_sq[i1 + 3 * aoffset] = 0.0;
      for (partialTrueCount = 0; partialTrueCount < 3; partialTrueCount++) {
        b_gam_sq[i1 + 3 * aoffset] += gam_sq * Wv[i1 + 3 * partialTrueCount] *
          B[partialTrueCount + 3 * aoffset];
      }
    }
  }

  for (i1 = 0; i1 < 9; i1++) {
    for (aoffset = 0; aoffset < 3; aoffset++) {
      A[aoffset + 12 * i1] = b_gam_sq[aoffset + 3 * i1];
    }

    memcpy(&A[i1 * 12 + 3], &Wu[i1 * 9], 9U * sizeof(double));
  }

  /*  Initial residual. */
  for (i1 = 0; i1 < 3; i1++) {
    c_gam_sq[i1] = 0.0;
    for (aoffset = 0; aoffset < 3; aoffset++) {
      c_gam_sq[i1] += gam_sq * Wv[i1 + 3 * aoffset] * v[aoffset];
    }
  }

  for (i1 = 0; i1 < 9; i1++) {
    b_Wu[i1] = 0.0;
    for (aoffset = 0; aoffset < 9; aoffset++) {
      b_Wu[i1] += Wu[i1 + 9 * aoffset] * ud[aoffset];
    }
  }

  for (i1 = 0; i1 < 3; i1++) {
    d_gam_sq[i1] = c_gam_sq[i1];
  }

  memcpy(&d_gam_sq[3], &b_Wu[0], 9U * sizeof(double));
  for (i1 = 0; i1 < 12; i1++) {
    b_A[i1] = 0.0;
    for (aoffset = 0; aoffset < 9; aoffset++) {
      b_A[i1] += A[i1 + 12 * aoffset] * u[aoffset];
    }

    d[i1] = d_gam_sq[i1] - b_A[i1];
  }

  /*  Determine indeces of free variables. */
  for (i = 0; i < 9; i++) {
    i_free[i] = (W[i] == 0.0);
  }

  /*  Iterate until optimum is found or maximum number of iterations */
  /*  is reached. */
  iter = 0;
  exitg1 = false;
  while ((!exitg1) && (iter <= (int)imax - 1)) {
    /*  ---------------------------------------- */
    /*   Compute optimal perturbation vector p. */
    /*  ---------------------------------------- */
    /*  Eliminate saturated variables. */
    trueCount = 0;
    for (i = 0; i < 9; i++) {
      if (i_free[i]) {
        trueCount++;
      }
    }

    partialTrueCount = 0;
    for (i = 0; i < 9; i++) {
      if (i_free[i]) {
        tmp_data[partialTrueCount] = i + 1;
        partialTrueCount++;
      }
    }

    A_free_size[0] = 12;
    A_free_size[1] = trueCount;
    for (i1 = 0; i1 < trueCount; i1++) {
      for (aoffset = 0; aoffset < 12; aoffset++) {
        A_free_data[aoffset + 12 * i1] = A[aoffset + 12 * (tmp_data[i1] - 1)];
      }
    }

    /*  Solve the reduced optimization problem for free variables. */
    if (trueCount == 0) {
      p_free_size[0] = 0;
    } else {
      qrsolve(A_free_data, A_free_size, d, p_free_data, p_free_size);
    }

    /*  Zero all perturbations corresponding to active constraints. */
    /*  Insert perturbations from p_free into free the variables. */
    partialTrueCount = 0;

    /*  ---------------------------- */
    /*   Is the new point feasible? */
    /*  ---------------------------- */
    aoffset = 0;
    for (i = 0; i < 9; i++) {
      gam_sq = 0.0;
      if (i_free[i]) {
        gam_sq = p_free_data[partialTrueCount];
        partialTrueCount++;
      }

      u_opt[i] = u[i] + gam_sq;
      if (i_free[i]) {
        aoffset++;
      }

      p[i] = gam_sq;
    }

    partialTrueCount = 0;
    for (i = 0; i < 9; i++) {
      if (i_free[i]) {
        b_tmp_data[partialTrueCount] = i + 1;
        partialTrueCount++;
      }

      b_u_opt[i] = ((u_opt[i] < umin[i]) || (u_opt[i] > umax[i]));
    }

    for (i1 = 0; i1 < aoffset; i1++) {
      x_data[i1] = b_u_opt[b_tmp_data[i1] - 1];
    }

    y = false;
    partialTrueCount = 1;
    exitg2 = false;
    while ((!exitg2) && (partialTrueCount <= aoffset)) {
      if (x_data[partialTrueCount - 1]) {
        y = true;
        exitg2 = true;
      } else {
        partialTrueCount++;
      }
    }

    if (!y) {
      /*  ---------------------------- */
      /*   Yes, check for optimality. */
      /*  ---------------------------- */
      /*  Update point and residual. */
      memcpy(&u[0], &u_opt[0], 9U * sizeof(double));
      if ((trueCount == 1) || (p_free_size[0] == 1)) {
        for (i1 = 0; i1 < 12; i1++) {
          d_gam_sq[i1] = 0.0;
          for (aoffset = 0; aoffset < trueCount; aoffset++) {
            d_gam_sq[i1] += A_free_data[i1 + 12 * aoffset] * p_free_data[aoffset];
          }
        }
      } else {
        memset(&d_gam_sq[0], 0, 12U * sizeof(double));
        for (partialTrueCount = 0; partialTrueCount < trueCount;
             partialTrueCount++) {
          if (p_free_data[partialTrueCount] != 0.0) {
            aoffset = partialTrueCount * 12;
            for (i = 0; i < 12; i++) {
              i1 = aoffset + i;
              d_gam_sq[i] += p_free_data[partialTrueCount] * A[i1 % 12 + 12 *
                (tmp_data[i1 / 12] - 1)];
            }
          }
        }
      }

      for (i1 = 0; i1 < 12; i1++) {
        d[i1] -= d_gam_sq[i1];
      }

      /*  Compute Lagrangian multipliers. */
      for (i1 = 0; i1 < 9; i1++) {
        b_Wu[i1] = 0.0;
        for (aoffset = 0; aoffset < 12; aoffset++) {
          b_Wu[i1] += A[aoffset + 12 * i1] * d[aoffset];
        }

        u_opt[i1] = W[i1] * b_Wu[i1];
      }

      /*  Are all lambda non-negative? */
      y = true;
      partialTrueCount = 1;
      exitg2 = false;
      while ((!exitg2) && (partialTrueCount < 10)) {
        if (!(u_opt[partialTrueCount - 1] >= -2.2204460492503131E-16)) {
          y = false;
          exitg2 = true;
        } else {
          partialTrueCount++;
        }
      }

      if (y) {
        /*  / ------------------------ \ */
        /*  | Optimum found, bail out. | */
        /*  \ ------------------------ / */
        exitg1 = true;
      } else {
        /*  -------------------------------------------------- */
        /*   Optimum not found, remove one active constraint. */
        /*  -------------------------------------------------- */
        /*  Remove constraint with most negative lambda from the */
        /*  working set. */
        gam_sq = u_opt[0];
        idx = 0;
        for (partialTrueCount = 0; partialTrueCount < 8; partialTrueCount++) {
          if (gam_sq > u_opt[partialTrueCount + 1]) {
            gam_sq = u_opt[partialTrueCount + 1];
            idx = partialTrueCount + 1;
          }
        }

        W[idx] = 0.0;
        i_free[idx] = true;
        iter++;
      }
    } else {
      /*  --------------------------------------- */
      /*   No, find primary bounding constraint. */
      /*  --------------------------------------- */
      /*  Compute distances to the different boundaries. Since alpha < 1 */
      /*  is the maximum step length, initiate with ones. */
      aoffset = 0;
      for (i = 0; i < 9; i++) {
        u_opt[i] = 1.0;
        y = (p[i] < 0.0);
        bv0[i] = (p[i] > 0.0);
        if (i_free[i] && y) {
          aoffset++;
        }

        b_u_opt[i] = y;
      }

      partialTrueCount = 0;
      for (i = 0; i < 9; i++) {
        if (i_free[i] && b_u_opt[i]) {
          c_tmp_data[partialTrueCount] = i + 1;
          partialTrueCount++;
        }
      }

      umin_size[0] = aoffset;
      for (i1 = 0; i1 < aoffset; i1++) {
        b_Wu[i1] = umin[c_tmp_data[i1] - 1] - u[c_tmp_data[i1] - 1];
      }

      for (i1 = 0; i1 < aoffset; i1++) {
        p_data[i1] = p[c_tmp_data[i1] - 1];
      }

      rdivide(b_Wu, umin_size, p_data, d_tmp_data, tmp_size);
      partialTrueCount = 0;
      aoffset = 0;
      for (i = 0; i < 9; i++) {
        if (i_free[i] && b_u_opt[i]) {
          u_opt[i] = d_tmp_data[partialTrueCount];
          partialTrueCount++;
        }

        if (i_free[i] && bv0[i]) {
          aoffset++;
        }
      }

      partialTrueCount = 0;
      for (i = 0; i < 9; i++) {
        if (i_free[i] && bv0[i]) {
          e_tmp_data[partialTrueCount] = i + 1;
          partialTrueCount++;
        }
      }

      umax_size[0] = aoffset;
      for (i1 = 0; i1 < aoffset; i1++) {
        b_Wu[i1] = umax[e_tmp_data[i1] - 1] - u[e_tmp_data[i1] - 1];
      }

      for (i1 = 0; i1 < aoffset; i1++) {
        p_data[i1] = p[e_tmp_data[i1] - 1];
      }

      rdivide(b_Wu, umax_size, p_data, d_tmp_data, tmp_size);
      partialTrueCount = 0;
      for (i = 0; i < 9; i++) {
        if (i_free[i] && bv0[i]) {
          u_opt[i] = d_tmp_data[partialTrueCount];
          partialTrueCount++;
        }
      }

      /*  Proportion of p to travel */
      gam_sq = u_opt[0];
      idx = 0;
      for (partialTrueCount = 0; partialTrueCount < 8; partialTrueCount++) {
        if (gam_sq > u_opt[partialTrueCount + 1]) {
          gam_sq = u_opt[partialTrueCount + 1];
          idx = partialTrueCount + 1;
        }
      }

      /*  Update point and residual. */
      for (i1 = 0; i1 < 9; i1++) {
        u[i1] += gam_sq * p[i1];
      }

      partialTrueCount = 12 * trueCount - 1;
      for (i1 = 0; i1 <= partialTrueCount; i1++) {
        A_free_data[i1] *= gam_sq;
      }

      if ((trueCount == 1) || (p_free_size[0] == 1)) {
        for (i1 = 0; i1 < 12; i1++) {
          d_gam_sq[i1] = 0.0;
          for (aoffset = 0; aoffset < trueCount; aoffset++) {
            d_gam_sq[i1] += A_free_data[i1 + 12 * aoffset] * p_free_data[aoffset];
          }
        }
      } else {
        memset(&d_gam_sq[0], 0, 12U * sizeof(double));
        for (partialTrueCount = 0; partialTrueCount < trueCount;
             partialTrueCount++) {
          if (p_free_data[partialTrueCount] != 0.0) {
            aoffset = partialTrueCount * 12;
            for (i = 0; i < 12; i++) {
              d_gam_sq[i] += p_free_data[partialTrueCount] * A_free_data[aoffset
                + i];
            }
          }
        }
      }

      for (i1 = 0; i1 < 12; i1++) {
        d[i1] -= d_gam_sq[i1];
      }

      /*  Add corresponding constraint to working set. */
      gam_sq = p[idx];
      if (p[idx] < 0.0) {
        gam_sq = -1.0;
      } else {
        if (p[idx] > 0.0) {
          gam_sq = 1.0;
        }
      }

      W[idx] = gam_sq;
      i_free[idx] = false;
      iter++;
    }
  }
}

/*
 * Arguments    : int n
 *                const double x_data[]
 *                int ix0
 * Return Type  : double
 */
static double xnrm2(int n, const double x_data[], int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  kend = (ix0 + n) - 1;
  for (k = ix0; k <= kend; k++) {
    absxk = fabs(x_data[k - 1]);
    if (absxk > scale) {
      t = scale / absxk;
      y = 1.0 + y * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * sqrt(y);
}

/*
 * Arguments    : int m
 *                int n
 *                int iv0
 *                double tau
 *                double C_data[]
 *                int ic0
 *                double work_data[]
 * Return Type  : void
 */
static void xzlarf(int m, int n, int iv0, double tau, double C_data[], int ic0,
                   double work_data[])
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  int jy;
  int j;
  int i2;
  int ia;
  int exitg1;
  double c;
  int ix;
  if (tau != 0.0) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C_data[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      i = ic0 + (lastc - 1) * 12;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C_data[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      for (i = 1; i <= lastc; i++) {
        work_data[i - 1] = 0.0;
      }

      i = 0;
      i2 = ic0 + 12 * (lastc - 1);
      for (jy = ic0; jy <= i2; jy += 12) {
        ix = iv0;
        c = 0.0;
        j = (jy + lastv) - 1;
        for (ia = jy; ia <= j; ia++) {
          c += C_data[ia - 1] * C_data[ix - 1];
          ix++;
        }

        work_data[i] += c;
        i++;
      }
    }

    if (!(-tau == 0.0)) {
      i = ic0 - 1;
      jy = 0;
      for (j = 1; j <= lastc; j++) {
        if (work_data[jy] != 0.0) {
          c = work_data[jy] * -tau;
          ix = iv0;
          i2 = lastv + i;
          for (ia = i; ia < i2; ia++) {
            C_data[ia] += C_data[ix - 1] * c;
            ix++;
          }
        }

        jy++;
        i += 12;
      }
    }
  }
}

/*
 * 修改qp_ca_sl源文件得到，删除其他算法的分支。矩阵在matlab转c语言后会变成数组，对应关系为：
 *  数组的元素依次为矩阵第一列、第二列、。。。、直到最后一列。
 *  Wrapper used in the QP control allocation Simulink block.
 * Arguments    : const double arg[12]
 *                const double B[27]
 *                const double plim[18]
 *                const double rlim[18]
 *                double T
 *                const double Wv[9]
 *                const double Wu[81]
 *                const double ud[9]
 *                double imax
 *                double gam
 *                boolean_T only_plim
 *                double u[9]
 * Return Type  : void
 */
void qp_ca_mch(const double arg[12], const double B[27], const double plim[18],
               const double rlim[18], double T, const double Wv[9], const double
               Wu[81], const double ud[9], double imax, double gam, boolean_T
               only_plim, double u[9])
{
  int i;
  double umin[9];
  double umax[9];
  double W0;
  double u0;
  double b_W0[9];

  /*  Dimensions */
  /* B行数 */
  /* B列数 */
  /*  Extract nonconstant input arguments */
  /*  Overall position limits */
  if (only_plim) {
    memcpy(&umin[0], &plim[0], 9U * sizeof(double));
    memcpy(&umax[0], &plim[9], 9U * sizeof(double));
  } else {
    for (i = 0; i < 9; i++) {
      W0 = arg[3 + i] + rlim[i] * T;
      u0 = plim[i];
      if (u0 > W0) {
        W0 = u0;
      }

      umin[i] = W0;
      W0 = arg[3 + i] + rlim[9 + i] * T;
      u0 = plim[9 + i];
      if (u0 < W0) {
        W0 = u0;
      }

      umax[i] = W0;
    }
  }

  /*  u0 = (umin+umax)/2; */
  /*  W0 = zeros(m,1); */
  /*  u0=uprev; */
  /*    u = wls_alloc(B,v,umin,umax,Wv,Wu,ud,gam,u0,W0,imax); */
  for (i = 0; i < 9; i++) {
    b_W0[i] = 0.0;
    if (arg[3 + i] <= umin[i]) {
      b_W0[i] = -1.0;
    }

    if (arg[3 + i] >= umax[i]) {
      b_W0[i] = 1.0;
    }

    u[i] = arg[i + 3];
  }

  wls_alloc_mch(B, *(double (*)[3])&arg[0], umin, umax, Wv, Wu, ud, gam, u, b_W0,
                imax);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void qp_ca_mch_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void qp_ca_mch_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for qp_ca_mch.c
 *
 * [EOF]
 */
