/*
 * File: wls_alloc_mch.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 08-Aug-2020 14:31:14
 */

/* Include Files */
#include <string.h>
#include <math.h>
#include "wls_alloc_mch.h"
sWls wls;
/* Function Declarations */
static void LSQFromQR(const double A_data[], const int A_size[2], const double
                      tau_data[], const int jpvt_data[], double B[11], int rankA,
                      double Y_data[], int Y_size[1]);
static void qrsolve(const double A_data[], const int A_size[2], const double B
                    [11], double Y_data[], int Y_size[1]);
static int rankFromQR(const double A_data[], const int A_size[2]);
static void rdivide(const double x_data[], const int x_size[1], const double
                    y_data[], double z_data[], int z_size[1]);
static double rt_hypotd(double u0, double u1);
static double xnrm2(int n, const double x_data[], int ix0);
static void xzlarf(int m, int n, int iv0, double tau, double C_data[], int ic0,
                   double work_data[]);

/* Function Definitions */

/*
 * Arguments    : const double A_data[]
 *                const int A_size[2]
 *                const double tau_data[]
 *                const int jpvt_data[]
 *                double B[11]
 *                int rankA
 *                double Y_data[]
 *                int Y_size[1]
 * Return Type  : void
 */
static void LSQFromQR(const double A_data[], const int A_size[2], const double
                      tau_data[], const int jpvt_data[], double B[11], int rankA,
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
      for (i = loop_ub + 1; i + 1 < 12; i++) {
        wj += A_data[i + A_size[0] * loop_ub] * B[i];
      }

      wj *= tau_data[loop_ub];
      if (wj != 0.0) {
        B[loop_ub] -= wj;
        for (i = loop_ub + 1; i + 1 < 12; i++) {
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
 *                const double B[11]
 *                double Y_data[]
 *                int Y_size[1]
 * Return Type  : void
 */
static void qrsolve(const double A_data[], const int A_size[2], const double B
                    [11], double Y_data[], int Y_size[1])
{
  int b_A_size[2];
  int n;
  double b_A_data[88];
  int b_n;
  int jpvt_data[8];
  int yk;
  int k;
  double b_B[11];
  double tau_data[8];
  double work_data[8];
  int i;
  double smax;
  double s;
  int i_i;
  int nmi;
  double absxk;
  double vn1_data[8];
  double vn2_data[8];
  double t;
  int pvt;
  b_A_size[0] = 11;
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
      for (n = k; n <= k + 10; n++) {
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
      k += 11;
    }

    for (i = 0; i < b_n; i++) {
      i_i = i + i * 11;
      nmi = b_n - i;
      if (nmi < 1) {
        n = -1;
      } else {
        n = 0;
        if (nmi > 1) {
          yk = i;
          smax = fabs(vn1_data[i]);
          for (k = 0; k + 2 <= nmi; k++) {
            yk++;
            s = fabs(vn1_data[yk]);
            if (s > smax) {
              n = k + 1;
              smax = s;
            }
          }
        }
      }

      pvt = i + n;
      if (pvt + 1 != i + 1) {
        yk = 11 * pvt;
        n = 11 * i;
        for (k = 0; k < 11; k++) {
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
      smax = xnrm2(10 - i, b_A_data, i_i + 2);
      if (smax != 0.0) {
        smax = rt_hypotd(b_A_data[i_i], smax);
        if (b_A_data[i_i] >= 0.0) {
          smax = -smax;
        }

        if (fabs(smax) < 1.0020841800044864E-292) {
          yk = 0;
          n = (i_i - i) + 11;
          do {
            yk++;
            for (k = i_i + 1; k < n; k++) {
              b_A_data[k] *= 9.9792015476736E+291;
            }

            smax *= 9.9792015476736E+291;
            s *= 9.9792015476736E+291;
          } while (!(fabs(smax) >= 1.0020841800044864E-292));

          smax = rt_hypotd(s, xnrm2(10 - i, b_A_data, i_i + 2));
          if (s >= 0.0) {
            smax = -smax;
          }

          absxk = (smax - s) / smax;
          s = 1.0 / (s - smax);
          n = (i_i - i) + 11;
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
          n = (i_i - i) + 11;
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
        xzlarf(11 - i, nmi - 1, i_i + 1, tau_data[i], b_A_data, (i + (i + 1) *
                11) + 1, work_data);
        b_A_data[i_i] = s;
      }

      for (yk = i + 1; yk < b_n; yk++) {
        if (vn1_data[yk] != 0.0) {
          smax = fabs(b_A_data[i + 11 * yk]) / vn1_data[yk];
          smax = 1.0 - smax * smax;
          if (smax < 0.0) {
            smax = 0.0;
          }

          s = vn1_data[yk] / vn2_data[yk];
          s = smax * (s * s);
          if (s <= 1.4901161193847656E-8) {
            vn1_data[yk] = xnrm2(10 - i, b_A_data, (i + 11 * yk) + 2);
            vn2_data[yk] = vn1_data[yk];
          } else {
            vn1_data[yk] *= sqrt(smax);
          }
        }
      }
    }
  }

  memcpy(&b_B[0], &B[0], 11U * sizeof(double));
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
    tol = 11.0 * fabs(A_data[0]) * 2.2204460492503131E-16;
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
  int i1;
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
      i = ic0 + (lastc - 1) * 11;
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
      i1 = ic0 + 11 * (lastc - 1);
      for (jy = ic0; jy <= i1; jy += 11) {
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
          i1 = lastv + i;
          for (ia = i; ia < i1; ia++) {
            C_data[ia] += C_data[ix - 1] * c;
            ix++;
          }
        }

        jy++;
        i += 11;
      }
    }
  }
}

/*
 * function [u] = wls_alloc_mch(v, u, umin, umax)
 *   [u] = wls_alloc_mch(v,u,p_limits,v_limits)
 *  WLS_ALLOC - Control allocation using weighted least squares.
 *   [u,W,iter] = wls_alloc(B,v,umin,umax,[Wv,Wu,ud,gamma,u0,W0,imax])
 *  Solves the weighted, bounded least-squares problem
 *    min ||Wu(u-ud)||^2 + gamma ||Wv(Bu-v)||^2
 *    subj. to  umin <= u <= umax
 *  using an active set method.
 *   Inputs:
 *   -------
 *  v     commanded virtual control (k x 1)
 *  u0    initial point (m x 1)
 *  W0    initial working set (m x 1) [empty]
 *  imax  max no. of iterations [100]
 *   Outputs:
 *   -------
 *  u     optimal control
 *  W     optimal active set
 *  iter  no. of iterations (= no. of changes in the working set + 1)
 *                             0 if u_i not saturated
 *  Working set syntax: W_i = -1 if u_i = umin_i
 *                            +1 if u_i = umax_i
 *  B         control effectiveness matrix (k x m).
 *  umin      lower position limits (m x 1).
 *  umax      upper position limits (m x 1).
 *  Wv        virtual control weighting matrix (k x k) [I].
 *  Wu        control weighting matrix (m x m) [I].
 *  gam       gamma weight (scalar) [1e6].
 *  ud        desired control (m x 1) [0].
 *  imax      maximum iterations.
 *  See also: WLSC_ALLOC, IP_ALLOC, FXP_ALLOC, QP_SIM.
 *  param of ducted fan
 *  使用力矩单位
 * Arguments    : const double v[3]
 *                double u[8]
 * Return Type  : void
 */
void wls_alloc_mch(const double v[3], double u[8])
{
  int i;
  double W[8];
  int aoffset;
  double a[11];
  double b_a[3];
  int partialTrueCount;
  double c_a[11];
  static const short d_a[9] = { 1000, 0, 0, 0, 1000, 0, 0, 0, 1000 };

  boolean_T i_free[8];
  double d[11];
  static const double e_a[88] = { -4671.7984908321587, 0.0, 4263.5174050298619,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -14963.185392197125,
    -6959.7182834274572, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    4671.7984908321587, 0.0, 4263.5174050298619, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 14963.185392197125, 15486.753093487183, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, -4671.7984908321587, 0.0, 4263.5174050298619, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -14963.185392197125, 15486.753093487183,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 4671.7984908321587, 0.0,
    4263.5174050298619, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    14963.185392197125, -6959.7182834274572, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0 };

  int iter;
  boolean_T exitg1;
  int trueCount;
  int A_free_size[2];
  int tmp_data[8];
  double A_free_data[88];
  double p_free_data[8];
  int p_free_size[1];
  int b_trueCount;
  double p;
  double u_opt[8];
  int b_tmp_data[8];
  boolean_T y;
  boolean_T x_data[8];
  boolean_T b_u_opt[8];
  double b_p[8];
  boolean_T exitg2;
  boolean_T bv0[8];
  int tmp_size[1];
  int c_tmp_data[8];
  double d_tmp_data[8];
  double p_data[8];
  double e_tmp_data[8];
  int b_tmp_size[1];
  static const double f_a[88] = { -4671.7984908321587, 0.0, 4671.7984908321587,
    0.0, -4671.7984908321587, 0.0, 4671.7984908321587, 0.0, 0.0,
    -14963.185392197125, 0.0, 14963.185392197125, 0.0, -14963.185392197125, 0.0,
    14963.185392197125, 4263.5174050298619, -6959.7182834274572,
    4263.5174050298619, 15486.753093487183, 4263.5174050298619,
    15486.753093487183, 4263.5174050298619, -6959.7182834274572, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0 };

  int c_tmp_size[1];
  int f_tmp_data[8];

  /* [rad  rad/s] */
  /* =============================================================== */
  /*  roll,pitch */
  /*  T */
  /*  yaw1 yaw3 */
  /*  yaw4 */
  /*  yaw2 */
  /*  =========================10 sureface================================= */
  /*  L=[-l_1 0 l_1 0 -l_1 0 l_1 0 l_2; */
  /*      0 -l_1 0 l_1 0 -l_1 0 l_1 0; */
  /*      l_3 -l_5 l_3 l_4 l_3 l_4 l_3 -l_5 0]; */
  /*  F=diag([kc kc kc kc kc kc kc kc 2*k_TS*speed],0);% */
  /*  */
  /* ======================================== */
  /*  */
  /* ========计算当前有效集============ */
  for (i = 0; i < 8; i++) {
    W[i] = 0.0;
    if (u[i] <= -0.3490658503988659) {
      W[i] = -1.0;
    }

    if (u[i] >= 0.3490658503988659) {
      W[i] = 1.0;
    }
  }

  /* ===============期望舵机位置======================== */
  /*  可作为输入 */
  /*  Number of variables. */
  /* ================================ */
  /*  加权系数 */
  /*  加权矩阵 */
  /* Wv1=Wv*0.5*y; */
  /*  迭代次数上限 */
  /*  Set default values of optional arguments. */
  /*  Initial residual. */
  for (aoffset = 0; aoffset < 3; aoffset++) {
    b_a[aoffset] = 0.0;
    for (partialTrueCount = 0; partialTrueCount < 3; partialTrueCount++) {
      b_a[aoffset] += (double)d_a[aoffset + 3 * partialTrueCount] *
        v[partialTrueCount];
    }

    a[aoffset] = b_a[aoffset];
  }

  memset(&a[3], 0, sizeof(double) << 3);
  for (aoffset = 0; aoffset < 11; aoffset++) {
    c_a[aoffset] = 0.0;
    for (partialTrueCount = 0; partialTrueCount < 8; partialTrueCount++) {
      c_a[aoffset] += e_a[aoffset + 11 * partialTrueCount] * u[partialTrueCount];
    }

    d[aoffset] = a[aoffset] - c_a[aoffset];
  }

  /*  Determine indeces of free variables. */
  for (i = 0; i < 8; i++) {
    i_free[i] = (W[i] == 0.0);
  }

  /*  Iterate until optimum is found or maximum number of iterations */
  /*  is reached. */
  iter = 0;
  exitg1 = false;
  while ((!exitg1) && (iter < 100)) {
    /*  ---------------------------------------- */
    /*   Compute optimal perturbation vector p. */
    /*  ---------------------------------------- */
    /*  Eliminate saturated variables. */
    trueCount = 0;
    for (i = 0; i < 8; i++) {
      if (i_free[i]) {
        trueCount++;
      }
    }

    partialTrueCount = 0;
    for (i = 0; i < 8; i++) {
      if (i_free[i]) {
        tmp_data[partialTrueCount] = i + 1;
        partialTrueCount++;
      }
    }

    A_free_size[0] = 11;
    A_free_size[1] = trueCount;
    for (aoffset = 0; aoffset < trueCount; aoffset++) {
      for (partialTrueCount = 0; partialTrueCount < 11; partialTrueCount++) {
        A_free_data[partialTrueCount + 11 * aoffset] = e_a[partialTrueCount + 11
          * (tmp_data[aoffset] - 1)];
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
    b_trueCount = 0;
    for (i = 0; i < 8; i++) {
      p = 0.0;
      if (i_free[i]) {
        p = p_free_data[partialTrueCount];
        partialTrueCount++;
      }

      u_opt[i] = u[i] + p;
      if (i_free[i]) {
        b_trueCount++;
      }

      b_p[i] = p;
    }

    partialTrueCount = 0;
    for (i = 0; i < 8; i++) {
      if (i_free[i]) {
        b_tmp_data[partialTrueCount] = i + 1;
        partialTrueCount++;
      }

      b_u_opt[i] = ((u_opt[i] < -0.3490658503988659) || (u_opt[i] >
        0.3490658503988659));
    }

    for (aoffset = 0; aoffset < b_trueCount; aoffset++) {
      x_data[aoffset] = b_u_opt[b_tmp_data[aoffset] - 1];
    }

    y = false;
    partialTrueCount = 1;
    exitg2 = false;
    while ((!exitg2) && (partialTrueCount <= b_trueCount)) {
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
      memcpy(&u[0], &u_opt[0], sizeof(double) << 3);
      if ((trueCount == 1) || (p_free_size[0] == 1)) {
        for (aoffset = 0; aoffset < 11; aoffset++) {
          a[aoffset] = 0.0;
          for (partialTrueCount = 0; partialTrueCount < trueCount;
               partialTrueCount++) {
            a[aoffset] += A_free_data[aoffset + 11 * partialTrueCount] *
              p_free_data[partialTrueCount];
          }
        }
      } else {
        memset(&a[0], 0, 11U * sizeof(double));
        for (partialTrueCount = 0; partialTrueCount < trueCount;
             partialTrueCount++) {
          if (p_free_data[partialTrueCount] != 0.0) {
            aoffset = partialTrueCount * 11;
            for (i = 0; i < 11; i++) {
              a[i] += p_free_data[partialTrueCount] * A_free_data[aoffset + i];
            }
          }
        }
      }

      for (aoffset = 0; aoffset < 11; aoffset++) {
        d[aoffset] -= a[aoffset];
      }

      /*  判别推出循环条件 */
      /*          if norm(p)<=eps */
      /*  Compute Lagrangian multipliers. */
      for (aoffset = 0; aoffset < 8; aoffset++) {
        b_p[aoffset] = 0.0;
        for (partialTrueCount = 0; partialTrueCount < 11; partialTrueCount++) {
          b_p[aoffset] += f_a[aoffset + (partialTrueCount << 3)] *
            d[partialTrueCount];
        }

        u_opt[aoffset] = W[aoffset] * b_p[aoffset];
      }

      /*  Are all lambda non-negative? */
      y = true;
      partialTrueCount = 1;
      exitg2 = false;
      while ((!exitg2) && (partialTrueCount < 9)) {
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
        p = u_opt[0];
        b_trueCount = 0;
        for (partialTrueCount = 0; partialTrueCount < 7; partialTrueCount++) {
          if (p > u_opt[partialTrueCount + 1]) {
            p = u_opt[partialTrueCount + 1];
            b_trueCount = partialTrueCount + 1;
          }
        }

        W[b_trueCount] = 0.0;
        i_free[b_trueCount] = true;

        /*          end */
        iter++;
      }
    } else {
      /*  --------------------------------------- */
      /*   No, find primary bounding constraint. */
      /*  --------------------------------------- */
      /*  Compute distances to the different boundaries. Since alpha < 1 is */
      /*  the maximum step length, initiate with ones. */
      b_trueCount = 0;
      for (i = 0; i < 8; i++) {
        u_opt[i] = 1.0;
        y = (b_p[i] < 0.0);
        bv0[i] = (b_p[i] > 0.0);
        if (i_free[i] && y) {
          b_trueCount++;
        }

        b_u_opt[i] = y;
      }

      partialTrueCount = 0;
      for (i = 0; i < 8; i++) {
        if (i_free[i] && b_u_opt[i]) {
          c_tmp_data[partialTrueCount] = i + 1;
          partialTrueCount++;
        }
      }

      tmp_size[0] = b_trueCount;
      for (aoffset = 0; aoffset < b_trueCount; aoffset++) {
        d_tmp_data[aoffset] = -0.3490658503988659 - u[c_tmp_data[aoffset] - 1];
      }

      for (aoffset = 0; aoffset < b_trueCount; aoffset++) {
        p_data[aoffset] = b_p[c_tmp_data[aoffset] - 1];
      }

      rdivide(d_tmp_data, tmp_size, p_data, e_tmp_data, b_tmp_size);
      partialTrueCount = 0;
      b_trueCount = 0;
      for (i = 0; i < 8; i++) {
        if (i_free[i] && b_u_opt[i]) {
          u_opt[i] = e_tmp_data[partialTrueCount];
          partialTrueCount++;
        }

        if (i_free[i] && bv0[i]) {
          b_trueCount++;
        }
      }

      partialTrueCount = 0;
      for (i = 0; i < 8; i++) {
        if (i_free[i] && bv0[i]) {
          f_tmp_data[partialTrueCount] = i + 1;
          partialTrueCount++;
        }
      }

      c_tmp_size[0] = b_trueCount;
      for (aoffset = 0; aoffset < b_trueCount; aoffset++) {
        d_tmp_data[aoffset] = 0.3490658503988659 - u[f_tmp_data[aoffset] - 1];
      }

      for (aoffset = 0; aoffset < b_trueCount; aoffset++) {
        p_data[aoffset] = b_p[f_tmp_data[aoffset] - 1];
      }

      rdivide(d_tmp_data, c_tmp_size, p_data, e_tmp_data, b_tmp_size);
      partialTrueCount = 0;
      for (i = 0; i < 8; i++) {
        if (i_free[i] && bv0[i]) {
          u_opt[i] = e_tmp_data[partialTrueCount];
          partialTrueCount++;
        }
      }

      /*  Proportion of p to travel */
      p = u_opt[0];
      b_trueCount = 0;
      for (partialTrueCount = 0; partialTrueCount < 7; partialTrueCount++) {
        if (p > u_opt[partialTrueCount + 1]) {
          p = u_opt[partialTrueCount + 1];
          b_trueCount = partialTrueCount + 1;
        }
      }

      /*  Update point and residual. */
      for (aoffset = 0; aoffset < 8; aoffset++) {
        u[aoffset] += p * b_p[aoffset];
      }

      /*  Update point and residual. */
      partialTrueCount = 11 * trueCount - 1;
      for (aoffset = 0; aoffset <= partialTrueCount; aoffset++) {
        A_free_data[aoffset] *= p;
      }

      if ((trueCount == 1) || (p_free_size[0] == 1)) {
        for (aoffset = 0; aoffset < 11; aoffset++) {
          a[aoffset] = 0.0;
          for (partialTrueCount = 0; partialTrueCount < trueCount;
               partialTrueCount++) {
            a[aoffset] += A_free_data[aoffset + 11 * partialTrueCount] *
              p_free_data[partialTrueCount];
          }
        }
      } else {
        memset(&a[0], 0, 11U * sizeof(double));
        for (partialTrueCount = 0; partialTrueCount < trueCount;
             partialTrueCount++) {
          if (p_free_data[partialTrueCount] != 0.0) {
            aoffset = partialTrueCount * 11;
            for (i = 0; i < 11; i++) {
              a[i] += p_free_data[partialTrueCount] * A_free_data[aoffset + i];
            }
          }
        }
      }

      for (aoffset = 0; aoffset < 11; aoffset++) {
        d[aoffset] -= a[aoffset];
      }

      /*  Add corresponding constraint to working set. */
      p = b_p[b_trueCount];
      if (b_p[b_trueCount] < 0.0) {
        p = -1.0;
      } else {
        if (b_p[b_trueCount] > 0.0) {
          p = 1.0;
        }
      }

      W[b_trueCount] = p;
      i_free[b_trueCount] = false;
      iter++;
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void wls_alloc_mch_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void wls_alloc_mch_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for wls_alloc_mch.c
 *
 * [EOF]
 */
