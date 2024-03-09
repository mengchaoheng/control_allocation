/*
 * File: wls_alloc_m.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 18-Aug-2020 19:20:00
 */

/* Include Files */
#include <string.h>
#include <math.h>
#include "wls_alloc_m.h"

/* Function Declarations */
static void LSQFromQR(const double A_data[], const int A_size[2], const double
                      tau_data[], const int jpvt_data[], double B[12], int rankA,
                      double Y_data[], int Y_size[1]);
static boolean_T any(const boolean_T x_data[], const int x_size[1]);
static void qrsolve(const double A_data[], const int A_size[2], const double B
                    [12], double Y_data[], int Y_size[1]);
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
      for (i = loop_ub + 1; i < 12; i++) {
        wj += A_data[i + A_size[0] * loop_ub] * B[i];
      }

      wj *= tau_data[loop_ub];
      if (wj != 0.0) {
        B[loop_ub] -= wj;
        for (i = loop_ub + 1; i < 12; i++) {
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
 * Arguments    : const boolean_T x_data[]
 *                const int x_size[1]
 * Return Type  : boolean_T
 */
static boolean_T any(const boolean_T x_data[], const int x_size[1])
{
  boolean_T y;
  int ix;
  boolean_T exitg1;
  y = false;
  ix = 1;
  exitg1 = false;
  while ((!exitg1) && (ix <= x_size[0])) {
    if (x_data[ix - 1]) {
      y = true;
      exitg1 = true;
    } else {
      ix++;
    }
  }

  return y;
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
  static double b_A_data[108];
  int b_n;
  int jpvt_data[9];
  int yk;
  int k;
  static double b_B[12];
  static double tau_data[9];
  static double work_data[9];
  int i;
  double smax;
  double s;
  int i_i;
  int nmi;
  double absxk;
  static double vn1_data[9];
  static double vn2_data[9];
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
      i1 = ic0 + 12 * (lastc - 1);
      for (jy = ic0; jy <= i1; jy += 12) {
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
        i += 12;
      }
    }
  }
}

/*
 * WLS_ALLOC - Control allocation using weighted least squares.
 *
 *   [u,W,iter] = wls_alloc(B,v,umin,umax,[Wv,Wu,ud,gamma,u0,W0,imax])
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
 *  B=[-4.6718         0    4.6718         0   -4.6718         0    4.6718         0    3.2055;
 *           0  -14.9632         0   14.9632         0  -14.9632         0   14.9632         0;
 *      4.2635   -6.9597    4.2635   15.4868    4.2635   15.4868    4.2635   -6.9597         0];
 *  Number of variables
 *  k_TS=9.9796018325697625989171178675552e-6;
 *  speed=1065;
 *  kc=3.157;
 *  l_1=0.17078793-0.09;% roll,pitch
 *  l_2=0.175;% T
 *  l_3=0.06647954;% yaw1 yaw3
 *  l_4=0.06647954+0.175;% yaw4
 *  l_5=0.175-0.06647954;% yaw2
 *  I_x=0.054593;
 *  I_y=0.017045;
 *  I_z=0.049226;
 *  I=[I_x 0 0;0 I_y 0;0 0 I_z];
 *  % =========================10 sureface=================================
 *  L=[-l_1 0 l_1 0 -l_1 0 l_1 0 l_2;
 *      0 -l_1 0 l_1 0 -l_1 0 l_1 0;
 *      l_3 -l_5 l_3 l_4 l_3 l_4 l_3 -l_5 0];
 *  F=diag([kc kc kc kc kc kc kc kc 1],0);%
 *  B=I\L*F;
 *  B=[-4.6718         0    4.6718         0   -4.6718         0    4.6718         0    3.2055;
 *           0  -14.9632         0   14.9632         0  -14.9632         0   14.9632         0;
 *      4.2635   -6.9597    4.2635   15.4868    4.2635   15.4868    4.2635   -6.9597         0];
 *  操纵面维数
 * Arguments    : const double B[27]
 *                const double v[3]
 *                double u[9]
 *                const double p_limits[2]
 *                const double v_limits[2]
 *                double T
 *                const double Wv[9]
 *                const double Wu[81]
 *                const double ud[9]
 *                double imax
 *                double gam
 *                boolean_T only_plim
 * Return Type  : void
 */
void wls_alloc_m(const double B[27], const double v[3], double u[9], const
                 double p_limits[2], const double v_limits[2], double T, const
                 double Wv[9], const double Wu[81], const double ud[9], double
                 imax, double gam, boolean_T only_plim)
{
  static double dist[9];
  static double umin[9];
  int i;
  static double u_opt[9];
  static double umax[9];
  int k;
  double gam_sq;
  double u0;
  static double W[9];
  static double b_gam_sq[27];
  int aoffset;
  static double A[108];
  double c_gam_sq[3];
  static double d_gam_sq[12];
  static double b_A[12];
  boolean_T i_free[9];
  static double d[12];
  int iter;
  boolean_T exitg1;
  int trueCount;
  int A_free_size[2];
  int tmp_data[9];
  static double A_free_data[108];
  static double p_free_data[9];
  int p_free_size[1];
  int u_opt_size[1];
  int b_tmp_data[9];
  boolean_T b_u_opt[9];
  boolean_T u_opt_data[9];
  static double p[9];
  static double c_tmp_data[12];
  boolean_T y;
  int umin_size[1];
  int d_tmp_data[9];
  static double p_data[9];
  static double e_tmp_data[9];
  int tmp_size[1];
  boolean_T exitg2;
  int umax_size[1];
  int f_tmp_data[9];
  int idx;
  static double y_data[108];

  /* 是否仅含幅值约束 */
  if (only_plim) {
    umin[8] = -p_limits[1];
    for (i = 0; i < 8; i++) {
      umin[i] = -p_limits[0];
      umax[i] = p_limits[0];
    }

    umax[8] = p_limits[1];
  } else {
    dist[8] = -p_limits[1];
    for (i = 0; i < 8; i++) {
      dist[i] = -p_limits[0];
      u_opt[i] = -v_limits[0] * T;
    }

    u_opt[8] = -v_limits[1] * T;
    for (k = 0; k < 9; k++) {
      gam_sq = u[k] + u_opt[k];
      u0 = dist[k];
      if (u0 > gam_sq) {
        gam_sq = u0;
      }

      umin[k] = gam_sq;
    }

    dist[8] = p_limits[1];
    for (i = 0; i < 8; i++) {
      dist[i] = p_limits[0];
      u_opt[i] = v_limits[0] * T;
    }

    u_opt[8] = v_limits[1] * T;
    for (k = 0; k < 9; k++) {
      gam_sq = u[k] + u_opt[k];
      u0 = dist[k];
      if (u0 < gam_sq) {
        gam_sq = u0;
      }

      umax[k] = gam_sq;
    }
  }

  /* ========计算当前有效集============ */
  /*  coder.varsize('infeasible1',[m 1]); */
  /*  coder.varsize('infeasible2',[m 1]); */
  for (i = 0; i < 9; i++) {
    W[i] = 0.0;
    if (u[i] <= umin[i]) {
      W[i] = -1.0;
    }

    if (u[i] >= umax[i]) {
      W[i] = 1.0;
    }
  }

  /* ===============期望舵机位置======================== */
  /*  ud=zeros(m,1);% 可作为输入 */
  /* ================================ */
  /*  加权系数 */
  /*  gam=1e6; */
  /*  加权矩阵 */
  /*  Wv=eye(3); */
  /*  Wu=eye(m);  */
  /*  迭代次数上限 */
  /*  imax=100; */
  /*  Set default values of optional arguments   */
  gam_sq = sqrt(gam);
  for (i = 0; i < 3; i++) {
    for (k = 0; k < 9; k++) {
      b_gam_sq[i + 3 * k] = 0.0;
      for (aoffset = 0; aoffset < 3; aoffset++) {
        b_gam_sq[i + 3 * k] += gam_sq * Wv[i + 3 * aoffset] * B[aoffset + 3 * k];
      }
    }
  }

  for (i = 0; i < 9; i++) {
    for (k = 0; k < 3; k++) {
      A[k + 12 * i] = b_gam_sq[k + 3 * i];
    }

    memcpy(&A[i * 12 + 3], &Wu[i * 9], 9U * sizeof(double));
  }

  /*  A=1.0e+04 *[ */
  /*       */
  /*   */
  /*     -0.4672         0    0.4672         0   -0.4672         0    0.4672         0    0.3206; */
  /*           0   -1.4963         0    1.4963         0   -1.4963         0    1.4963         0; */
  /*      0.4264   -0.6960    0.4264    1.5487    0.4264    1.5487    0.4264   -0.6960         0; */
  /*      0.0001         0         0         0         0         0         0         0         0; */
  /*           0    0.0001         0         0         0         0         0         0         0; */
  /*           0         0    0.0001         0         0         0         0         0         0; */
  /*           0         0         0    0.0001         0         0         0         0         0; */
  /*           0         0         0         0    0.0001         0         0         0         0; */
  /*           0         0         0         0         0    0.0001         0         0         0; */
  /*           0         0         0         0         0         0    0.0001         0         0; */
  /*           0         0         0         0         0         0         0    0.0001         0; */
  /*           0         0         0         0         0         0         0         0    0.0001 */
  /*  ]; */
  /*  Initial residual. */
  for (i = 0; i < 3; i++) {
    c_gam_sq[i] = 0.0;
    for (k = 0; k < 3; k++) {
      c_gam_sq[i] += gam_sq * Wv[i + 3 * k] * v[k];
    }
  }

  for (i = 0; i < 9; i++) {
    u_opt[i] = 0.0;
    for (k = 0; k < 9; k++) {
      u_opt[i] += Wu[i + 9 * k] * ud[k];
    }
  }

  for (i = 0; i < 3; i++) {
    d_gam_sq[i] = c_gam_sq[i];
  }

  memcpy(&d_gam_sq[3], &u_opt[0], 9U * sizeof(double));
  for (i = 0; i < 12; i++) {
    b_A[i] = 0.0;
    for (k = 0; k < 9; k++) {
      b_A[i] += A[i + 12 * k] * u[k];
    }

    d[i] = d_gam_sq[i] - b_A[i];
  }

  /*  Determine indeces of free variables. */
  for (i = 0; i < 9; i++) {
    i_free[i] = (W[i] == 0.0);
  }

  /*  Iterate until optimum is found or maximum number of iterations */
  /*  is reached. */
  /*    A_free=0; */
  /*    p_free=0; */
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

    k = 0;
    for (i = 0; i < 9; i++) {
      if (i_free[i]) {
        tmp_data[k] = i + 1;
        k++;
      }
    }

    A_free_size[0] = 12;
    A_free_size[1] = trueCount;
    for (i = 0; i < trueCount; i++) {
      for (k = 0; k < 12; k++) {
        A_free_data[k + 12 * i] = A[k + 12 * (tmp_data[i] - 1)];
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
    k = 0;

    /*  ---------------------------- */
    /*   Is the new point feasible? */
    /*  ---------------------------- */
    aoffset = 0;
    for (i = 0; i < 9; i++) {
      gam_sq = 0.0;
      if (i_free[i]) {
        gam_sq = p_free_data[k];
        k++;
      }

      u_opt[i] = u[i] + gam_sq;
      if (i_free[i]) {
        aoffset++;
      }

      p[i] = gam_sq;
    }

    k = 0;
    for (i = 0; i < 9; i++) {
      if (i_free[i]) {
        b_tmp_data[k] = i + 1;
        k++;
      }

      b_u_opt[i] = ((u_opt[i] < umin[i]) || (u_opt[i] > umax[i]));
    }

    u_opt_size[0] = aoffset;
    for (i = 0; i < aoffset; i++) {
      u_opt_data[i] = b_u_opt[b_tmp_data[i] - 1];
    }

    if (!any(u_opt_data, u_opt_size)) {
      /*  ---------------------------- */
      /*   Yes, check for optimality. */
      /*  ---------------------------- */
      /*  Update point and residual. */
      memcpy(&u[0], &u_opt[0], 9U * sizeof(double));
      if ((trueCount == 1) || (p_free_size[0] == 1)) {
        for (i = 0; i < 12; i++) {
          for (k = 0; k < 1; k++) {
            c_tmp_data[i] = 0.0;
            for (aoffset = 0; aoffset < trueCount; aoffset++) {
              c_tmp_data[i] += A_free_data[i + 12 * aoffset] *
                p_free_data[aoffset];
            }
          }
        }
      } else {
        memset(&c_tmp_data[0], 0, 12U * sizeof(double));
        for (k = 0; k < trueCount; k++) {
          if (p_free_data[k] != 0.0) {
            aoffset = k * 12;
            for (i = 0; i < 12; i++) {
              c_tmp_data[i] += p_free_data[k] * A_free_data[aoffset + i];
            }
          }
        }
      }

      for (i = 0; i < 12; i++) {
        d[i] -= c_tmp_data[i];
      }

      /*  Compute Lagrangian multipliers. */
      for (i = 0; i < 9; i++) {
        u_opt[i] = 0.0;
        for (k = 0; k < 12; k++) {
          u_opt[i] += A[k + 12 * i] * d[k];
        }

        dist[i] = W[i] * u_opt[i];
      }

      /*  Are all lambda non-negative? */
      y = true;
      k = 1;
      exitg2 = false;
      while ((!exitg2) && (k < 10)) {
        if (!(dist[k - 1] >= -2.2204460492503131E-16)) {
          y = false;
          exitg2 = true;
        } else {
          k++;
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
        gam_sq = dist[0];
        idx = 0;
        for (k = 0; k < 8; k++) {
          if (gam_sq > dist[k + 1]) {
            gam_sq = dist[k + 1];
            idx = k + 1;
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
        dist[i] = 1.0;
        y = (p[i] < 0.0);
        b_u_opt[i] = (p[i] > 0.0);
        if (i_free[i] && y) {
          aoffset++;
        }

        u_opt_data[i] = y;
      }

      k = 0;
      for (i = 0; i < 9; i++) {
        if (i_free[i] && u_opt_data[i]) {
          d_tmp_data[k] = i + 1;
          k++;
        }
      }

      umin_size[0] = aoffset;
      for (i = 0; i < aoffset; i++) {
        u_opt[i] = umin[d_tmp_data[i] - 1] - u[d_tmp_data[i] - 1];
      }

      for (i = 0; i < aoffset; i++) {
        p_data[i] = p[d_tmp_data[i] - 1];
      }

      rdivide(u_opt, umin_size, p_data, e_tmp_data, tmp_size);
      k = 0;
      aoffset = 0;
      for (i = 0; i < 9; i++) {
        if (i_free[i] && u_opt_data[i]) {
          dist[i] = e_tmp_data[k];
          k++;
        }

        if (i_free[i] && b_u_opt[i]) {
          aoffset++;
        }
      }

      k = 0;
      for (i = 0; i < 9; i++) {
        if (i_free[i] && b_u_opt[i]) {
          f_tmp_data[k] = i + 1;
          k++;
        }
      }

      umax_size[0] = aoffset;
      for (i = 0; i < aoffset; i++) {
        u_opt[i] = umax[f_tmp_data[i] - 1] - u[f_tmp_data[i] - 1];
      }

      for (i = 0; i < aoffset; i++) {
        p_data[i] = p[f_tmp_data[i] - 1];
      }

      rdivide(u_opt, umax_size, p_data, e_tmp_data, tmp_size);
      k = 0;
      for (i = 0; i < 9; i++) {
        if (i_free[i] && b_u_opt[i]) {
          dist[i] = e_tmp_data[k];
          k++;
        }
      }

      /*  Proportion of p to travel */
      gam_sq = dist[0];
      idx = 0;
      for (k = 0; k < 8; k++) {
        if (gam_sq > dist[k + 1]) {
          gam_sq = dist[k + 1];
          idx = k + 1;
        }
      }

      /*  Update point and residual. */
      for (i = 0; i < 9; i++) {
        u[i] += gam_sq * p[i];
      }

      k = 12 * trueCount;
      for (i = 0; i < k; i++) {
        y_data[i] = A_free_data[i] * gam_sq;
      }

      if ((trueCount == 1) || (p_free_size[0] == 1)) {
        for (i = 0; i < 12; i++) {
          for (k = 0; k < 1; k++) {
            c_tmp_data[i] = 0.0;
            for (aoffset = 0; aoffset < trueCount; aoffset++) {
              c_tmp_data[i] += y_data[i + 12 * aoffset] * p_free_data[aoffset];
            }
          }
        }
      } else {
        memset(&c_tmp_data[0], 0, 12U * sizeof(double));
        for (k = 0; k < trueCount; k++) {
          if (p_free_data[k] != 0.0) {
            aoffset = k * 12;
            for (i = 0; i < 12; i++) {
              c_tmp_data[i] += p_free_data[k] * y_data[aoffset + i];
            }
          }
        }
      }

      for (i = 0; i < 12; i++) {
        d[i] -= c_tmp_data[i];
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
 * Arguments    : void
 * Return Type  : void
 */
void wls_alloc_m_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void wls_alloc_m_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for wls_alloc_m.c
 *
 * [EOF]
 */
