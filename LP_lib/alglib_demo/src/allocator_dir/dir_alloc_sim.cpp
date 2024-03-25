/*
 * File: dir_alloc_sim.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 25-Apr-2021 15:48:18
 */

/* Include Files */
#include <string.h>
#include "dir_alloc_sim.h"

/* Function Declarations */
static void Simplex_loop_C(double B[13], double A[260], double b[13], double c
  [20], double *z, double x[20], double *iters);
static void inv_mch(double A[169], double col, double A_inv[169]);

/* Function Definitions */

/*
 * (c) mengchaoheng
 *  �������޽������
 *  Last edited 2019-11
 *    min z=c*x   subj. to  A*x (=�� >=�� <=) b
 *    x
 *     %% Initialization
 * Arguments    : double B[13]
 *                double A[260]
 *                double b[13]
 *                double c[20]
 *                double *z
 *                double x[20]
 *                double *iters
 * Return Type  : void
 */
static void Simplex_loop_C(double B[13], double A[260], double b[13], double c
  [20], double *z, double x[20], double *iters)
{
  int exitg1;
  boolean_T flag;
  int e;
  int L;
  int i;
  boolean_T exitg2;
  double MIN;
  double delta;
  int k;
  double a_je;

  /*  Iterate through simplex algorithm main loop */
  memset(&x[0], 0, 20U * sizeof(double));
  *iters = 0.0;
  *z = 0.0;

  /*      [m,n] = size(A); */
  /*      while ~all(c>=0)                      % 3.~isempty(c(c(N)<0)) */
  /*      e = find(c < 0, 1, 'first'); % ������������    % 4. e = N(find(c(N)<0,1)) */
  do {
    exitg1 = 0;
    flag = false;
    e = -1;
    L = -1;
    i = 0;
    exitg2 = false;
    while ((!exitg2) && (i < 20)) {
      if (c[i] < 0.0) {
        flag = true;
        e = i;
        exitg2 = true;
      } else {
        i++;
      }
    }

    if (flag) {
      /*              a_ie=A(:,e); */
      /*              ip=a_ie>(1/tol); */
      /*              delta=tol*ones(m,1); */
      /*              if ~isempty(ip) */
      /*                  delta(ip)=b(ip)./a_ie(ip); */
      /*              end */
      MIN = 1.0E+8;
      for (i = 0; i < 13; i++) {
        if (A[i + 13 * e] > 1.0E-8) {
          delta = b[i] / A[i + 13 * e];
        } else {
          delta = 1.0E+8;
        }

        if (delta < MIN) {
          L = i;
          MIN = delta;
        }
      }

      /*              [~,L]=min(delta);%ѡ����� (�����B�����е�������) */
      /*          li = B(L);    % �����������                */
      /*              if delta(L) >= tol     */
      if (MIN >= 1.0E+8) {
        exitg1 = 1;
      } else {
        /*  ��ʱһ����һ��L */
        /*  (c) mengchaoheng */
        /*  Last edited 2019-11 */
        /*    min z=c*x   subj. to  A*x (=�� >=�� <=) b */
        /*    x  */
        /*     %% Compute the coefficients of the equation for new basic variabLe x_e. */
        /*      [m, n] = size(A); */
        /*  row of Leaving var    L = find(B==li,1); */
        /*  Perform pivot operation, exchanging L-row with e-coLumn variabLe */
        b[L] /= A[L + 13 * e];

        /*  4. */
        MIN = A[L + 13 * e];
        for (i = 0; i < 20; i++) {
          A[L + 13 * i] /= MIN;
        }

        /*     %% Compute the coefficients of the remaining constraints. */
        /*      i=[1:L-1 L+1:m];     %  i = find(B~=li); */
        /*      if ~isempty(i) */
        /*          b(i) = b(i) - A(i,e)*b(L);                                                */
        /*          A(i,1:n) = A(i,1:n) - A(i,e)*A(L,1:n);	 */
        /*      end */
        MIN = b[L];
        for (i = 0; i < 13; i++) {
          delta = b[i];
          if (1 + i != L + 1) {
            delta = b[i] - A[i + 13 * e] * MIN;
            a_je = A[i + 13 * e];
            for (k = 0; k < 20; k++) {
              A[i + 13 * k] -= a_je * A[L + 13 * k];
            }
          }

          b[i] = delta;
        }

        /*     %% Compute the objective function */
        MIN = c[e];
        *z -= c[e] * b[L];
        for (k = 0; k < 20; k++) {
          c[k] -= MIN * A[L + 13 * k];
        }

        /*      c(1:n) = c(1:n) - c_e * A(L,1:n);       */
        /*     %% Compute new sets of basic and nonbasic variabLes. */
        /*  N(find(N==e,1)) = li;  */
        B[L] = (double)e + 1.0;

        /*   B(find(B==li,1)) = e; */
        /* �����������г����б任 */
        (*iters)++;
      }
    } else {
      for (i = 0; i < 13; i++) {
        x[(int)B[i] - 1] = b[i];
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

/*
 * �Ծ�����г����б任������
 * Arguments    : double A[169]
 *                double col
 *                double A_inv[169]
 * Return Type  : void
 */
static void inv_mch(double A[169], double col, double A_inv[169])
{
  int k;
  int i;
  double div_i;
  int jj;

  /*  function Ad_eye=inv_mvh(B_inv,Ad) */
  /*  Ad_eye=B_inv\Ad;% ���� */
  /*  [row, col] = size(A); */
  /*  BΪ��λ���� */
  memset(&A_inv[0], 0, 169U * sizeof(double));
  for (k = 0; k < 13; k++) {
    A_inv[k + 13 * k] = 1.0;
  }

  for (i = 0; i < 13; i++) {
    /*  ���ν��Խ��е�Ԫ�ع�һ�� */
    div_i = A[i + 13 * i];
    for (k = 0; k < (int)col; k++) {
      A[i + 13 * k] /= div_i;
      A_inv[i + 13 * k] /= div_i;
    }

    for (k = 0; k < 13; k++) {
      div_i = -A[k + 13 * i] / A[i + 13 * i];
      if (1 + i == 1 + k) {
        div_i = 0.0;
      }

      /*  �����б任 */
      for (jj = 0; jj < (int)col; jj++) {
        A[k + 13 * jj] += div_i * A[i + 13 * jj];
        A_inv[k + 13 * jj] += div_i * A_inv[i + 13 * jj];
      }
    }
  }
}

/*
 * (c) mengchaoheng
 *  Last edited 2019-11
 *    min z=c*x   subj. to  A*x (=�� >=�� <=) b
 *    x
 *  ԭ����
 *  Performs direct control allocation by solving the LP
 *    max z=a   subj. to  Bu = av
 *    a,u               umin <= u <= umax
 *  If a > 1, set u = u/a.
 *  Note: This function has not been optimized for speed.
 *   Inputs:
 *   -------
 *  B     control effectiveness matrix (k x m)
 *  v     commanded virtual control (k x 1)
 *  umin  lower position limits (m x 1)
 *  umax  upper position limits (m x 1)
 *   Outputs:
 *   -------
 *  u     optimal control (m x 1)
 *  a     scaling factor
 *  �����
 *    min z=[0 -1]x   subj. to  [B -v]x = 0
 *    x                       [I 0;-I 0]x <= [umax; -umin]
 *    ���� x=[u; a]
 *  ��Ӧ��͹�Ż���p139,��Ϊ
 *    min z=c*x   subj. to  Aeq*x = beq
 *    x                     G*x <= h
 *  �ϲ�
 *    min z=c*x   subj. to  [Aeq; G]*x (=��<=) [beq;h]
 *    x
 *  ��֤x>=0������
 *    min z=[c -c]*X   subj. to  [Aeq -Aeq;G -G]*X (=��<=) [beq;h]
 *     X
 *  ���� X=[x^+; x^-]
 * Arguments    : const double v[3]
 *                const double umin[4]
 *                const double umax[4]
 *                double u[4]
 *                double *z
 *                double *iters
 * Return Type  : void
 */
void dir_alloc_sim(const double v[3], const double umin[4], const double umax[4],
                   double u[4], double *z, double *iters)
{
  int i;
  int i0;
  static double Aeq[15];
  static const double dv0[12] = { -0.5, 0.0, 0.25, 0.0, -0.5, 0.25, 0.5, 0.0,
    0.25, 0.0, 0.5, 0.25 };

  static double Ad[130];
  static const signed char iv0[100] = { 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1 };

  static double b_Ad[169];
  static double P_inv[169];
  static const signed char iv1[100] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1 };

  double dv1[13];
  static const signed char iv2[13] = { 1, 2, 3, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 20 };

  static double Ad_eye[130];
  int i1;
  double dv2[13];
  static double b_Ad_eye[260];
  static double dv3[20];
  static const double dv4[20] = { 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static double x[20];
  for (i = 0; i < 4; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      Aeq[i0 + 3 * i] = dv0[i0 + 3 * i];
    }
  }

  for (i = 0; i < 3; i++) {
    Aeq[12 + i] = -v[i];
  }

  /* b������Թ滮 */
  /*  �������Թ滮��׼�� */
  /*  Convert free variables to positively constrained variables */
  for (i = 0; i < 5; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      Ad[i0 + 13 * i] = Aeq[i0 + 3 * i];
      Ad[i0 + 13 * (i + 5)] = -Aeq[i0 + 3 * i];
    }
  }

  for (i = 0; i < 10; i++) {
    for (i0 = 0; i0 < 10; i0++) {
      Ad[(i0 + 13 * i) + 3] = iv0[i0 + 10 * i];
    }
  }

  /*  Ad=[B -v -B v; eye(5) -eye(5);-eye(5) eye(5)]; */
  /*  Ad=[-0.5000         0    0.5000         0         0    0.5000         0   -0.5000         0         0; */
  /*            0   -0.5000         0    0.5000         0         0    0.5000         0   -0.5000         0; */
  /*       0.2500    0.2500    0.2500    0.2500         0   -0.2500   -0.2500   -0.2500   -0.2500   	    0; */
  /*       1.0000         0         0         0         0   -1.0000         0         0         0         0; */
  /*            0    1.0000         0         0         0         0   -1.0000         0         0         0; */
  /*            0         0    1.0000         0         0         0         0   -1.0000         0         0; */
  /*            0         0         0    1.0000         0         0         0         0   -1.0000         0; */
  /*            0         0         0         0    1.0000         0         0         0         0   -1.0000; */
  /*      -1.0000         0         0         0         0    1.0000         0         0         0         0; */
  /*            0   -1.0000         0         0         0         0    1.0000         0         0         0; */
  /*            0         0   -1.0000         0         0         0         0    1.0000         0         0; */
  /*            0         0         0   -1.0000         0         0         0         0    1.0000         0; */
  /*            0         0         0         0   -1.0000         0         0         0         0    1.0000]; */
  /*  Adֻ�е�5����10�и���v��ͬ����ͬ�������̶����� */
  /*  [mad,~]= size(Ad); */
  /*  �Ȱ�ǰ������ʽ�Ļ��ҵ��������� */
  /*  P=[Ad(1:3,1:3) zeros(3,10);Ad(4:mad,1:3) eye(10)]; */
  /*  P=[Ad(1:3,1:3) zeros(3,10);Ad(4:13,1:3) eye(10)]; */
  /*  ���� */
  /*  Ad_eye=P\Ad;% ���� */
  /*  �޹��е�����P_inv�ǳ����� */
  /*  P_inv=[-1     1     2     0     0     0     0     0     0     0     0     0     0; */
  /*       0    -2     0     0     0     0     0     0     0     0     0     0     0; */
  /*       1     1     2     0     0     0     0     0     0     0     0     0     0; */
  /*       1    -1    -2     1     0     0     0     0     0     0     0     0     0; */
  /*       0     2     0     0     1     0     0     0     0     0     0     0     0; */
  /*      -1    -1    -2     0     0     1     0     0     0     0     0     0     0; */
  /*       0     0     0     0     0     0     1     0     0     0     0     0     0; */
  /*       0     0     0     0     0     0     0     1     0     0     0     0     0; */
  /*      -1     1     2     0     0     0     0     0     1     0     0     0     0; */
  /*       0    -2     0     0     0     0     0     0     0     1     0     0     0; */
  /*       1     1     2     0     0     0     0     0     0     0     1     0     0; */
  /*       0     0     0     0     0     0     0     0     0     0     0     1     0; */
  /*       0     0     0     0     0     0     0     0     0     0     0     0     1]; */
  for (i = 0; i < 3; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      b_Ad[i0 + 13 * i] = Ad[i0 + 13 * i];
    }
  }

  for (i = 0; i < 10; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      b_Ad[i0 + 13 * (i + 3)] = 0.0;
    }
  }

  for (i = 0; i < 3; i++) {
    memcpy(&b_Ad[i * 13 + 3], &Ad[i * 13 + 3], 10U * sizeof(double));
  }

  for (i = 0; i < 10; i++) {
    for (i0 = 0; i0 < 10; i0++) {
      b_Ad[(i0 + 13 * (i + 3)) + 3] = iv1[i0 + 10 * i];
    }
  }

  inv_mch(b_Ad, 13.0, P_inv);

  /*  Ad_eye=[1     0     0     1     0    -1     0     0    -1     0; */
  /*          0     1     0    -1     0     0    -1     0     1     0; */
  /*          0     0     1     1     0     0     0    -1    -1     0; */
  /*          0     0     0    -1     0     0     0     0     1     0; */
  /*          0     0     0     1     0     0     0     0    -1     0; */
  /*          0     0     0    -1     0     0     0     0     1     0; */
  /*          0     0     0     1     0     0     0     0    -1     0; */
  /*          0     0     0     0     0     0     0     0     0     0; */
  /*          0     0     0     1     0     0     0     0    -1     0; */
  /*          0     0     0    -1     0     0     0     0     1     0; */
  /*          0     0     0     1     0     0     0     0    -1     0; */
  /*          0     0     0    -1     0     0     0     0     1     0; */
  /*          0     0     0     0     0     0     0     0     0     0]; */
  /*  �������Ϸ�����P_inv*Adֻ�е�5����10����v��ͬ���仯�� */
  /*  for i=1:13 */
  /*      temp1=0; */
  /*      temp2=0; */
  /*      for k=1:13 */
  /*          temp1=temp1 + P_inv(i,k)*Ad(k,5); */
  /*          temp2=temp2 + P_inv(i,k)*Ad(k,10); */
  /*      end */
  /*      Ad_eye(i,5)=temp1; */
  /*      Ad_eye(i,10)=temp2; */
  /*  end */
  /*  �����ɳڱ�����Ӧ�Ļ� */
  /*  A��Ad_eye�����䣬��5����10����P_inv*Ad�仯�Ĳ������йأ������ǳ��� */
  /*  A=[1     0     0     1     0    -1     0     0    -1     0     0     0     0     0     0     0     0     0     0     0; */
  /*     0     1     0    -1     0     0    -1     0     1     0     0     0     0     0     0     0     0     0     0     0; */
  /*     0     0     1     1     0     0     0    -1    -1     0     0     0     0     0     0     0     0     0     0     0; */
  /*     0     0     0    -1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0; */
  /*     0     0     0     1     0     0     0     0    -1     0     0     1     0     0     0     0     0     0     0     0; */
  /*     0     0     0    -1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0; */
  /*     0     0     0     1     0     0     0     0    -1     0     0     0     0     1     0     0     0     0     0     0; */
  /*     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0; */
  /*     0     0     0     1     0     0     0     0    -1     0     0     0     0     0     0     1     0     0     0     0; */
  /*     0     0     0    -1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0; */
  /*     0     0     0     1     0     0     0     0    -1     0     0     0     0     0     0     0     0     1     0     0; */
  /*     0     0     0    -1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0; */
  /*     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1]; */
  /*  for i=1:13 */
  /*      temp1=0; */
  /*      temp2=0; */
  /*      for k=1:13 */
  /*          temp1=temp1 + P_inv(i,k)*Ad5(k); */
  /*          temp2=temp2 + P_inv(i,k)*Ad10(k); */
  /*      end */
  /*      A(i,5)=temp1; */
  /*      A(i,10)=temp2; */
  /*  end */
  /*  תC��Ҫ�ر�ע���±������ */
  /*  Simplex algorithm */
  /*  Iterate through simplex algorithm main loop */
  *z = 0.0;
  for (i = 0; i < 13; i++) {
    for (i0 = 0; i0 < 10; i0++) {
      Ad_eye[i + 13 * i0] = 0.0;
      for (i1 = 0; i1 < 13; i1++) {
        Ad_eye[i + 13 * i0] += P_inv[i + 13 * i1] * Ad[i1 + 13 * i0];
      }
    }

    dv1[i] = iv2[i];
  }

  for (i = 0; i < 10; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      b_Ad_eye[i0 + 13 * i] = Ad_eye[i0 + 13 * i];
      b_Ad_eye[i0 + 13 * (i + 10)] = 0.0;
    }

    for (i0 = 0; i0 < 10; i0++) {
      b_Ad_eye[(i0 + 13 * i) + 3] = Ad_eye[(i0 + 13 * i) + 3];
      b_Ad_eye[(i0 + 13 * (i + 10)) + 3] = iv1[i0 + 10 * i];
    }
  }

  for (i = 0; i < 3; i++) {
    dv2[i] = 0.0;
  }

  dv2[7] = 20.0;
  for (i = 0; i < 4; i++) {
    dv2[i + 3] = umax[i];
    dv2[i + 8] = -umin[i];
  }

  dv2[12] = 0.0;
  memcpy(&dv3[0], &dv4[0], 20U * sizeof(double));
  Simplex_loop_C(dv1, b_Ad_eye, dv2, dv3, z, x, iters);

  /*  ���Թ滮�����η� */
  for (i = 0; i < 4; i++) {
    u[i] = x[i] - x[i + 5];
  }

  if (*z > 1.0) {
    /*  �Ŵ��˱������ٻ�ԭ����С��1�����ʾ��Ҫ��С��x�Ѿ���Ȼ����߽� */
    for (i = 0; i < 4; i++) {
      u[i] /= *z;
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void dir_alloc_sim_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void dir_alloc_sim_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for dir_alloc_sim.c
 *
 * [EOF]
 */
