/*
 * File: dir_alloc_mch.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Nov-2019 20:09:11
 */

/* Include Files */
#include <string.h>
#include "dir_alloc_mch.h"

/* Function Declarations */
static boolean_T all(const boolean_T x[20]);
static void eye(double I[100]);
static void inv_mch(double A[169], double A_inv[169]);
static void rdivide(const double x_data[], const int x_size[1], const double
                    y_data[], double z_data[], int z_size[1]);

/* Function Definitions */

/*
 * Arguments    : const boolean_T x[20]
 * Return Type  : boolean_T
 */
static boolean_T all(const boolean_T x[20])
{
  boolean_T y;
  int k;
  boolean_T exitg1;
  y = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 20)) {
    if (!x[k]) {
      y = false;
      exitg1 = true;
    } else {
      k++;
    }
  }

  return y;
}

/*
 * Arguments    : double I[100]
 * Return Type  : void
 */
static void eye(double I[100])
{
  int k;
  memset(&I[0], 0, 100U * sizeof(double));
  for (k = 0; k < 10; k++) {
    I[k + 10 * k] = 1.0;
  }
}

/*
 * �Ծ�����г����б任������
 * Arguments    : double A[169]
 *                double A_inv[169]
 * Return Type  : void
 */
static void inv_mch(double A[169], double A_inv[169])
{
  int k;
  int i;
  double div_i;
  int jj;

  /*  function Ad_eye=inv_mch(B_inv,Ad) */
  /*  Ad_eye=B_inv\Ad;% ���� */
  /*  BΪ��λ���� */
  memset(&A_inv[0], 0, 169U * sizeof(double));
  for (k = 0; k < 13; k++) {
    A_inv[k + 13 * k] = 1.0;
  }

  for (i = 0; i < 13; i++) {
    /*  ���ν��Խ��е�Ԫ�ع�һ�� */
    div_i = A[i + 13 * i];
    for (k = 0; k < 13; k++) {
      A[i + 13 * k] /= div_i;
      A_inv[i + 13 * k] /= div_i;
    }

    for (k = 0; k < 13; k++) {
      div_i = -A[k + 13 * i] / A[i + 13 * i];
      if (1 + i == 1 + k) {
        div_i = 0.0;
      }

      /*  �����б任 */
      for (jj = 0; jj < 13; jj++) {
        A[k + 13 * jj] += div_i * A[i + 13 * jj];
        A_inv[k + 13 * jj] += div_i * A_inv[i + 13 * jj];
      }
    }
  }
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
  int i1;
  z_size[0] = x_size[0];
  loop_ub = x_size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    z_data[i1] = x_data[i1] / y_data[i1];
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
 *  ������
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
 *
 *  B=[-0.5   0       0.5   0;
 *       0  -0.5    0       0.5;
 *      0.25   0.25   0.25   0.25];
 * Arguments    : const double v[3]
 *                const double umin[4]
 *                const double umax[4]
 *                double u[4]
 * Return Type  : void
 */
void dir_alloc_mch(const double v[3], const double umin[4], const double umax[4],
                   double u[4])
{
  int i;
  double b[13];
  int i0;
  double Ad[130];
  static const double dv0[12] = { -0.5, 0.0, 0.25, 0.0, -0.5, 0.25, 0.5, 0.0,
    0.25, 0.0, 0.5, 0.25 };

  static const double dv1[12] = { 0.5, -0.0, -0.25, -0.0, 0.5, -0.25, -0.5, -0.0,
    -0.25, -0.0, -0.5, -0.25 };

  static const signed char iv0[50] = { 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1 };

  static const signed char iv1[50] = { -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 };

  double b_Ad[169];
  double dv2[169];
  static const signed char iv2[100] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1 };

  double dv3[100];
  double Ad_eye[130];
  int partialTrueCount;
  double c[20];
  static const signed char iv3[20] = { 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0 };

  double A[260];
  double z;
  signed char basis[13];
  static const signed char iv4[13] = { 1, 2, 3, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 20 };

  int e;
  int exitg1;
  boolean_T b_c[20];
  boolean_T exitg2;
  int trueCount;
  double delta[13];
  int b_size[1];
  signed char tmp_data[13];
  double b_data[13];
  double A_data[13];
  double b_tmp_data[13];
  int tmp_size[1];
  double ex;
  int idx;
  signed char y_data[12];
  signed char b_y_data[12];
  double i_data[24];
  signed char c_tmp_data[24];
  double b_b_data[24];
  double b_A_data[480];

  /*  Aeq=[B -v]; */
  /*  beq=zeros(3,1); */
  /*  G=[eye(5);-eye(5)]; */
  /*  h=[umax; 20; -umin; 0]; */
  /*  ������Թ滮 */
  /*  b=[beq;h]; */
  for (i = 0; i < 3; i++) {
    b[i] = 0.0;
  }

  b[7] = 20.0;
  b[12] = 0.0;

  /*  �������Թ滮��׼�� */
  /*  Convert free variables to positively constrained variables */
  /*  Ad=[Aeq -Aeq; G -G]; */
  for (i = 0; i < 4; i++) {
    b[i + 3] = umax[i];
    b[i + 8] = -umin[i];
    for (i0 = 0; i0 < 3; i0++) {
      Ad[i0 + 13 * i] = dv0[i0 + 3 * i];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    Ad[52 + i0] = -v[i0];
  }

  for (i0 = 0; i0 < 4; i0++) {
    for (i = 0; i < 3; i++) {
      Ad[i + 13 * (i0 + 5)] = dv1[i + 3 * i0];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    Ad[117 + i0] = v[i0];
  }

  for (i0 = 0; i0 < 10; i0++) {
    for (i = 0; i < 5; i++) {
      Ad[(i + 13 * i0) + 3] = iv0[i + 5 * i0];
      Ad[(i + 13 * i0) + 8] = iv1[i + 5 * i0];
    }
  }

  /*  [mad,~]= size(Ad); */
  /*  mad=13; */
  /*  �Ȱ�ǰ������ʽ�Ļ��ҵ��������� */
  /*  B_inv=[Ad(1:3,1:3) zeros(3,mad-3);Ad(4:mad,1:3) eye(mad-3)]; */
  /*  B_inv=[Ad(1:3,1:3) zeros(3,10);Ad(4:mad,1:3) eye(10)]; */
  /*  ���� */
  /*  Ad_eye=P\Ad;% ���� */
  for (i0 = 0; i0 < 3; i0++) {
    for (i = 0; i < 3; i++) {
      b_Ad[i + 13 * i0] = Ad[i + 13 * i0];
    }
  }

  for (i0 = 0; i0 < 10; i0++) {
    for (i = 0; i < 3; i++) {
      b_Ad[i + 13 * (i0 + 3)] = 0.0;
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    memcpy(&b_Ad[i0 * 13 + 3], &Ad[i0 * 13 + 3], 10U * sizeof(double));
  }

  for (i0 = 0; i0 < 10; i0++) {
    for (i = 0; i < 10; i++) {
      b_Ad[(i + 13 * (i0 + 3)) + 3] = iv2[i + 10 * i0];
    }
  }

  inv_mch(b_Ad, dv2);
  for (i0 = 0; i0 < 13; i0++) {
    for (i = 0; i < 10; i++) {
      Ad_eye[i0 + 13 * i] = 0.0;
      for (partialTrueCount = 0; partialTrueCount < 13; partialTrueCount++) {
        Ad_eye[i0 + 13 * i] += dv2[i0 + 13 * partialTrueCount] *
          Ad[partialTrueCount + 13 * i];
      }
    }
  }

  /*  �����ɳڱ�����Ӧ�Ļ� */
  eye(dv3);
  for (i0 = 0; i0 < 10; i0++) {
    for (i = 0; i < 3; i++) {
      A[i + 13 * i0] = Ad_eye[i + 13 * i0];
      A[i + 13 * (i0 + 10)] = 0.0;
    }

    for (i = 0; i < 10; i++) {
      A[(i + 13 * i0) + 3] = Ad_eye[(i + 13 * i0) + 3];
      A[(i + 13 * (i0 + 10)) + 3] = dv3[i + 10 * i0];
    }
  }

  for (i0 = 0; i0 < 20; i0++) {
    c[i0] = iv3[i0];
  }

  for (i0 = 0; i0 < 13; i0++) {
    basis[i0] = iv4[i0];
  }

  z = 0.0;

  /*  Simplex algorithm */
  /*  Iterate through simplex algorithm main loop */
  /*  [x,z]=Simplex_loop_mch(basis, A, b, c, z); % ���Թ滮�����η� */
  /*  Initialization */
  e = -1;

  /*  [m,n] = size(A); */
  /*  m=13; */
  /*  n=20; */
  do {
    exitg1 = 0;
    for (i0 = 0; i0 < 20; i0++) {
      b_c[i0] = (c[i0] >= 0.0);
    }

    if (!all(b_c)) {
      /*  3.~isempty(c(c(N)<0)) */
      /*      e = find(c < 0, 1, 'first'); % ������������    % 4. e = N(find(c(N)<0,1)) */
      i = 0;
      exitg2 = false;
      while ((!exitg2) && (i < 20)) {
        if (c[i] < 0.0) {
          e = i;
          exitg2 = true;
        } else {
          i++;
        }
      }

      trueCount = 0;
      for (i = 0; i < 13; i++) {
        delta[i] = 1.0E+6;
        if (A[i + 13 * e] > 0.0) {
          trueCount++;
        }
      }

      partialTrueCount = 0;
      for (i = 0; i < 13; i++) {
        if (A[i + 13 * e] > 0.0) {
          tmp_data[partialTrueCount] = (signed char)(i + 1);
          partialTrueCount++;
        }
      }

      b_size[0] = trueCount;
      for (i0 = 0; i0 < trueCount; i0++) {
        b_data[i0] = b[tmp_data[i0] - 1];
      }

      for (i0 = 0; i0 < trueCount; i0++) {
        A_data[i0] = A[(tmp_data[i0] + 13 * e) - 1];
      }

      rdivide(b_data, b_size, A_data, b_tmp_data, tmp_size);
      partialTrueCount = 0;
      for (i = 0; i < 13; i++) {
        if (A[i + 13 * e] > 0.0) {
          delta[i] = b_tmp_data[partialTrueCount];
          partialTrueCount++;
        }
      }

      ex = delta[0];
      idx = 0;
      for (partialTrueCount = 0; partialTrueCount < 12; partialTrueCount++) {
        if (ex > delta[partialTrueCount + 1]) {
          ex = delta[partialTrueCount + 1];
          idx = partialTrueCount + 1;
        }
      }

      /* ѡ����� ������ */
      /*          li = basis(L);    % �����������                */
      if (delta[idx] == 1.0E+6) {
        /*  disp('System is unbounded!'); */
        exitg1 = 1;
      } else {
        /*         %% Perform pivot operation, exchanging L-row with e-coLumn variabLe */
        /*          [basis,A,b,c,z] = pivot_mch(basis,A,b,c,z,L,e); %�����������г����б任 */
        /*  Compute the coefficients of the equation for new basic variabLe x_e. */
        b[idx] /= A[idx + 13 * e];

        /*  4. */
        ex = A[idx + 13 * e];
        for (i0 = 0; i0 < 20; i0++) {
          A[idx + 13 * i0] /= ex;
        }

        /*  Compute the coefficients of the remaining constraints. */
        if (idx < 1) {
          i = 0;
        } else {
          i = idx;
          partialTrueCount = idx - 1;
          for (i0 = 0; i0 <= partialTrueCount; i0++) {
            y_data[i0] = (signed char)(1 + i0);
          }
        }

        if (13 < idx + 2) {
          trueCount = 0;
        } else {
          trueCount = 12 - idx;
          partialTrueCount = 11 - idx;
          for (i0 = 0; i0 <= partialTrueCount; i0++) {
            b_y_data[i0] = (signed char)((idx + i0) + 2);
          }
        }

        partialTrueCount = i + trueCount;
        for (i0 = 0; i0 < i; i0++) {
          i_data[i0] = y_data[i0];
        }

        for (i0 = 0; i0 < trueCount; i0++) {
          i_data[i0 + i] = b_y_data[i0];
        }

        /*   i = find(B~=li); */
        if (!(partialTrueCount == 0)) {
          for (i0 = 0; i0 < partialTrueCount; i0++) {
            c_tmp_data[i0] = (signed char)i_data[i0];
          }

          for (i0 = 0; i0 < partialTrueCount; i0++) {
            b_b_data[i0] = b[(int)i_data[i0] - 1] - A[((int)i_data[i0] + 13 * e)
              - 1] * b[idx];
          }

          for (i0 = 0; i0 < partialTrueCount; i0++) {
            b[c_tmp_data[i0] - 1] = b_b_data[i0];
          }

          for (i0 = 0; i0 < partialTrueCount; i0++) {
            b_b_data[i0] = A[((int)i_data[i0] + 13 * e) - 1];
          }

          for (i0 = 0; i0 < partialTrueCount; i0++) {
            for (i = 0; i < 20; i++) {
              b_A_data[i0 + partialTrueCount * i] = A[((int)i_data[i0] + 13 * i)
                - 1] - b_b_data[i0] * A[idx + 13 * i];
            }
          }

          for (i0 = 0; i0 < 20; i0++) {
            for (i = 0; i < partialTrueCount; i++) {
              A[((int)i_data[i] + 13 * i0) - 1] = b_A_data[i + partialTrueCount *
                i0];
            }
          }
        }

        /*  Compute the objective function */
        z -= c[e] * b[idx];
        ex = c[e];
        for (i0 = 0; i0 < 20; i0++) {
          c[i0] -= ex * A[idx + 13 * i0];
        }

        /*  Compute new sets of basic and nonbasic variabLes. */
        /*  N(find(N==e,1)) = li;  */
        basis[idx] = (signed char)(e + 1);

        /*   B(find(B==li,1)) = e; */
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  memset(&c[0], 0, 20U * sizeof(double));

  /*  [13,16]. */
  for (i0 = 0; i0 < 13; i0++) {
    c[basis[i0] - 1] = b[i0];
  }

  /*  ת���� */
  for (i0 = 0; i0 < 4; i0++) {
    u[i0] = c[i0] - c[5 + i0];
  }

  if (z > 1.0) {
    /*  �Ŵ��˱������ٻ�ԭ����С��1�����ʾ��Ҫ��С��x�Ѿ���Ȼ����߽� */
    for (i0 = 0; i0 < 4; i0++) {
      u[i0] /= z;
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void dir_alloc_mch_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void dir_alloc_mch_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for dir_alloc_mch.c
 *
 * [EOF]
 */
