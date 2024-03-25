#include <stdio.h>
#include "lp_tiny.h"

int main(int argc, char *argv[]){
	{ printf("Testing n=2 infeasible\n");
		// This problem geometrically is equivalent to finding a point in
		// a null region of the plane; the region is the positive quadrant,
		// but our one constraint says that it must lie below the line
		// x + y == -1, which is impossible. Our objective function is
		// just to minimize x.
		lp_tiny lp;
		lp_tiny_status status;
		double x[2];
		lp_tiny_init(&lp, 2, 1);
		
		lp.c[0] = 1;
		lp.A[0] = 1;
		lp.A[1] = 1;
		lp.b[0] = -1;
		
		lp_tiny_solve(&lp, x, NULL, NULL, &status, NULL);
		printf(" Status: %s\n", lp_tiny_status_string(status));
		
		lp_tiny_destroy(&lp);
	}
	{ printf("Testing n=2 unbounded sublevel set\n");
		// This problem geometrically is equivalent to finding a point 
		// along y == x in the first quadrant that minimizes x. Even
		// though the origin is a solution, the sublevel sets are
		// unbounded, so this should be an unbounded problem.
		lp_tiny lp;
		lp_tiny_status status;
		double x[2];
		lp_tiny_init(&lp, 2, 1);
		
		lp.c[0] = 1;
		lp.A[0] = 1;
		lp.A[1] = -1;
		lp.b[0] = 0;
		
		lp_tiny_solve(&lp, x, NULL, NULL, &status, NULL);
		printf(" Status: %s\n", lp_tiny_status_string(status));
		
		lp_tiny_destroy(&lp);
	}
	{ printf("Testing n=2 unbounded\n");
		// This problem geometrically is equivalent to finding a point 
		// along y == x in the first quadrant that maximizes x. Not
		// only are the sublevel sets unbounded, the problem is also
		// unbounded.
		lp_tiny lp;
		lp_tiny_status status;
		double x[2];
		lp_tiny_init(&lp, 2, 1);
		
		lp.c[0] = -1;
		lp.A[0] = 1;
		lp.A[1] = -1;
		lp.b[0] = 0;
		
		lp_tiny_solve(&lp, x, NULL, NULL, &status, NULL);
		printf(" Status: %s\n", lp_tiny_status_string(status));
		
		lp_tiny_destroy(&lp);
	}
	{ printf("Testing n=2 feasibility\n");
		// This problem geometrically is equivalent to finding a point on
		// the line x+y == 1. We are only interested in feasibility,
		// so the object is just zero.
		lp_tiny lp;
		lp_tiny_status status;
		double x[2];
		lp_tiny_init(&lp, 2, 1);
		
		lp.c[0] = 0;
		lp.A[0] = 1;
		lp.A[1] = 1;
		lp.b[0] = 1;
		
		lp_tiny_solve(&lp, x, NULL, NULL, &status, NULL);
		printf(" Status: %s\n", lp_tiny_status_string(status));
		printf(" Solution: %g, %g\n", x[0], x[1]);
		
		lp_tiny_destroy(&lp);
	}
	{ printf("Testing n=2 feasibility\n");
		// This problem geometrically is equivalent to finding a point on
		// the line x+y == 1 that minimizes x + 2y.
		lp_tiny lp;
		lp_tiny_status status;
		double x[2];
		lp_tiny_init(&lp, 2, 1);
		
		lp.c[0] = 1;
		lp.c[1] = 2;
		lp.A[0] = 1;
		lp.A[1] = 1;
		lp.b[0] = 1;
		
		lp_tiny_solve(&lp, x, NULL, NULL, &status, NULL);
		printf(" Status: %s\n", lp_tiny_status_string(status));
		printf(" Solution: %g, %g\n", x[0], x[1]);
		
		lp_tiny_destroy(&lp);
	}
	{ printf("Testing n=9 feasibility\n");
		// % we use the Standard Forms for Linear Programming Problems
		// % min c'x subj. to A*x =b
		// %                  0 <= x
		// %% so we have to reformula the direction-preserving control allocation
		// % problem to:
		// % min z=[0; -1]'[u; a]   s.t.  [B -v][u; a] = 0
		// %                                umin <= u <= umax
		// %                                   0 <= a
		// % and set x=u-umin, then
		// % min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
		// %                                0 <= x <= umax-umin
		// %                                0 <= a
		// % add slack to converted inequalities into equalities of x:
		// % min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
		// %                                     x + s = umax-umin
		// %                                         0 <= x 
		// %                                         0 <= a
		// %                                         0 <= s
		// % set X=[x; a; s], that is:
		// % min z=[0; -1; 0]'[x; a; s]   s.t.  [B -v 0; I 0 I][x; a; s] = [-B*umin; umax-umin] 
		// %                                         0 <= x 
		// %                                         0 <= a
		// %                                         0 <= s
		// % A=[B -v 0; I 0 I]; b=[-B*umin; umax-umin]; c=[0; -1; 0];
		// [k,m] = size(B);
		// Aeq=[B -v zeros(k,m); eye(m) zeros(m,1) eye(m)]; 
		// beq=[-B*umin; umax-umin]; 
		// c=[zeros(m,1); -1; zeros(m,1)];
	// 	Aeq =
	// 	-0.4440         0    0.4440         0   -0.2000         0         0         0         0
    //      0   -0.4440         0    0.4440   -0.1000         0         0         0         0
    // 0.2070    0.2070    0.2070    0.2070   -0.1000         0         0         0         0
    // 1.0000         0         0         0         0    1.0000         0         0         0
    //      0    1.0000         0         0         0         0    1.0000         0         0
    //      0         0    1.0000         0         0         0         0    1.0000         0
    //      0         0         0    1.0000         0         0         0         0    1.0000


	// 	beq =
	// 			0
	// 			0
	// 		0.2890
	// 		0.6981
	// 		0.6981
	// 		0.6981
	// 		0.6981

	// 	c =
	// 		0
	// 		0
	// 		0
	// 		0
	// 		-1
	// 		0
	// 		0
	// 		0
	// 		0
		lp_tiny lp;
		lp_tiny_status status;
		double x[9];
		lp_tiny_init(&lp, 9, 7);
		
		lp.c[4] = -1;
		lp.b[2] = 0.2890;
		lp.b[3] = 0.6981;
		lp.b[4] = 0.6981;
		lp.b[5] = 0.6981;
		lp.b[6] = 0.6981;
		// i.e. A_{i,j} with n, m, is located at A[i+j*m]
		lp.A[0] = -0.4440;
		lp.A[2] = 0.2070;
		lp.A[3] = 1;

		lp.A[8] = -0.4440;
		lp.A[9] = 0.2070;
		lp.A[11] = 1;

		lp.A[14] = 0.4440;
		lp.A[16] = 0.2070;
		lp.A[19] = 1;

		lp.A[22] = 0.4440;
		lp.A[23] = 0.2070;
		lp.A[27] = 1;
		// v
		lp.A[28] = -0.2;
		lp.A[29] = -0.1;
		lp.A[30] = -0.1;

		lp.A[38] = 1;
		lp.A[46] = 1;
		lp.A[54] = 1;
		lp.A[62] = 1;



		lp_tiny_solve(&lp, x, NULL, NULL, &status, NULL);
		printf(" Status: %s\n", lp_tiny_status_string(status));
		double u[4];
		for(int i=0;i<4;i++)
		{
			u[i]=x[i]-0.3491;
		}
		if(x[4]>1)
		{
			for(int i=0;i<4;i++)
			{
				u[i]=u[i]/x[4];
			}
		}
		printf(" Solution: %g, %g, %g, %g, f: %g\n", u[0], u[1], u[2], u[3], x[4]);
		
		lp_tiny_destroy(&lp);
	}
	return 0;
}
