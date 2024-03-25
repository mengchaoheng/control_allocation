#include <stdio.h>
#include "lp_tiny.h"
#include <stdlib.h>
#include <math.h>
double Z, X, Y, phi, theta;
void generate_uniform_sphere_points(int num_points) {
    for (int i = 0; i < num_points; ++i) {
        
        phi = acos(2.0 * rand() / RAND_MAX - 1); // 位于[0, pi]
        theta = 2.0 * M_PI * rand() / RAND_MAX; // 位于[0, 2pi)
        X = cos(phi) * sin(theta);
        Y = sin(phi);
        Z = cos(phi) * cos(theta);
 
        printf("Point %d: (%f, %f, %f)\n", i, X, Y, Z);
    }
}
int main(int argc, char *argv[]){

    srand(123); // 用固定的种子以保证可重复性
    int num=1000;
    generate_uniform_sphere_points(num); // 生成10个点

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
    // -v
    lp.A[28] = 0;
    lp.A[29] = 0;
    lp.A[30] = 0;

    lp.A[38] = 1;
    lp.A[46] = 1;
    lp.A[54] = 1;
    lp.A[62] = 1;
    printf("sphere test\n");
	for(int i=0;i<num;i++)
	{
        // -v
		lp.A[28] = -X;
		lp.A[29] = -Y;
		lp.A[30] = -Z;
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
	}
    lp_tiny_destroy(&lp);
	return 0;
}
