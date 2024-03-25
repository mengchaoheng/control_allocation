#include <stdio.h>
#include <sys/time.h>
#include "LinAlg.h"

double counter()
{
    struct timeval now;
    alglib_impl::ae_int64_t r, v;
    gettimeofday(&now, NULL);
    v = now.tv_sec;
    r = v*1000;
    v = now.tv_usec/1000;
    r = r+v;
    return 0.001*r;
}

int main()
{
    alglib::real_2d_array a, b, c;
    int n = 2000;
    int i, j;
    double timeneeded, flops;
    
    // Initialize arrays
    a.setlength(n, n);
    b.setlength(n, n);
    c.setlength(n, n);
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
        {
            a[i][j] = alglib::randomreal()-0.5;
            b[i][j] = alglib::randomreal()-0.5;
            c[i][j] = 0.0;
        }
    
    // Set global threading settings (applied to all ALGLIB functions);
    // default is to perform serial computations, unless parallel execution
    // is activated. Parallel execution tries to utilize all cores; this
    // behavior can be changed with alglib::setnworkers() call.
    alglib::setglobalthreading(alglib::parallel);
    
    // Perform matrix-matrix product.
    flops = 2*pow((double)n, (double)3);
    timeneeded = counter();
    alglib::rmatrixgemm(
        n, n, n,
        1.0,
        a, 0, 0, 0,
        b, 0, 0, 1,
        0.0,
        c, 0, 0);
    timeneeded = counter()-timeneeded;
    
    // Evaluate performance
    printf("Performance is %.1f GFLOPS\n", (double)(1.0E-9*flops/timeneeded));
    
    return 0;
}
