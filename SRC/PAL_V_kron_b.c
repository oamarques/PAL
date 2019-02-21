#include <stdio.h>
#include <malloc.h>
#include <complex.h>
#include "PAL_lineig.h"
/******************************************************************
PAL_V_kron_b computes V_kron = kron(eye(n),b)*V, where b is p-by-1,
V is n-by-m and V_kron is (p*n)-by-m. V contains the R factor from 
[Q,R] = qr(W) as computed by LAPACK's dgeqrf/dorgqr. V(:,1:m) is  
copied into the appropriate row of V_kron and then scaled.
*******************************************************************
Note that if W is rank-deficient R is a trapezoidal matrix.
******************************************************************/
void PAL_V_kron_b( int m, int n, int p, double *b, double *V,  
                   double **V_kron_ptr )
{
    int i, j, nxp, off_V_kron;
    double b_j;
    double *V_kron;

    V_kron = (double *)malloc(sizeof(double)*(m*n*p));

    nxp = n*p;
    off_V_kron = 0;
    for ( i=0; i<n; i++ ) {
        for ( j=0; j<p; j++ ) {
            b_j = b[j];
            dcopy_( &m, V+i, &m, V_kron+off_V_kron, &nxp );
            dscal_( &m, &b_j, V_kron+off_V_kron, &nxp );
            off_V_kron += 1;
        }
    } 

    *V_kron_ptr = V_kron;

    return;
}
