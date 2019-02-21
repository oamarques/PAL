#include <stdio.h>
#include <malloc.h>
#include <complex.h>
#include "PAL_lineig.h"
/*******************************************************************
PAL_V_kron_a computes V_kron = V*(kron(eye(n),a), where V is m-by-n, 
a is 1-by-p and V_kron is m-by-(p*n). V contains the Q factor from 
[Q,R] = qr(W) as computed by LAPACK's dgeqrf/dorgqr. V(1:m,:) is 
copied into the appropriate column of V_kron and then scaled.
*******************************************************************/
void PAL_V_kron_a( int m, int n, int p, double *a, double *V, 
                   double **V_kron_ptr )
{
    int i_one=1;
    int i, j, off_V, off_V_kron;
    double a_j;
    double *V_kron;

    V_kron = (double *)malloc(sizeof(double)*(m*n*p));

    off_V = 0;
    off_V_kron = 0;
    for ( i=0; i<n; i++ ) {
        for ( j=0; j<p; j++ ) {
            a_j = a[j];
            dcopy_( &m, V+off_V, &i_one, V_kron+off_V_kron, &i_one );
            dscal_( &m, &a_j, V_kron+off_V_kron, &i_one );
            off_V_kron += m;
        }
        off_V += m;
    } 

    *V_kron_ptr = V_kron;

    return;
}
