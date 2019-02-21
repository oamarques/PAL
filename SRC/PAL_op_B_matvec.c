#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <complex.h>
#include "PAL_lineig.h"
/******************************************************************
PAL_op_B_matvec computes y = B*x, where B, x, y are complex, and 
   B = [ M      ]
       [   -i*D ]
where M is n_M-by-n_M and stored in CSR format, D is block diagonal 
kron(eye(rank_i),D_i), and D_i is p_i-by-p_i.
******************************************************************/
void PAL_op_B_matvec( int n_M, int *ptrrow_M, int *indrow_M, 
                      double *val_M, int n_RT, PAL_RT *RT, 
                      double complex *x, double complex *y )
{
    int i, j, k, l, off_x, pd, rank;
    double *D;
    double complex zero={0.0+0.0*I};

    k = n_M;

    /* multiplication by M */

    PAL_csr_matvec( n_M, ptrrow_M, indrow_M, val_M, x, y );

    /* multiplication by D */

    off_x = n_M;
    for ( l=0; l<n_RT; l++ ) {
        pd = (RT+l)->pd;
        rank = (RT+l)->rank;
        for ( k=0; k<pd*rank; k++ ) {
            y[off_x+k] = zero;
        }
        for ( k=0; k<rank; k++ ) {
            D = (RT+l)->D;
            for ( j=0; j<pd; j++ ) {
                for ( i=0; i<pd; i++ ) {
                    y[off_x+i] += CMPLX(0.0,-D[i])*x[off_x+j]; 
                }
                D += pd;
            } 
        off_x += pd;  
        }
    }

    return;
}
