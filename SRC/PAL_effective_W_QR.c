#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "PAL_lineig.h"
/***********************************************************************
PAL_effective_W_QR computes a QR factorization of the effective (packed) 
representation of the matrix W by calling LAPACK's dgeqp3 and dorgqr.
It uses a cutoff (on the diagonal of R) to determine the rank of W.
***********************************************************************/
int PAL_effective_W_QR( int n_eff_W, int *rank_W, double *eff_W, 
                        double **Q_eff_W, double cutoff )
{
    int j, k, lwork, off_eff_W, off_Q, rank, status;
    int *jpvt;
    double zero=0.0;
    double *Q, *tau, *work;

    lwork = 2*n_eff_W + (n_eff_W+1)*128;
    jpvt  = (int *)malloc(sizeof(int)*(n_eff_W));
    work  = (double *)malloc(sizeof(double)*(lwork));
    tau   = (double *)malloc(sizeof(double)*(n_eff_W));
    Q     = (double *)malloc(sizeof(double)*(n_eff_W*n_eff_W));

    /* dgeqp3 computes the QR factorization of eff_W */
    for ( k=0; k<n_eff_W; k++ ) { jpvt[k] = 0; }
    dgeqp3_( &n_eff_W, &n_eff_W, eff_W, &n_eff_W, jpvt, tau, 
             work, &lwork, &status );
    if ( status != 0 ) {
       printf("PAL_effective_W_QR: dgeqp3(eff_W), status = %6i\n",status);
       return 1;
    }
    /* check the rank of eff_W */
    rank = 0;
    for ( k=0; k<n_eff_W; k++ ) {
        if ( fabs( eff_W[k*(n_eff_W+1)] ) > cutoff )
           rank += 1;
    }
    /* copy eff_W into Q (to be used in zungqr) */
    for ( k=0; k<n_eff_W*n_eff_W; k++ ) {
        Q[k] = eff_W[k];
    }
    /* apply (the permutation) jpvt to R */
    for ( k=0; k<n_eff_W; k++ ) {
        off_Q = k*n_eff_W;
        off_eff_W = (jpvt[k]-1)*n_eff_W;
        for ( j=0; j<k+1; j++ ) {
            eff_W[off_eff_W+j] = Q[off_Q+j];
        }
        for ( j=k+1; j<n_eff_W; j++ ) {
            eff_W[off_eff_W+j] = zero;
        }      
    }
    /* dorgqr computes Q from the reflectors returned by dgeqp3 */
    dorgqr_( &n_eff_W, &n_eff_W, &rank, Q, &n_eff_W, tau, 
             work, &lwork, &status );
    if ( status != 0 ) {
       printf("PAL_effective_W_QR: dorgqr(eff_W), status = %6i\n",status);
       return 1;
    }

    *rank_W = rank;
    *Q_eff_W = Q;
    free(jpvt);
    free(work);
    free(tau); 

    return 0;
}
