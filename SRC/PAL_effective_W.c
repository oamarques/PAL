#include <stdio.h>
#include <malloc.h>
#include "PAL_lineig.h"
/************************************************************************
PAL_effective_W forms an effective (packed) representation of the matrix 
W. It is assumed that the matrix W shares a subset of the degrees of 
freedom (DOF) of K and M. To speed up the computations, PAL_effective_W 
creates a mask to map those DOF into a compact representation for W. 
Subsequently, it computes a QR factorization of this representation, and 
maps it back into the original DOF.
************************************************************************/
void PAL_effective_W( PAL_Matrix *W, int *n_eff_W, int **mask_eff_W_ptr, 
                      double **eff_W_ptr )
{
    int j, k, n_eff;
    int frst_row_W, last_row_W;
    int *mask_W, *mask_eff_W; 
    double zero=0.0;
    double *eff_W;

    /* first and last rows of W */

    frst_row_W = W->n;
    last_row_W = 0;
    for (k = 0; k < W->nnz; k++) {
        frst_row_W = MIN( frst_row_W, W->indrow[k] );
        last_row_W = MAX( last_row_W, W->indrow[k] );
    }

    /* mask uses 1-based indexing (row index) */

    n_eff = last_row_W-frst_row_W+1;
    mask_W = (int *)malloc(sizeof(int)*(n_eff));
    mask_W[0] = frst_row_W;
    for ( k=1; k<n_eff; k++ ) {
        mask_W[k] = mask_W[k-1] + 1;
    }
    eff_W = (double *)malloc(sizeof(double)*(n_eff*n_eff));
    for ( k=0; k<n_eff*n_eff; k++ ) {
        eff_W[k] = zero;
    }
    for ( k=1; k<(W->n)+1; k++ ) {
        for ( j=W->ptrrow[k-1]; j<W->ptrrow[k]; j++ ) {
            eff_W[(W->indrow[j]-frst_row_W)+
                  (W->indcol[j]-frst_row_W)*n_eff ] = W->val[j];
        }
    }

    j = 0;
    mask_eff_W = (int *)malloc(sizeof(int)*(n_eff));
    mask_eff_W[j] = frst_row_W;
    for ( k=1; k<n_eff; k++ ) {
        if ( mask_W[k] > mask_W[k-1] ) {
           j++;
           mask_eff_W[j] = frst_row_W + k;
        }   
    }

    *mask_eff_W_ptr = mask_eff_W;
    *eff_W_ptr = eff_W;
    *n_eff_W = n_eff;
    free(mask_W);

    return;
}
