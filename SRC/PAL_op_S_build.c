#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <complex.h>
#include "PAL_lineig.h"
/*******************************************************************
PAL_op_S_build builds the operators used in the solution of A*x = b,
where  

   A  = [ K-s*M   L   ] =
        [   U   C-s*D ]

        [ I    L   ] * [ H       L     ] * [ I    0   ]
        [ 0  C-s*D ]   [ I  inv(C-s*D) ]   [ U  C-s*D ]

and
   H = ( K-s*M ) - L*inv(C-s*D)*U
     = ( K-s*M ) - sum_i h_i*L0_i*U0_i
where h_i is the (Pade) approximation for the nonlinear term 
associated with W_i = L0_i*U0_i (rank-revealing factorization). 
Then, the solution of H*x_1 = b_1, x_1=x(1:n), b_1=b(1:n), can 
be obtained with the Sherman–Morrison–Woodbury formula,
   inv( H + X*Y' )*b_1 = ( inv(H) - X*inv(S)*Y'*inv(H) )*b_1
where
   S  = I + Y'*X
   X  = inv( K-s*M )*[ h_1*L0_1 h_2*L0_2 h_3*L0_3 ... ]
   Y' = [ U0_1 U0_2 U0_3 ... ]'
   G  = C-s*D (for each nonlinear term)
*******************************************************************
- The matrices L and U associated with each nonlinear term 
  are stored in packed form, where
      L is n-by-(p*rank)
      U is (p*rank)-by-n
- The array mask (for L and U) uses 1-based indexing, therefore
  the -1 when mask is used
*******************************************************************
n_K     : (in)  dimension of K and M
n_S     : (in)  dimension of S
n_RT    : (in)  number of rational (nonlinear) terms
RT      : (in)  structure for the rational terms
factors : (in)  SuperLU structure for K-s*M = L*U
S       : (out) the operator S = ( I + Y'*X ) = L_s*U_s
ispiv   : (out) pivot indices for S
s       : (in)  shift s
*******************************************************************/
void PAL_op_S_build( int n_K, int n_S, int n_RT, PAL_RT *RT, 
                     fptr *factors, double complex *S, 
                     int *ispiv, double s )
{ 
    int i_one=1;
    int info, j, k, l, mask_j, n, n_Kxrank, p, pxp, rank;
    int *igpiv, *mask;
    double d, h, zeta;
    double one=1.0, zero=0.0;
    double *a, *b, *C, *D, *L, *ptr_X, *U, *X, *Z, *temp1, *temp2;
    double complex c_one={0.0+1.0*I}, r_one={1.0+0.0*I}, c_zero={0.0+0.0*I};
    double complex f; 
    double complex *ptr_S;

    /* Form right hand side of ( K-s*M )*X = L*h */

    X = (double *)malloc(sizeof(double)*(n_K*n_S));
    
    ptr_X = X;
    dlaset_( "A", &n_K, &n_S, &zero, &zero, X, &n_K );
    for ( l=0; l<n_RT; l++ ) {
        n = (RT+l)->n;
        p = (RT+l)->pd;
        rank = (RT+l)->rank; 
        mask = (RT+l)->mask;
        pxp = p*p;
        /* Compute h = d - a'*inv(C-s*D)*b */ 
        a = (RT+l)->a;
        b = (RT+l)->b;
        C = (RT+l)->C;
        D = (RT+l)->D;
        h = (RT+l)->d; /* initialization */
        temp1 = (double *)malloc(sizeof(double)*pxp);
        temp2 = (double *)malloc(sizeof(double)*p);
        igpiv = (int *)malloc(sizeof(int)*p);
        /* Form G = C-s*D */
        for ( k=0; k<pxp; k++ ) {
            temp1[k] = C[k] - s*D[k];
        }
        /* Solve G*x = b */
        for ( k=0; k<p; k++ ) {
            temp2[k] = b[k];
        }
        dgesv_( &p, &i_one, temp1, &p, igpiv, temp2, &p, &info ); 
        /* Dot product */  
        for ( k=0; k<p; k++ ) {     
            h -= a[k]*temp2[k];
        }
        /* Compute right-hand side, X = L*h */
        L = (RT+l)->L0;
        for ( k=0; k<rank; k++ ) {
            for ( j=0; j<n; j++ ) {
                mask_j = mask[j]-1;
                ptr_X[mask_j] = L[j]*h;
            }
            ptr_X += n_K;
            L += n;
        }
        free(temp1);
        free(temp2);
        free(igpiv);
    }

    /* Solve ( K-s*M )*X = L*h, copy X into RT */

    PAL_SuperLU_dgstrs( &n_K, &n_S, X, &n_K, &factors, &info );

    ptr_X = X;
    for ( l=0; l<n_RT; l++ ) {
        rank = (RT+l)->rank; 
        n_Kxrank = n_K*rank;
        dcopy_( &n_Kxrank, ptr_X, &i_one, (RT+l)->X, &i_one );
        ptr_X += n_Kxrank;
    }

    /* Compute S = I + Y'*X, X=inv(K-s*M)*L*h, Y'=i*sqrt(zeta)*U */

    ptr_S = S;
    zlaset_( "A", &n_S, &n_S, &c_zero, &r_one, S, &n_S );
    for ( l=0; l<n_RT; l++ ) {
        n = (RT+l)->n;
        U = (RT+l)->U0;
        zeta = (RT+l)->zeta;
        rank = (RT+l)->rank; 
        mask = (RT+l)->mask;
        ptr_X = X+(mask[0]-1); 
        f = c_one*csqrt( CMPLX(zeta,zero) );
        Z = (double *)malloc(sizeof(double)*(rank*n_S));
        dgemm_( "N", "N", &rank, &n_S, &n, &one, U, &rank, 
                ptr_X, &n_K, &zero, Z, &rank );
        for ( k=0; k<n_S; k++ ) {
            for ( j=0; j<rank; j++ ) {
                ptr_S[j+k*n_S] += f*Z[j+k*rank];
            } 
        }
        ptr_S += rank;
        free(Z);         
    }

    /* Factor S=L_S*U_S */

    zgetrf_( &n_S, &n_S, S, &n_S, ispiv, &info );

    free(X);

    return;
}
