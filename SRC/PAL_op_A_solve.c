#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <complex.h>
#include "PAL_lineig.h"
/******************************************************************
PAL_op_A_solve solves A*x = b, where

   A  = [ K-s*M   L   ] =
        [   U   C-s*D ]

        [ I    L   ] * [ H       L     ] * [ I    0   ]
        [ 0  C-s*D ]   [ I  inv(C-s*D) ]   [ U  C-s*D ]

by computing
   x_2_1 = inv(C-s*D)*b_2
   x_1_1 = b_1 - L*x_2_1
   x_1 = inv(H)*x_1_1
   x_2 = x_2_1 - inv(C-s*D)*U*x_1
where 
   x_1 = x(  1:n), b_1 = b(  1:n)
   x_2 = x(n+1: ), b_2 = b(n+1: )
   H = (K-s*M) - L*inv(C-s*D)*U= (K-s*M) - sum_i h_i*L0_i*U0_i
and inv(H)*x_1_1 is computed through the Sherman–Morrison–Woodbury 
formula. See op_S_build. Note that (K-s*M) is real, and therefore
PAL_SuperLU_dgstrs uses real arithmetic with two right hand sides.
*******************************************************************
- n_S is the dimension of S
- The structure RT stores information for each rational term,
  i.e. the components of L, U, C and D, where n_RT is the number 
  of terms
- On output, b is overwritten by x
*******************************************************************
n_K     : (in) dimension of K and M
n_S     : (in) dimension of S
n_V     : (in) dimension of C and D
n_RT    : (in) number of rational (nonlinear) terms
RT      : (in) structure for the rational terms
factors : (in) SuperLU structure for K-s*M = L*U
S       : (in) the factors S = ( I + Y'*X ) = L_s*U_s
ispiv   : (in) pivot indices for S
b       : (in/out) right-hand side in A*x = b (overwritten by x)
*******************************************************************/
void PAL_op_A_solve( int n_K, int n_S, int n_V, int n_RT, PAL_RT *RT, 
                     fptr *factors, double complex *S, int *ispiv,
                     double s, double complex *b )
{
    int i_one=1, nrhs_1=1, nrhs_2=2;
    int info, j, k, k1, k2, l, mask_0, n, p, pxp, pxrank, 
        off_L, off_U, off_v, rank;
    int *igpiv, *mask;
    double zeta;
    double zero=0.0;
    double *C, *D;
    double complex c_one={0.0+1.0*I}, c_zero={0.0+0.0*I};
    double complex c, f;
    double *L, *U, *X, *v1, *v1_i, *v1_r;
    double complex *G, *v2;
    double complex *b_1, *b_2;

    v1 = (double *)malloc(sizeof(double)*(n_K*2));
    v2 = (double complex *)malloc(sizeof(double complex)*n_V);
    igpiv = (int *)malloc(sizeof(int)*n_V);

    b_1 = b;
    b_2 = b+n_K;
    v1_r = v1;
    v1_i = v1+n_K;

    /* x_2_1 = inv(G)*b_2, G = C-s*D */

    off_v = 0;
    zcopy_( &n_V, b_2, &i_one, v2, &i_one );
    for ( l=0; l<n_RT; l++ ) {
        p = (RT+l)->pd;
        rank = (RT+l)->rank;
        C = (RT+l)->C;
        D = (RT+l)->D;
        pxp = p*p;
        igpiv = (int *)malloc(sizeof(int)*p);
        G = (double complex *)malloc(sizeof(double complex)*pxp);
        for ( k=0; k<pxp; k++ ) { G[k] = -c_one*( C[k] - s*D[k]); }
        zgesv_( &p, &rank, G, &p, igpiv, v2+off_v, &p, &info ); 
        off_v += p*rank;
        free(igpiv);
        free(G);
    }

    /* Compute x_1_1 = x_1 - L*x_2_1 */

    k1 = 0;
    for ( k=0; k<n_K; k++ ) { 
        v1_r[k] = creal(b_1[k]); 
        v1_i[k] = cimag(b_1[k]); 
    }
    for ( l=0; l<n_RT; l++ ) {
        off_L = 0;
        n = (RT+l)->n;
        p = (RT+l)->pd;
        L = (RT+l)->L;
        rank = (RT+l)->rank;
        mask = (RT+l)->mask;
        k2 = k1 + p*rank;
        for ( k=k1; k<k2; k++ ) {
            for ( j=0; j<n; j++ ) {
                c = L[off_L+j]*v2[k];
                v1_r[mask[j]-1] -= creal(c);
                v1_i[mask[j]-1] -= cimag(c);
            }
            off_L += n;
        }
        k1 = k2;
    }

    /* Compute x_1 = inv(H)*x_1_1 using Sherman–Morrison–Woodbury */ 

    /* step1: v1 = x_1_2 = inv(K-s*M)*x_1_1 */
    PAL_SuperLU_dgstrs( &n_K, &nrhs_2, v1, &n_K, &factors, &info );
    for ( k=0; k<n_K; k++ ) {
        b_1[k] = v1_r[k] + v1_i[k]*I; 
    }
    /* step2: v2 = U*v1 */
    off_v = 0;
    zlaset_( "A", &n_S, &i_one, &c_zero, &c_zero, v2, &n_S );
    for ( l=0; l<n_RT; l++ ) {
        off_U = 0;
        n = (RT+l)->n;
        U = (RT+l)->U0;
        zeta = (RT+l)->zeta;
        rank = (RT+l)->rank;
        mask = (RT+l)->mask;
        f = c_one*csqrt( CMPLX(zeta,0.0) );
        mask_0 = mask[0]-1;
        for ( k=0; k<n; k++ ) {
            c = ( v1_r[mask_0+k] + v1_i[mask_0+k]*I )*f;
            for ( j=0; j<rank; j++ ) {
                v2[off_v+j] += U[off_U+j]*c;
            }
            off_U += rank;
        }
        off_v += rank; ;
    }
    /* step3: v2 := inv(S)*v2 */
    zgetrs_( "N", &n_S, &nrhs_1, S, &n_S, ispiv, v2, &n_S, &info );
    /* step4: v1 = X*v2 */
    off_v = 0;
    dlaset_( "A", &n_K, &i_one, &zero, &zero, v1_r, &n_K );
    dlaset_( "A", &n_K, &i_one, &zero, &zero, v1_i, &n_K );
    for ( l=0; l<n_RT; l++ ) {
        rank = (RT+l)->rank;
        X = (RT+l)->X;
        for ( k=0; k<rank; k++ ) {
            for ( j=0; j<n_K; j++ ) { v1_r[j] += X[j]*creal(v2[off_v+k]); }
            for ( j=0; j<n_K; j++ ) { v1_i[j] += X[j]*cimag(v2[off_v+k]); }
            X += n_K;
        }
        off_v += rank;
    }
    /* step5: x_1 = x_1_2 - v1 */
    for ( k=0; k<n_K; k++ ) {
        b_1[k] -= v1_r[k] + v1_i[k]*I;
    }

    /* Compute x_2 = x_2_1 - inv(C-s*D)*U*x_1 =
               inv(G)*b_2 - inv(G)*U*x_1 =
               inv(G)*( b_2 - U*x_1 )    */

    off_v = 0;
    zcopy_( &n_V, b_2, &i_one, v2, &i_one );  
    for ( l=0; l<n_RT; l++ ) {
        off_U = 0;
        n = (RT+l)->n;
        p = (RT+l)->pd;
        U = (RT+l)->U;
        zeta = (RT+l)->zeta;
        rank = (RT+l)->rank;
        mask = (RT+l)->mask; 
        f = csqrt( CMPLX(zeta,0.0) ); 
        mask_0 = mask[0]-1;
        pxrank = p*rank;
        pxp = p*p;
        /* v2 = b_2 - U*x_1 */
        for ( k=0; k<n; k++ ) {
            c = -b_1[mask_0+k]*f;
            for ( j=0; j<pxrank; j++ ) {
                v2[off_v+j] += U[off_U+j]*c;
            }
            off_U += pxrank;
        }      
        /* x_2 = inv(G)*v2 */
        C = (RT+l)->C;
        D = (RT+l)->D;
        igpiv = (int *)malloc(sizeof(int)*p);
        G = (double complex *)malloc(sizeof(double complex)*pxp);
        for ( k=0; k<pxp; k++ ) { G[k] = -c_one*( C[k] - s*D[k]); }
        zgesv_( &p, &rank, G, &p, igpiv, v2+off_v, &p, &info );
        off_v += pxrank;
        free(igpiv);
        free(G);
    }
    zcopy_( &n_V, v2, &i_one, b_2, &i_one );

    free(v1);
    free(v2);

    return;
}
