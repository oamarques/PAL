#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <complex.h>
#include "slu_zdefs.h"
#include "PAL_wtime.h"
#include "PAL_lineig.h"
/*=================================================================
This version of PAL_lineig solves the nonlinear eigenvalue problem
   ( K - lambda*M +
         i*sqrt(lambda-sigma_1**2 )*W_1 +
         i*sqrt(lambda-sigma_2**2 )*W_2 )*x = 0, i = sqrt(-1),
through linearization, i.e an associated generalized eigenvalue 
problem
   A*x = lambda*B*x
where
   A = [ K L ; U' C ]
   B = block_diagonal( M, D ). 
It uses Pade approximants to represent the nonlinear terms and 
exploits low-rank representations of W_1 and W_2. The generalized 
eigenvalue problem is then solved by PAL_geneig with ARPACK and 
SuperLU (i.e. shift-and-invert mode). Complex arithmetic is used
only in ARPACK.
===================================================================*/
int PAL_lineig( PAL_Matrix *K, PAL_Matrix *M, PAL_Matrix *W, 
                double *sigma_W, double sigma_P, double sigma_0, 
                double cutoff, int *pd, int n_RT, int nev )
{
    int k, ncv, n_A, n_eff_W, n_S, n_V, pdeg, rank_W;
    int *mask_eff_W;
    double d_mr, t_0, t_geneig, t_RT, t_resid, zeta;
    double *a_mr, *b_mr, *C_mr, *D_mr, *eff_W, *Q_eff_W, *L, *U;
    double complex *val, *vec;
    PAL_RT *RT_ptr, *RT;

    /****************************************
     Temporary restriction: sigma_0 = sigma_P 
     ****************************************/

    sigma_0 = sigma_P;
    
    /*****************************************************
     Form an effective (packed) representation for each 
     W, compute a QR factorization of this representation,
     and compute a Pade approximant for the corresponding
     nonlinear term. Note that different approximants can 
     implemented with a switch. 
     *****************************************************/

    RT = (PAL_RT *)malloc(sizeof(PAL_RT)*n_RT);

    /* Loop over the nonlinear (rational) terms */

    n_S = 0;
    n_V = 0;
    t_0 = PAL_wtime();
    for ( k=0; k<n_RT; k++ ) {
        pdeg = *(pd+k);
        zeta = sigma_P - pow( *(sigma_W+k), 2.0 );
        /* Form effective W */
        PAL_effective_W( W+k, &n_eff_W, &mask_eff_W, &eff_W );
        /* QR factorization of effective W */
        PAL_effective_W_QR( n_eff_W, &rank_W, eff_W, &Q_eff_W, cutoff ); 
        /* Pade approximant for the nonlinear term */
        a_mr = (double *)malloc(sizeof(double)*(pdeg+1));
        b_mr = (double *)malloc(sizeof(double)*(pdeg+1));
        C_mr = (double *)malloc(sizeof(double)*(pdeg*pdeg));
        D_mr = (double *)malloc(sizeof(double)*(pdeg*pdeg));
        if ( PAL_pade_approx( pdeg, *(sigma_W+k), sigma_P,
                              a_mr, b_mr, C_mr, D_mr, &d_mr ) ) {
           printf("PAL_lineig: PAL_pade_approx failed for term %1i\n",k+1);
           return -1;
        }
        n_S += rank_W;
        n_V += pdeg*rank_W;
        /* Kronecker product to form L and U */
        PAL_V_kron_a( n_eff_W, rank_W, pdeg, a_mr, Q_eff_W, &L );
        PAL_V_kron_b( n_eff_W, rank_W, pdeg, b_mr, eff_W, &U );
        /* Populate the RT structure */
        printf("PAL_lineig: nonlinear term %1i, rank(W) = %4i\n",k+1,rank_W);
        RT_ptr = RT+k;
        RT_ptr->X  = (double *)malloc(sizeof(double)*(rank_W*(K->n)));
        RT_ptr->L0 = (double *)malloc(sizeof(double)*(rank_W*n_eff_W));
        RT_ptr->U0 = (double *)malloc(sizeof(double)*(rank_W*n_eff_W));
        dlacpy_( "A", &n_eff_W, &rank_W, Q_eff_W, &n_eff_W, RT_ptr->L0, &n_eff_W );
        dlacpy_( "A", &rank_W, &n_eff_W, eff_W, &n_eff_W, RT_ptr->U0, &rank_W ); 
        RT_ptr->n = n_eff_W;
        RT_ptr->pd = pdeg;
        RT_ptr->rank = rank_W;
        RT_ptr->mask = mask_eff_W;
        RT_ptr->zeta = zeta;
        RT_ptr->a = a_mr;
        RT_ptr->b = b_mr;
        RT_ptr->d = d_mr;
        RT_ptr->C = C_mr;
        RT_ptr->D = D_mr;
        RT_ptr->L = L;
        RT_ptr->U = U;
        free(eff_W); 
        free(Q_eff_W);
    }
    t_RT = PAL_wtime() - t_0;
   
    t_0 = PAL_wtime();

    /*****************************************
     The following adjustment is not needed if 
     row index starts at 0. 
     *****************************************/
    for ( k=0; k<K->nnz; k++ )
        K->indrow[k] -= 1;
    for ( k=0; k<M->nnz; k++ )
        M->indrow[k] -= 1;
    for ( k=0; k<(W+0)->nnz; k++ )
        (W+0)->indrow[k] -= 1;
    for ( k=0; k<(W+1)->nnz; k++ )
        (W+1)->indrow[k] -= 1;
    
    /***************************************
     Solve the linearized eigenvalue problem 
     ***************************************/

    n_A = K->n + n_V;

    printf("PAL_lineig: non linear eigenproblem size = %5i\n",K->n);
    printf("PAL_lineig: linearized eigenproblem size = %5i\n",n_A);
    printf("----------------------------------------------------\n");
    
    ncv = MIN(nev*2+20,n_A); 
    val = (double complex*)malloc(sizeof(double complex)*(ncv)); 
    vec = (double complex*)malloc(sizeof(double complex)*(n_A*ncv));
 
    t_0 = PAL_wtime();
    PAL_geneig( n_S, n_V, n_RT, K, M, RT, nev, ncv, sigma_P, sigma_0, val, vec );
    t_geneig = PAL_wtime() - t_0;

    /* Check residuals of nonlinear eigenvalue problem */

    t_0 = PAL_wtime();
    PAL_residual( n_A, nev, sigma_W, K, M, W, val, vec ); 
    t_resid = PAL_wtime() - t_0;

    free(val); 
    free(vec);

    /* Timings */
    printf("----------------------------------------------------\n");
    printf("PAL_lineig: timings\n");
    printf("   Nonlinear (rational) terms     :%10.2e\n",t_RT);
    printf("   Generalized eigenvalue problem :%10.2e\n",t_geneig);
    printf("   Residuals of nonlinear problem :%10.2e\n",t_resid);

    for ( k=0; k<n_RT; k++ ) PAL_free_RT(RT+k);
    free(RT);
   
    return 0;
}
