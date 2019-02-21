#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <complex.h>
#include "PAL_wtime.h"
#include "PAL_lineig.h"
/********************************************************
PAL_geneig solves the linearized eigenvalue problem
       A*x = lambda*B*x
where
       A = [ K L ; U' C ]
       B = block_diagonal( M, D )
by using ARPACK and SuperLU (shift-and-invert mode).
*********************************************************/
void PAL_geneig( int n_S, int n_V, int n_RT, PAL_Matrix *K, 
                 PAL_Matrix *M, PAL_RT *RT, int nev, int ncv, 
                 double sigma_P, double sigma_0, double complex *d, 
                 double complex *v )
{
    int i_one=1; 
    int ido, info, iparam[11], ipntr[14], j, 
        lworkl, k, nnz, n_A, ptrrow;
    int *ispiv, *select;
    char *bmat="I", *which="LM";
    double f, normv, tol, t_0, t_factor, t_znaupd, t_zneupd;
    double *rwork, *val_A;
    double complex sigma;
    double complex *ptr_v, *ptr_workd_0, *ptr_workd_1, *resid, *S, 
                   *workd, *workl, *workv;
    fptr *factors;

    n_A = K->n + n_V;

    val_A = (double *)malloc( sizeof(double)*(K->nnz) );
    S = (double complex *)malloc( sizeof(double complex)*(n_S*n_S) );
    ispiv = (int *)malloc( sizeof(int)*(n_S) );

    lworkl = 3*ncv*ncv+10*ncv;
    rwork  = (double *)malloc( sizeof(double)*(n_A) );
    resid  = (double complex*)malloc( sizeof(double complex)*(n_A) ); 
    workd  = (double complex*)malloc( sizeof(double complex)*(n_A*3) );
    workl  = (double complex*)malloc( sizeof(double complex)*(lworkl) ); 
    workv  = (double complex*)malloc( sizeof(double complex)*(ncv*2) );
    select = (int *)malloc( sizeof(int)*(ncv) );

    iparam[0] =  1;  // iparam(1), exact shifts
    iparam[2] = 30;  // iparam(3), maximum number of Arnoldi iterations
    iparam[6] =  1;  // iparam(7), OP = inv[A-sigma*B]*B*x

    tol  = 0.0;
    ido  = 0;
    info = 0;

    /* Factor K - sigma_0*M = (L_f)*(U_f) */

    for ( k=1; k<=K->n; k++ ) {
        ptrrow = K->ptrrow[k-1];
        nnz = K->ptrrow[k] - K->ptrrow[k-1];
        for ( j=0; j<nnz; j++ ) {
            val_A[ptrrow] = K->val[ptrrow] - sigma_0*(M->val[ptrrow]);
            ptrrow++;
        }
    }    

    t_0 = PAL_wtime();
    PAL_SuperLU_dgstrf( &K->n, &K->nnz, val_A, K->indrow, K->ptrrow, &factors, &info );
    t_factor = PAL_wtime() - t_0;

    /* Build the operators used in the solution of ( A - sigma_0*B )*x = b */

    sigma = sigma_0 - sigma_P;

    PAL_op_S_build( K->n, n_S, n_RT, RT, factors, S, ispiv, sigma );

    /* Reverse communication with ARPACK's znaupd */

    printf("PAL_geneig: znaupd, reverse communication start\n");
    t_0 = PAL_wtime();
    do {
       znaupd_( &ido, bmat, &n_A, which, &nev, &tol, resid, &ncv,
                v, &n_A, iparam, ipntr, workd, workl, &lworkl,
                rwork, &info, 1, 2 );
       switch( ido ) {
          case -1:
          case  1:
             /* 
             Compute y := OP*x = inv[A-sigma*B]*B*x
             x = workd(ipntr(1)), ipntr(1) --> ipntr[0]
             y = workd(ipntr(2)), ipntr(2) --> ipntr[1] 
             */
             ptr_workd_0 = workd+ipntr[0]-1;
             ptr_workd_1 = workd+ipntr[1]-1;
             PAL_op_B_matvec( M->n, M->ptrrow, M->indrow, M->val, n_RT, RT,
                          ptr_workd_0, ptr_workd_1 );
             PAL_op_A_solve( K->n, n_S, n_V, n_RT, RT, factors, S, ispiv,
                         sigma, ptr_workd_1 );  
             break;
          default:
             ido = 12345;
             break;
       }
    } while ( ido != 12345 );
    printf("PAL_geneig: znaupd, reverse communication end\n");
    t_znaupd = PAL_wtime() - t_0;

    /* Compute eigenvectors with ARPACK's zneupd */

    sigma = sigma_0;
    t_0 = PAL_wtime();
    if ( info < 0 ) {
       printf("PAL_geneig: znaupd exit, info  =%4i (failure)\n",info);
       printf("            1: maximum number of iterations allowed\n");
       printf("            3: no shifts could be applied in a cycle\n");
    }
    else {
       printf("PAL_geneig: znaupd exit, info  =%4i\n",info);
       zneupd_( &i_one, "A", select, d, v, &n_A, &sigma, workv, bmat, &n_A, which, 
                &nev, &tol, resid, &ncv, v, &n_A, iparam, ipntr, workd, 
                workl, &lworkl, rwork, &info, 1, 1, 2 );
       if ( info != 0 ) {
          printf("PAL_geneig: zneupd exit, info  =%4i (failure)\n",info);
          printf("            -14: no eigenvalues to sufficient accuracy\n");
          printf("            -15: inconsistent count for converged Ritz values\n");
       }
       else {
          printf("PAL_geneig: zneupd exit, nconv =%4i\n",iparam[4]);
          printf("PAL_geneig: zneupd [d]=\n");
          for ( k=0; k<iparam[4]; k++ ) {
              d[k] = 1.0/d[k] + sigma_0;
              printf("%3i%16.7e%15.7e*i%14.4e\n",k+1,
                     creal(d[k]),cimag(d[k]),fabs(creal(workl[ipntr[10]+k]))); 
          }
          /* normalize the eigenvectors (1:K->n)*/
          ptr_v = v;
          for ( k=0; k<iparam[4]; k++ ) {
              normv = dznrm2_( &K->n, ptr_v, &i_one ); f = 1.0/normv;
              for ( j=0; j<K->n; j++ ) { ptr_v[j] *= f; };
              ptr_v += n_A;
          }
       }
    }
    t_zneupd = PAL_wtime() - t_0;
    PAL_SuperLU_free( &factors );

    printf("PAL_geneig: timing breakdown\n");
    printf("   SuperLU factor :%10.2e\n",t_factor);
    printf("   ARPACK  znaupd :%10.2e\n",t_znaupd);
    printf("   ARPACK  zneupd :%10.2e\n",t_zneupd);  

    free(val_A);
    free(S);
    free(ispiv);   
    free(rwork);
    free(resid); 
    free(workd); 
    free(workv); 
    free(workl); 
    free(select); 

    return;
}
