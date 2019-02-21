#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <complex.h>
#include "PAL_lineig.h"
/**************************************************************
Given an approximate eigenpair (sigma,x), PAL_residual computes
      resd = [ K - lambda*M + 
               i*sqrt(lambda-sigma_1^2)*W_1 + 
               i*sqrt(lambda-sigma_2^2)*W_2 ]*x
and resd/f, where
      f = norm(K,1) + 
          abs(lambda)*norm(M,1) +
          sqrt(abs(lambda-sigma_1^2))*norm(W_1,1) +
          sqrt(abs(lambda-sigma_2^2))*norm(W_2,1)
**************************************************************/
void PAL_residual( int n_A, int nev, double *sigma_W, 
                   PAL_Matrix *K, PAL_Matrix *M, PAL_Matrix *W, 
                   double complex *d, double complex *v ) 
{
     int i_one=1;
     int k, n;
     double f, normr;
     double complex s;
     double complex *ptr_v, *resd, *work;

     n = K->n;
     ptr_v = v;

     resd = (double complex*)malloc(sizeof(double complex)*n); 
     work = (double complex*)malloc(sizeof(double complex)*n); 

     printf("PAL_residual: eig=sqrt(d) and residuals\n");
     printf("        real(eig)      imag(eig)  ");
     printf("        ||r||        ||r||/f\n");

     for ( k=0; k<nev; k++ ) {
         /* resd = K*v */
         PAL_csr_matvec( K->n, K->ptrrow, K->indrow, K->val, ptr_v, resd );
         /* resd = resd - lambda*M*v */
         s = -d[k];
         PAL_csr_matvec( M->n, M->ptrrow, M->indrow, M->val, ptr_v, work );
         zaxpy_( &n, &s, work, &i_one, resd, &i_one );
         /* resd = resd + i*sqrt(lambda - sigma_1^2)*W_1*v */
         s = I*csqrt( d[k] - pow( *(sigma_W+0), 2.0 ) );
         PAL_csr_matvec( (W+0)->n, (W+0)->ptrrow, (W+0)->indrow, (W+0)->val, ptr_v, work );
         zaxpy_( &n, &s, work, &i_one, resd, &i_one );
         /* resd = resd + i*sqrt(lambda - sigma_2^2)*W_2*v */
         s = I*csqrt( d[k] - pow( *(sigma_W+1), 2.0 ) );
         PAL_csr_matvec( (W+1)->n, (W+1)->ptrrow, (W+1)->indrow, (W+1)->val, ptr_v, work );
         zaxpy_( &n, &s, work, &i_one, resd, &i_one );
         /* norm(resd) and norm(resd)/f */
         f = K->norm + cabs(d[k])*M->norm +
             sqrt(cabs(d[k]-pow(*(sigma_W+0),2.0)))*((W+0)->norm) +
             sqrt(cabs(d[k]-pow(*(sigma_W+1),2.0)))*((W+0)->norm);
         normr = dznrm2_( &n, resd, &i_one );
         printf("%3i%16.7e%15.7e*i%14.4e%13.4e\n",
                k+1,creal(csqrt(d[k])),cimag(csqrt(d[k])),normr,normr/f); 
         ptr_v += n_A;
     }

     printf("norm(K)   =%16.7e\n",K->norm);
     printf("norm(M)   =%16.7e\n",M->norm);
     printf("norm(W_1) =%16.7e\n",(W+0)->norm);
     printf("norm(W_2) =%16.7e\n",(W+1)->norm);

     free(resd); 
     free(work);
     return;
}
