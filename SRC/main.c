#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <complex.h>
#include "PAL_wtime.h"
#include "PAL_lineig.h"
/*===============================================================
This program solves the nonlinear eigenvalue problem
   [ K - lambda*M +
         i*sqrt(lambda-(sigma_1)^2)*W_1 +
         i*sqrt(lambda-(sigma_2)^2)*W_2 ]*x = 0, 
where i = sqrt(-1), and K, M, W_1 and W_2 are real matrices.
Pade approximants are used for the construction of an equivalent
linear eigenvalue problem.
================================================================*/

int main()
{
    int k, nev, n_RT, status; 
    int *pd;
    double cutoff, sigma_P, *sigma_W, sigma_0;
    PAL_Matrix *K, *M, *W;
    char K_matrix[32], M_matrix[32], W1_matrix[32], W2_matrix[32];
    FILE *file;

    n_RT = 2;        /* number of nonlinear terms */
    cutoff = 1.0e-6; /* cutoff for the ranks of W_1 and W_2 */

    K = (PAL_Matrix *)malloc(sizeof(PAL_Matrix));
    M = (PAL_Matrix *)malloc(sizeof(PAL_Matrix));
    W = (PAL_Matrix *)malloc(sizeof(PAL_Matrix)*2);
    sigma_W = (double *)malloc(sizeof(double)*2);
    pd = (int *)malloc(sizeof(int)*2);

    /* The input file should contain:
       name of the file storing K in COO format
       name of the file file storing M in COO format
       name of the file file storing W_1 in COO format
       name of the file file storing W_2 in COO format
       sigma_1
       sigma_2
       sigma_P = expansion point for Pade approximants
       sigma_0 = shift of origin (currently set to sigma_P)
       pd_1 = degree of the approximant for sqrt(lambda-(sigma_1)^2)
       pd_2 = degree of the approximant for sqrt(lambda-(sigma_2)^2)
       nev = number of required eigenvalues
    */

    file = stdin;
    status = fscanf(file,"%s",&K_matrix);
    status = fscanf(file,"%s",&M_matrix);
    status = fscanf(file,"%s",&W1_matrix);
    status = fscanf(file,"%s",&W2_matrix);
    status = fscanf(file,"%lf",&sigma_W[0]);
    status = fscanf(file,"%lf",&sigma_W[1]);
    status = fscanf(file,"%lf",&sigma_P);
    status = fscanf(file,"%lf",&sigma_0);
    status = fscanf(file,"%i",&pd[0]); 
    status = fscanf(file,"%i",&pd[1]);
    status = fscanf(file,"%i",&nev);
    fclose(file);

    printf("====================================================\n");
    printf("Input parameters of the nonlinear eigenvalue problem\n");
    printf("====================================================\n");
    printf("   Matrix K  : %s\n",K_matrix);
    printf("   Matrix M  : %s\n",M_matrix);
    printf("   Matrix W1 : %s\n",W1_matrix);
    printf("   Matrix W2 : %s\n",W2_matrix);
    printf("   sigma_1 =%16.7e\n",sigma_W[0]);
    printf("   sigma_2 =%16.7e\n",sigma_W[1]);
    printf("   sigma_P =%16.7e (Pade expansion point)\n",sigma_P);
    printf("   sigma_0 =%16.7e (shift of origin)\n",sigma_0);
    printf("   pd_1 =%3i (pol. deg. 1)\n",pd[0]);
    printf("   pd_2 =%3i (pol. deg. 2)\n",pd[1]);
    printf("   nev  =%3i\n",nev);
    printf("----------------------------------------------------\n");

    /************************************************************
     Read real matrices K, M, W_1 and W_2. The assumption is that 
     the matrices are stored in coordinate format, and the first 
     line of each file contains the dimension of the matrix and 
     the number of values to be read. After the matrices are 
     read they are converted to CSR format.
     ************************************************************/

    /* read K */

    if ( PAL_read_matrix( K_matrix, K ) ) {
         printf("PAL_read_matrix failed to read matrix K\n");
         return -1;
    }

    /* read M */

    if ( PAL_read_matrix( M_matrix, M ) ){
         printf("PAL_read_matrix failed to read matrix M\n");
         return -1;
    }

    /* read W_1 */

    if ( PAL_read_matrix( W1_matrix, W+0 ) ) {
         printf("PAL_read_matrix failed to read matrix W_1\n");
         return -1;
    }

    /* read W_2 */

    if ( PAL_read_matrix( W2_matrix, W+1 ) ) {
         printf("PAL_read_matrix failed to read matrix W_1\n");
         return -1;
    }

    printf("Matrices K, M, W_1 and W_2 read successfully\n");

    /* build and solve the linearized eigenvalue problem */

    if ( PAL_lineig( K, M, W, sigma_W, sigma_P, sigma_0, 
                     cutoff, pd, n_RT, nev ) ) {
         printf("PAL_lineig failed to solve the linearized problem\n");
    }

    PAL_free_Matrix( K );
    PAL_free_Matrix( M );
    for ( k=0; k<n_RT; k++ ) PAL_free_Matrix(W+k);
    free( K );
    free( M );
    free( W );

    return 0;
}
