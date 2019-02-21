#include <stdio.h>
#include <malloc.h>
#include <complex.h>
#include "PAL_lineig.h"
/********************************************************************
PAL_read_matrix reads a matrix in COO format and calls PAL_coo_to_csr 
to obtain a CSR representation. It also computes the 1-norm of the
matrix.
********************************************************************/
int PAL_read_matrix( char *file, PAL_Matrix *A )
{
    int i, j, k, n, nnz, status;
    int *indcol, *indrow, *ptrrow;
    double norm=0.0, v;
    double *sum, *val;
    FILE *fp;

    fp = fopen(file, "r");
    status = fscanf(fp, "%d%d", &n, &nnz);
    if (!status) return -1;
    sum = (double *)malloc(sizeof(double)*(n));
    if (!sum) {
        fprintf(stderr, "Cannot allocate sum for %s!\n",file);
        fflush(stderr); return -1;
    }
    for ( k=0; k<n; k++ ) { sum[k] = 0.0; }
    indrow = (int *)malloc(sizeof(int)*(nnz));
    if (!indrow) {
        fprintf(stderr, "Cannot allocate indrow for %s!\n",file);
        fflush(stderr); return -1;
    }
    indcol = (int *)malloc(sizeof(int)*(nnz));
    if (!indcol) {
        fprintf(stderr, "Cannot allocate indcol for %s!\n",file);
        fflush(stderr); return -1;
    }
    val = (double *)malloc(sizeof(double)*(nnz));
    if (!val) {
        fprintf(stderr, "Cannot allocate val for %s!\n",file);
        fflush(stderr); return -1;
    }
    for ( k=0; k<nnz; k++ ) {
        status = fscanf(fp, "%d%d%lf", &i, &j, &v); 
        if (!status) return -1;
        indrow[k] = i;
        indcol[k] = j;
        val[k] = v;
        sum[j] += fabs(v);
    }
    ptrrow = (int *)malloc(sizeof(int)*(n+1));
    if (!ptrrow) {
        fprintf(stderr, "Cannot allocate ptrrow for %s!\n",file);
        fflush(stderr); return -1;
    }
    fclose(fp);

    /* Get the norm and a CSR representation */

    for ( k=0; k<n; k++ ) { 
        if ( sum[k] > norm ) { norm = sum[k]; }
    }

    PAL_coo_to_csr( n, nnz, indrow, ptrrow  );

    /* Store data in PAL_Matrix structure */

    A->norm = norm;
    A->indrow = indrow;
    A->indcol = indcol;
    A->ptrrow = ptrrow;
    A->val = val;
    A->nnz = nnz;
    A->n = n;

    free(sum);
    return 0;
}
