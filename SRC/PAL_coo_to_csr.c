/**************************************************************
PAL_coo_to_csr constructs a CSR representation from the (input) 
COO format by generating pointers to the beginning of the rows,
i.e. ptrrow. ptrrow uses 0-based indexing.
**************************************************************/
void PAL_coo_to_csr( int n, int nnz, int *indrow, int *ptrrow )
{
    int j, k;

    for ( j=0; j<n; j++ ) {
        ptrrow[j] = 0;
    }
    for ( j=0; j<nnz; j++ ) {
        k = indrow[j];
        ptrrow[k] += 1;
    }
    for ( j=1; j<n+1; j++ ) {                        
        ptrrow[j] += ptrrow[j-1];
    }

    return;
}
