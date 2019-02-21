#include "PAL_lineig.h"
/**********************************************************************
PAL_SuperLU_dgstrf is an interface to SuperLU's factorization phase. It
is based on SuperLU's FORTRAN/c_fortran_dgssv.
**********************************************************************/
void PAL_SuperLU_dgstrf( int *n, int *nnz, double *values, int *rowind, 
                         int *colptr, fptr **factors, int *info )
{
    SuperMatrix A, AC;
    SuperMatrix *L, *U;
    int panel_size, permc_spec, relax;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    superlu_options_t options;
    SuperLUStat_t stat;
    factors_t *LUfactors;
    GlobalLU_t Glu;

    double t1, t2, SuperLU_timer_();
    double   *utime;
    flops_t  *ops;

    /* Set the default input options. */
    set_default_options(&options);
    /* Symmetric mode */ 
    options.SymmetricMode = YES;
    options.ColPerm = MMD_AT_PLUS_A;
    options.DiagPivotThresh = 0.0; /* or 0.001, 0.01, etc. */

    /* Initialize the statistics variables. */
    StatInit(&stat);

    utime = (&stat)->utime;
    ops   = (&stat)->ops;

    dCreate_CompCol_Matrix( &A, *n, *n, *nnz, values, rowind, colptr,
	           	    SLU_NC, SLU_D, SLU_GE );
    L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
    U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
    if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(etree = intMalloc(*n)) ) ABORT("Malloc fails for etree[].");
    // dPrint_CompCol_Matrix("Matrix", &A );

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering 
     *   permc_spec = 1: minimum degree on structure of A'*A
     *   permc_spec = 2: minimum degree on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */    	
    permc_spec = options.ColPerm;
    get_perm_c(permc_spec, &A, perm_c);
    sp_preorder(&options, &A, perm_c, etree, &AC);
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);

    t1 = SuperLU_timer_();
    dgstrf( &options, &AC, relax, panel_size, etree, NULL, 0, 
            perm_c, perm_r, L, U, &Glu, &stat, info );
    t2 = SuperLU_timer_();
    utime[FACT] = t2 - t1;

    printf("PAL_SuperLU_dgstrf: factor time  = %11.4e\n",utime[FACT]);
    printf("PAL_SuperLU_dgstrf: factor flops = %11.4e\n",ops[FACT]);

    if ( *info == 0 ) {
       Lstore = (SCformat *) L->Store;
       Ustore = (NCformat *) U->Store;
       printf("    number of nonzeros in factor L = %d\n", Lstore->nnz);
       printf("    number of nonzeros in factor U = %d\n", Ustore->nnz);
       printf("    number of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
    } else {
       printf("SuperLU dgstrf() error returns INFO= %d\n", *info);
    }

    /* Save the LU factors in the factors handle */
    LUfactors = (factors_t *) SUPERLU_MALLOC(sizeof(factors_t));
    LUfactors->L = L;
    LUfactors->U = U;
    LUfactors->perm_c = perm_c;
    LUfactors->perm_r = perm_r;
    *factors = (fptr) LUfactors;

    SUPERLU_FREE(etree);
    Destroy_SuperMatrix_Store(&A);
    Destroy_CompCol_Permuted(&AC);
    StatFree(&stat);

    return;
}
