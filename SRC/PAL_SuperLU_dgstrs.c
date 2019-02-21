#include "PAL_lineig.h"
#define HANDLE_SIZE  8
/*****************************************************************
PAL_SuperLU_dgstrs is an interface to SuperLU's solve phase. It is
based on SuperLU's FORTRAN/c_fortran_dgssv.
*****************************************************************/
void PAL_SuperLU_dgstrs( int *n, int *nrhs, double *b, int *ldb,
                         fptr *factors, int *info )
{
    SuperMatrix B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    trans_t trans;
    SuperLUStat_t stat;
    factors_t *LUfactors;
    double SuperLU_timer_();

    trans = NOTRANS;

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* Extract the LU factors in the factors handle */
    LUfactors = (factors_t *) *factors;
    L = LUfactors->L;
    U = LUfactors->U;
    perm_c = LUfactors->perm_c;
    perm_r = LUfactors->perm_r;

    dCreate_Dense_Matrix( &B, *n, *nrhs, b, *ldb, SLU_DN, SLU_D, SLU_GE );
    
    /* Solve the system A*X=B, overwriting B with X. */
    dgstrs( trans, L, U, perm_c, perm_r, &B, &stat, info );

    Destroy_SuperMatrix_Store(&B);
    StatFree(&stat);

    return;
}
