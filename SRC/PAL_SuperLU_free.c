#include <complex.h>
#include "PAL_lineig.h"
#define HANDLE_SIZE  8
/***********************************************************
PAL_SuperLU_free is an interface to SuperLU's freeing phase.
It is based on SuperLU's FORTRAN/c_fortran_dgssv.
***********************************************************/
void PAL_SuperLU_free( fptr *factors ) 
{
    factors_t *LUfactors;
    LUfactors = (factors_t *) *factors;
    SUPERLU_FREE (LUfactors->perm_r);
    SUPERLU_FREE (LUfactors->perm_c);
    Destroy_SuperNode_Matrix(LUfactors->L);
    Destroy_CompCol_Matrix(LUfactors->U);
    SUPERLU_FREE (LUfactors->L);
    SUPERLU_FREE (LUfactors->U);
    SUPERLU_FREE (LUfactors);
    return;
}
