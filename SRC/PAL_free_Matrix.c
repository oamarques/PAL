#include <malloc.h>
#include "PAL_lineig.h"
/*************************************************
PAL_free_Matrix frees the PAL_Matrix structure.
*************************************************/
void PAL_free_Matrix( PAL_Matrix *A )
{
    free(A->indrow);
    free(A->indcol);
    free(A->val);
    return;
}
