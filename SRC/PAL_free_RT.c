#include <malloc.h>
#include "PAL_lineig.h"
/*************************************************
PAL_free_RT frees the PAL_RT structure.
*************************************************/
void PAL_free_RT( PAL_RT *RT )
{
    free(RT->mask);
    free(RT->a);
    free(RT->b);
    free(RT->C);
    free(RT->D);
    free(RT->L);
    free(RT->U);
    free(RT->L0);
    free(RT->U0);
    free(RT->X); 
    return;
}
