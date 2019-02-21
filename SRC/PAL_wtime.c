#ifdef unix
#include <sys/time.h>
#endif
#include <time.h>
#include "PAL_wtime.h"
/***************************************
Timing function provided by Meiyue Shao.
****************************************/
static double now;
double PAL_wtime()
{
#ifdef unix
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1.0e-6;
#else
    return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

double tic()
{
    return now = PAL_wtime();
}

double toc()
{
    return PAL_wtime()-now;
}
