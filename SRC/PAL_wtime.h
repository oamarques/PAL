#ifndef MY_WTIME_H
#define MY_WTIME_H

#define PAL_wtime PAL_wtime_
#define tic tic_
#define toc toc_

#ifdef __cplusplus
extern "C" {
#endif

double PAL_wtime();
double tic();
double toc();

#ifdef __cplusplus
}
#endif

#endif
