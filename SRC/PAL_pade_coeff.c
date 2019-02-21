#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <complex.h>
#include "PAL_lineig.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/***********************************************************************
PAL_pade_coeff computes the coefficients of the polynomials for the Pade 
approximants of sqrt(z+1) ~ p(z)/q(z). Here we use an explicit
formula. Alternatively, the coefficients could be obtained 
with (e.g.) Maxima, through
    d=7: pade(taylor(sqrt(x+1),x,0,14),7,7);
    d=8: pade(taylor(sqrt(x+1),x,0,16),8,8);
    etc
************************************************************************
Example (Maxima):
pade(taylor(sqrt(z+1),z,0,14),7,7);
                          2          3          4         5        6       7
 16384 + 61440 z + 92160 z  + 70400 z  + 28800 z  + 6048 z  + 560 z  + 15 z
[---------------------------------------------------------------------------]
   7        6         5          4          3          2
  z  + 112 z  + 2016 z  + 13440 z  + 42240 z  + 67584 z  + 53248 z + 16384
then (normalizing by the independent term),
p[ 0] = 9.15527343750000e-04;
p[ 1] = 3.41796875000000e-02;
p[ 2] = 3.69140625000000e-01;
p[ 3] = 1.75781250000000e+00;
p[ 4] = 4.29687500000000e+00;
p[ 5] = 5.62500000000000e+00;
p[ 6] = 3.75000000000000e+00;
p[ 7] = 1.00000000000000e+00;
q[ 0] = 6.10351562500000e-05;
q[ 1] = 6.83593750000000e-03;
q[ 2] = 1.23046875000000e-01;
q[ 3] = 8.20312500000000e-01;
q[ 4] = 2.57812500000000e+00;
q[ 5] = 4.12500000000000e+00;
q[ 6] = 3.25000000000000e+00;
q[ 7] = 1.00000000000000e+00;
***********************************************************************/
int PAL_pade_coeff( int m, double *c_m, double *a_m, double *b_m )
{
    int k;
    double a, c, gamma, xi;

    c = (double)( 2*m + 1 );
    for ( k=0; k<m; k++ ) {
        a = (k+1)*M_PI / c;
        gamma = ( 2.0/c )*pow( sin( a ), 2.0 );
        xi = pow( cos( a ), 2.0 );
        a_m[k] = sqrt( gamma/xi );
        b_m[k] = xi;
    }
    *c_m = c;

    return 0;
}
