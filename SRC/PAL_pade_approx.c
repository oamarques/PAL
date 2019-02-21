#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <complex.h>
#include "PAL_lineig.h"
/*******************************************************************
PAL_pade_approx computes the Pade approximants associated with the 
nonlinear terms (W) and the minimal realization 
    p(x)/q(x) = a'*(C-xD)^(-1)*b + d
The Pade approximant of order [m/m] for sqrt(x+1), see variable 
substituion below, can be obtained with Maxima, e.g.,
    case 7: pade(taylor(sqrt(x+1),x,0,14),7,7);
    case 8: pade(taylor(sqrt(x+1),x,0,16),8,8);
    etc
*******************************************************************
Sea also Matlab: [A,B,C,D] = tf2ss(p,q)
*******************************************************************/
int PAL_pade_approx( int m, double sigma_W, double sigma_P,
                     double *a, double *b, double *C, double *D, 
                     double *d )
{
    int i_one=1;
    int k, mxm;
    double zero=0.0, one=1.0;
    double s;
    double *p, *q;

    p = (double *)malloc(sizeof(double)*(m));
    q = (double *)malloc(sizeof(double)*(m));

    if ( PAL_pade_coeff( m, d, p, q ) ) {
       return -1;
    } 

    /* Minimal realization: p(x)/q(x) = -a'*(C-xD)^(-1)*b + d */

    mxm = m*m;
    for ( k=0; k<m; k++ ){
        a[k] = p[k];
        b[k] = p[k];
    }
    dlaset_( "A", &m, &m, &zero, &zero, C, &m );
    dlaset_( "A", &m, &m, &zero, &zero, D, &m );
    for ( k=0; k<m; k++ ) {
        C[k*(m+1)] = one;
        D[k*(m+1)] = -q[k];
    }

    /*
    Variable substitution: 
        x = (lambda - sigma_P)/(sigma_P - sigma_W^2) 
    Mapping:
        sqrt(lambda-sigma_W^2) --> a'*(C-lambda*D)^(-1)*b + d
    Scaling:
        C and D are typically very sparse and the scaling below
        could be done in a more efficient way,
    */
    
    s = sigma_P - pow( sigma_W, 2.0 );
    s = one/s;
    dscal_( &mxm, &s, D, &i_one );
    
    free(p);
    free(q);
    return 0;
}
