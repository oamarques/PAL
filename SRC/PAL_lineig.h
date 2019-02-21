#include <complex.h>
#include "slu_ddefs.h"

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

typedef long long int fptr;
/* SuperLU structure */
typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;
/* Structure to store the matrices in the nonlinear eigenvalue problem */
typedef struct {
    int n;
    int nnz;
    int *indrow;
    int *indcol;
    int *ptrrow;
    double *val;
    double norm;
} PAL_Matrix;
/* Structure to store information about the nonlinear terms */
typedef struct {
    int n;         /* number of rows in L */
    int pd;        /* polynomial degree */
    int rank;      /* rank of W associated to the rational term */
    int *mask;     /* mask for mapping L and U into rows and columns */
    double zeta;   /* lambda-sigma^2 */
    double d;      /* d in the minimal realization a'*(C-s*D)^(-1)*b + d */
    double *a;     /* a in the minimal realization a'*(C-s*D)^(-1)*b + d */
    double *b;     /* b in the minimal realization a'*(C-s*D)^(-1)*b + d */
    double *C;     /* C in the minimal realization a'*(C-s*D)^(-1)*b + d */
    double *D;     /* D in the minimal realization a'*(C-s*D)^(-1)*b + d */
    double *L;     /* L in A = [ K-s*M  L ; U' C-s*D ] */
    double *U;     /* U in A = [ K-s*M  L ; U' C-s*D ] */
    double *L0;    /* L0 in the rank-revealing decomposition W = L0*U0 */
    double *U0;    /* U0 in the rank-revealing decomposition W = L0*U0 */
    double *X;     /* ( K-s*M )*X = L*h */
} PAL_RT;

/* BLAS */

extern void   daxpy_( int *N, double *DA, double *DX, int *INCX, double *DY, int *INCY );
extern void   dcopy_( int *N, double *DX, int *INCX, double *DY, int *INCY );
extern void   dscal_( int *N, double *DA, double *DX, int *INCX );

extern double dznrm2_( int *N, double complex *DX, int *INCX );
extern void   zaxpy_( int *N, double complex *DA, double complex *DX, int *INCX, 
                    double complex *DY, int *INCY );
extern void   zcopy_( int *N, double complex *DX, int *INCX, double complex *DY, int *INCY );
extern double complex zdotc_( int *N, double complex *DX, int *INCX, double complex *DY, int *INCY );
extern void   zscal_( int *N, double complex *DA, double complex *DX, int *INCX );

/* zgemm is defined in slu_zdefs.h */
/* zgemv is defined in slu_zdefs.h */
/* ztrsm is defined in slu_zdefs.h */

/* LAPACK */ 

extern void   dlacpy_( char *UPLO, int *M, int *N, double *A, int *LDA,
                    double *B, int *LDB );
extern double dlange_( char *NORM, int *M, int *N, double *A, int *LDA, double *WORK );
extern void   dlaset_( char *UPLO, int *M, int *N, double *ALPHA, double *BETA,
                    double *A, int *LDA );
extern void   zlacpy_( char *UPLO, int *M, int *N, double complex *A, int *LDA,
                    double complex *B, int *LDB );
extern void   zlaset_( char *UPLO, int *M, int *N, double complex *ALPHA, double complex *BETA,
                    double complex *A, int *LDA );
extern void   dgeqp3_( int *M, int *N, double *A, int *LDA, int *JPVT, double *TAU, 
                    double *WORK, int *LWORK, int *INFO );
extern void   dorgqr_( int *M, int *N, int *K, double *A, int *LDA, double *TAU, 
                    double *WORK, int *LWORK, int *INFO );
extern void   dgesv_( int *M, int *N, double *A, int *LDA, int *IPIV, double *B, 
                    int *LDB, int *INFO );
extern void   zgetrf_( int *M, int *N, double complex *A, int *LDA, int *IPIV, int *INFO );
extern void   zgetrs_( char *TRANS, int *N, int *NRHS, double complex *A, int *LDA, int *IPIV, 
                    double complex *B, int *LDB, int *INFO );

/* ARPACK */

extern void   znaupd_( int *ido, char *bmat, int *n, char *which, int *nev, double *tol, 
                    double complex *resid, int *ncv, double complex *v, int *ldv, 
                    int *iparam, int *ipntr, double complex *workd, double complex *workl, 
                    int *lworkl, double *rwork, int *info, int, int );
extern void   zneupd_( int *rvec, char *howmany, int *select, double complex *d, double complex *z, 
                    int *ldz, double complex *sigma, double complex *workv, char *bmat, int *n, 
                    char *which, int *nev, double *tol, double complex *resid, int *ncv, 
                    double complex *v, int *ldv, int *iparam, int *ipntr, 
                    double complex *workd, double complex *workl, int *lworkl, 
                    double *rwork, int *ierr, int, int, int );
/* lineig */

void PAL_coo_to_csr( int n, int nnz, int *indrow, int *ptrrow);
void PAL_csr_matvec( int n, int *ptrrow, int *indrow, double *B, 
             double complex *x, double complex *y);
void PAL_effective_W( PAL_Matrix *W, int *n_eff_W, int **mask_eff_W, double **eff_W );
int  PAL_effective_W_QR( int n_eff_W, int *rank_W, double *eff_W, 
             double **Q_eff_W, double cutoff );
void PAL_free_Matrix( PAL_Matrix *A );
void PAL_free_RT( PAL_RT *RT );
void PAL_geneig( int n_S, int n_V, int n_RT, PAL_Matrix *K, PAL_Matrix *M, 
             PAL_RT *RT, int nev, int ncv, double sigma_P, double sigma_0, 
             double complex *d, double complex *v );
int  PAL_lineig( PAL_Matrix *K, PAL_Matrix *M, PAL_Matrix *W, double *sigma_W, 
             double sigma_P, double sigma_0, double cutoff, 
             int *pd, int n_RT, int nev );
int  PAL_pade_approx( int m, double sigma, double sigma_W,
             double *a, double *b, double *C, double *D, double *d );
int  PAL_pade_coeff( int d_m, double *c_m, double *a_m, double *b_m );
void PAL_V_kron_a( int m, int n, int p, double *a, double *V, double **V_kron );
void PAL_V_kron_b( int m, int n, int p, double *b, double *V, double **V_kron );
int  PAL_read_matrix( char *file, PAL_Matrix *A );

void PAL_op_B_matvec( int n_M, int *ptrrow_M, int *indrow_M, double *val_M, 
             int n_RT, PAL_RT *RT, double complex *x, double complex *y );
void PAL_residual( int n_A, int nev, double *sigma_W, PAL_Matrix *K, 
             PAL_Matrix *M, PAL_Matrix *W, double complex *d, 
             double complex *v );

/* solver */

void PAL_op_A_solve( int n_K, int n_S, int n_V, int n_RT, PAL_RT *RT, 
             fptr *factors, double complex *S, int *ipiv,
             double s, double complex *b );
void PAL_op_S_build( int n_K, int n_S, int n_RT, PAL_RT *RT, fptr *factors,
             double complex *S, int *ispiv, double s );
void PAL_SuperLU_dgstrf( int *n, int *nnz, double *values, int *rowind, 
             int *colptr, fptr **factors, int *info );
void PAL_SuperLU_dgstrs( int *n, int *nrhs, double *b, int *ldb,
             fptr *factors, int *info );
void PAL_SuperLU_free( fptr *factors );
