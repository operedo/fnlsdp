#include "declarations.h"

#define DEBUG_FNLSDP 1

/*
 * General matrix struct for manipulate compleib computations
 *
 * Added by Oscar Peredo 17 Jul 2009
 */

typedef struct _genmatrix {
  int nrows;
  int ncols;
  double *mat;
} genmatrix;

/*
 Functions prototypes
*/

/*read_data.c*/
void load_compleib(const char *code, const char *dir, genmatrix *A, genmatrix *B1, genmatrix *B, genmatrix *C1, genmatrix *C, genmatrix *D11, genmatrix *D12, genmatrix *D21, int *nx, int *nw, int *nu, int *nz, int *ny);
void load_initial_point(const char *code, const char *id, const char *dir, genmatrix *F0, genmatrix *Q0, genmatrix *V0, int nx, int nu, int ny);
void free_initial_point(genmatrix *F0, genmatrix *Q0, genmatrix *V0);
void free_compleib(genmatrix *A, genmatrix *B1, genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21);

/*print_mats.c*/
void print_genmatrix(genmatrix *A);
void print_blockmatrix(struct blockmatrix *A);
void print_a(int nx, int ny, int nu, double *a);
void print_constraintmatrix(int nx, int ny, int nu, struct constraintmatrix *constraints);

/*build_mats.c*/
void copy_mat_gen(genmatrix *A,genmatrix *B);
void copy_mat_gen_trans(genmatrix *A,genmatrix *B);
void mat_mult_raw_gen1(double scale1,double scale2,genmatrix *ap,genmatrix *bp,genmatrix *cp);
void mat_mult_raw_gen2(double scale1,double scale2,genmatrix *ap,genmatrix *bp,genmatrix *cp);
void mat_mult_raw_gen3(double scale1,double scale2,genmatrix *ap,genmatrix *bp,genmatrix *cp);
void build_aux_mats1(genmatrix *data,genmatrix *out, int nrows, int ncols);
void build_aux_mats2(genmatrix *out, int nrows, int ncols);
void make_i_gen(genmatrix *I,int nrows, int ncols);
void make_diag(genmatrix *M, struct blockmatrix *diagM);
void build_C(struct blockmatrix *pC, int nx, int ny, int nu, genmatrix* Q0, genmatrix *V0, genmatrix *B1, genmatrix *AF,double rho);
void build_a(double **a, int nx, int ny, int nu);
void build_constraints(struct constraintmatrix **pconstraints, int nx, int ny, int nz,int nu, genmatrix* Q0, genmatrix *V0, genmatrix *B, genmatrix *C,genmatrix *D12, genmatrix *AF,genmatrix * CF);
void free_constraints(int nx,int ny,int nu,struct constraintmatrix *constraints);
void free_mat_gen(genmatrix *A);


/*solve_qp.c*/
int solve_qp(const char *code, const char *tag, genmatrix *F0,genmatrix *Q0, genmatrix *V0,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz);

