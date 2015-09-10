/*
  Other important includes that we need.
 */

/* This file was modified compared to version 5.0
*  A number of functions have new and/or extra arguments 
*  and new functions are declared.
*  
*  Modified by Ivan Ivanov  31 Jan 2007
*/

#include "index.h"
#include "blockmat.h"
#include "parameters.h"

/*  
*  Includes needed by MPI, BLACS and ScaLAPACK.
*/

#include "mpi.h"
#include "PBblacs.h"
#include "scalapack.h"
/*
  Our own routines.
  */

void triu(struct blockmatrix A);
void store_packed(struct blockmatrix A, struct blockmatrix B);
void store_unpacked(struct blockmatrix A, struct blockmatrix B);
void alloc_mat_packed(struct blockmatrix A, struct blockmatrix *pB);
void free_mat_packed(struct blockmatrix A);
int structnnz(int n, int k, struct blockmatrix C, struct constraintmatrix *constraints);
int actnnz(int n, int lda, double *A);
int bandwidth(int n, int lda, double *A);
void qreig(int n, double *maindiag, double *offdiag);
void sort_entries(int k, struct blockmatrix C, struct constraintmatrix *constraints);
double norm2(int n, double *x);
double norm1(int n, double *x);
double norminf(int n, double *x);
double Fnorm(struct blockmatrix A);
double Knorm(struct blockmatrix A);
double mat1norm(struct blockmatrix A);
double matinfnorm(struct blockmatrix A);
double calc_pobj(struct blockmatrix C, struct blockmatrix X, 
		 double constant_offset);
double calc_dobj(int k, double *a, double *y, double constant_offset);
double trace_prod(struct blockmatrix A, struct blockmatrix B);
double linesearch(int n, struct blockmatrix dX,
		  struct blockmatrix work1, struct blockmatrix work2, 
		  struct blockmatrix work3, struct blockmatrix cholinv, 
		  double *q, double *oldq, double *z, double *workvec,
		  double stepfrac,double start, int printlevel, int id);
double pinfeas(int k, struct constraintmatrix *constraints,
	       struct blockmatrix X, double *a, double *workvec,
		int *prank1);
double dinfeas(int k, struct blockmatrix C, 
	       struct constraintmatrix *constraints, double *y, 
	       struct blockmatrix Z, struct blockmatrix work1, 
		int *prank1);
double dimacserr3(int k, struct blockmatrix C, 
	       struct constraintmatrix *constraints, double *y, 
	       struct blockmatrix Z, struct blockmatrix work1,
		 int *prank1);
void op_a(int k, struct constraintmatrix *constraints,
	  struct blockmatrix X, double *result, int *prank1);
void op_at(int k, double *y, struct constraintmatrix *constraints,
	   struct blockmatrix result, int *prank1);
void makefill(int k, struct blockmatrix C, 
	      struct constraintmatrix *constraints, 
	      struct constraintmatrix *pfill, struct blockmatrix work1, 
	      int printlevel,int *prank1, int id);
void op_o(int k, struct constraintmatrix *constraints,
	  struct sparseblock **byblocks, struct blockmatrix Zi, 
          struct blockmatrix X, double *O, double *O_d, int *descOd, 
	  struct blockmatrix work1, struct blockmatrix work2,struct blockmatrix work3,
	  struct scalapackpar scapack, int *prank1);
void addscaledmat(struct blockmatrix A, double scale, struct blockmatrix B,
		  struct blockmatrix C);
void zero_mat(struct blockmatrix A);
void add_mat(struct blockmatrix A,struct blockmatrix B);
void sym_mat(struct blockmatrix A);
void make_i(struct blockmatrix A);
void copy_mat(struct blockmatrix A, struct blockmatrix B);
void mat_mult(double scale1, double scale2, struct blockmatrix A,
	      struct blockmatrix B, struct blockmatrix C);
void mat_multspa(double scale1, double scale2, struct blockmatrix A,
		 struct blockmatrix B, struct blockmatrix C, 
		 struct constraintmatrix fill);
void mat_multspb(double scale1, double scale2, struct blockmatrix A,
		 struct blockmatrix B, struct blockmatrix C, 
		 struct constraintmatrix fill);
void mat_multspc(double scale1, double scale2, struct blockmatrix A,
		 struct blockmatrix B, struct blockmatrix C, 
		 struct constraintmatrix fill);
void mat_mult_raw(int n, double scale1, double scale2, double *ap,
		  double *bp, double *cp);
void matvec(struct blockmatrix A, double *x, double *y);
void alloc_mat(struct blockmatrix A, struct blockmatrix *pB);
void free_mat(struct blockmatrix A);
void initparams(struct paramstruc *params, int *pprintlevel);
void initsoln(int n, int k, struct blockmatrix C, double *a, 
	      struct constraintmatrix *constraints, struct blockmatrix *pX0,
	      double **py0, struct blockmatrix *pZ0,int *prank1);
void trans(struct blockmatrix A);
void chol_inv(struct blockmatrix A, struct blockmatrix B,struct scalapackpar scapack);
int chol(struct blockmatrix A,struct  scalapackpar scapack);
int solvesys(int m, int ldam, int *descOd, double *O_d, double *RHS_d, 
		double *rhs, struct scalapackpar scapack);
int read_sol(char *fname, int n, int k, struct blockmatrix C, 
	     struct blockmatrix *pX, double **py, struct blockmatrix *pZ);
int read_prob(char *fname, int *pn, int *pk, struct blockmatrix *pC,
	      double **pa, struct constraintmatrix **pconstraints,
	      int printlevel, int *prank1, int id);
int write_prob(char *fname, int n, int k, struct blockmatrix C,
	       double *a, struct constraintmatrix *constraints);
int write_sol(char *fname, int n, int k, struct blockmatrix X,
	      double *y, struct blockmatrix Z);
void free_prob(int n, int k, struct blockmatrix C, double *a, 
	       struct constraintmatrix *constraints, struct blockmatrix X,
	       double *y, struct blockmatrix Z);
int sdp(int n, int k, struct blockmatrix C, double *a, double constant_offset,
	struct constraintmatrix *constraints, struct sparseblock **byblocks,
	struct constraintmatrix fill, struct blockmatrix X, double *y, 
	struct blockmatrix Z, struct blockmatrix cholxinv,
	struct blockmatrix cholzinv, double *pobj, double *dobj, 
	struct blockmatrix work1, struct blockmatrix work2,     
	struct blockmatrix work3,
	double *workvec1, double *workvec2,
	double *workvec3, double *workvec4, double *workvec5,
	double *workvec6, double *workvec7, double *workvec8, 
	double *diagO, struct blockmatrix bestx, double *besty,
	struct blockmatrix bestz, struct blockmatrix Zi, double *O,
        double *O_d, double *RHS_d, int *descOd, double *mind_tmp, 
	double *mind_gen, double *rhs, struct blockmatrix dZ,
        struct blockmatrix dX, double *dy, double *dy1, double *Fp, 
        int nfree, int *freevarblocks, int *freevarin, int printlevel,
        struct paramstruc parameters, struct scalapackpar scapack);
int easy_sdp(int n, int k, struct blockmatrix C, double *a, 
	     struct constraintmatrix *constraints, double constant_offset,
	     struct blockmatrix *pX, double **py, struct blockmatrix *pZ,
	     double *ppobj, double *pdobj, struct scalapackpar scapack,
		 struct paramstruc params, int printlevel);
void tweakgap(int n, int k, double *a, struct constraintmatrix *constraints,
	      double gap, struct blockmatrix Z, struct blockmatrix dZ, 
	      double *y, double *dy, struct blockmatrix work1, 
	      struct blockmatrix work2, struct blockmatrix work3, 
	      struct blockmatrix work4, double *workvec1, double *workvec2,
	      double *workvec3, double *workvec4, int printlevel,int *prank1,
		int id);
int bisect_(int *n, double *eps1, double *d, double *e, double *e2,
	    double *lb, double *ub, int *mm, int *m, double *w, int *ind, 
	    int *ierr, double *rv4, double *rv5);
/* my new functions */

void vec_mult_mat_raw(int n, double scale1, double scale2, double *ap, double *bp, double *cp);
double vec_mult_vec(int n, double *bp, double *cp);

/* MPI stuff declarations */

void mpiend (void);
int max (int a, int b);
void ddistr (int ictxt, int n, int numc, int nb, double *Aseq, double *Apar, int *descApar);
void dgather (int ictxt, int n, int numc, int nb,  double *Aseq, double *Apar, int *descApar);


/*
  fnlsdp declarations 
*/

/*read_data.c*/
void load_genmatrix(const char* filenameA, genmatrix *A, int nrows, int ncols, int pp);
void load_compleib(const char *code, const char *dir, genmatrix *A, genmatrix *B1, genmatrix *B, genmatrix *C1, genmatrix *C, genmatrix *D11, genmatrix *D12, genmatrix *D21, int *nx, int *nw, int *nu, int *nz, int *ny, int pp);
void load_initial_point(const char *code, const char *id, const char *dir, genmatrix *F0, genmatrix *Q0, genmatrix *V0, int nx, int nu, int ny,int pp);
void free_initial_point(genmatrix *F0, genmatrix *Q0, genmatrix *V0,int pp);
void free_compleib(genmatrix *A, genmatrix *B1, genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21,int pp);

/*print_mats.c*/
void print_genmatrix(genmatrix *A);
void fix_genmatrix(genmatrix *A);
void print_blockmatrix(struct blockmatrix *A);
void print_a(int nx, int ny, int nu, double *a);
void print_constraintmatrix(int nx, int ny, int nu, struct constraintmatrix *constraints);

/*build_mats.c*/
void copy_mat_gen(genmatrix *A,genmatrix *B,int np);
void copy_mat_gen_trans(genmatrix *A,genmatrix *B,int np);
void copy_mat_gen_inverseadd(genmatrix *A,genmatrix *B,int np);
void mat_mult_raw_gen1(double scale1,double scale2,genmatrix *ap,genmatrix *bp,genmatrix *cp,int np);
void mat_mult_raw_gen2(double scale1,double scale2,genmatrix *ap,genmatrix *bp,genmatrix *cp,int np);
void mat_mult_raw_gen3(double scale1,double scale2,genmatrix *ap,genmatrix *bp,genmatrix *cp,int np);
void build_aux_mats1(genmatrix *data,genmatrix *out, int nrows, int ncols,int np);
void build_aux_mats2(genmatrix *out, int nrows, int ncols,int np);
void make_i_gen(genmatrix *I,int nrows, int ncols,int np);
void make_diag(genmatrix *M, struct blockmatrix *diagM,int np);
void build_C(struct blockmatrix *pC, int nx, int ny, int nu, genmatrix* Q0, genmatrix *V0, genmatrix *B1, genmatrix *AF,double rho,int np);
void build_a(double **a, int nx, int ny, int nu,int np);
void build_constraints(struct constraintmatrix **pconstraints, int nx, int ny, int nz,int nu, genmatrix* Q0, genmatrix *V0, genmatrix *B, genmatrix *C,genmatrix *D12, genmatrix *AF,genmatrix * CF,int np);
void free_constraints(int nx,int ny,int nu,struct constraintmatrix *constraints,int np);
void free_mat_gen(genmatrix *A,int np);
void mats2vec(double *X, genmatrix *F0, genmatrix *Q0, genmatrix *V0,int nu, int ny, int nx);
void vec2mats(double *X, genmatrix *F0, genmatrix *Q0, genmatrix *V0,int nu, int ny, int nx);


/*solve_qp.c*/
int solve_qp(const char *code, const char *tag, double *step_x,genmatrix *F0,genmatrix *Q0, genmatrix *V0,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz,struct scalapackpar scapack, struct paramstruc parameters, int printlevel,int np,char* solutionQP);
int solve_qp_vec(const char *code, const char *tag, double *step_dx, double *X/* genmatrix *F0,genmatrix *Q0, genmatrix *V0*/,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np, char* solutionQP);

/*filter.c*/
void initialize_filter(filter *Fil, double ut, double uf);
void free_filter(filter *Fil);
int acceptable(filter *Fil,double t, double f, double beta, double gamma);
void add(filter *Fil, double t, double f);
void extract(filter *Fil, double t, double f);
void print_filter(filter *Fil);

/*objective_theta.c*/
double eval_f(genmatrix *F0,genmatrix *Q0, genmatrix *V0,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np);
double eval_f_vec(double *X,struct outputinfo *out,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np);
double trace(genmatrix *A,int np);
double eval_theta(genmatrix *F0,genmatrix *Q0, genmatrix *V0,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np);
double eval_theta_vec(double *X,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np);
double eval_nabla_f_vec(double *dX,double *X/*genmatrix *F0,genmatrix *Q0, genmatrix *V0*/,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np);
double lambda1(genmatrix *a, int n, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np);

/*restoration.c*/
void nelmin ( double fn ( double *x ,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np), int n, double *start, double *xmin, double *ynewlo, double reqmin, double *step, int konvge, int kcount, int *icount, int *numres, int *ifault ,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np);
void test_nelmin (genmatrix *F0, genmatrix *Q0, genmatrix *V0,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np);
void fminsearch (double *xmin,genmatrix *F0, genmatrix *Q0, genmatrix *V0,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np);
void fminsearch_vec(double *step_x,double *x_inicial/*genmatrix *F0, genmatrix *Q0, genmatrix *V0*/,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np);

/*algorithm.c*/
void algorithm(const char *code, const char *tag, const char *dir_compleib,const char *dir_points, int max_iter, double romax, double beta, double gamma, double sigma, double tolerancia_ro, int step_1_fail,  struct scalapackpar scapack,struct paramstruc params, int printlevel,int id);

/*outputer.c*/
void initialize_outputer(struct outputinfo *outt,int nu, int ny, int nx, int id);
//void free_outputer(struct outputinfo *outt, int id);
void print_outputer(int flag, struct outputinfo *outt, int id);

/*
  BLAS and LINPACK stuff.
  */

/*
  First, BLAS routines.
 */

#ifdef CAPSBLAS
#ifdef NOUNDERBLAS
double DNRM2();
double DASUM();
double DDOT();
int IDAMAX();
void DGEMM();
void DGEMV();
void DGER();
void DTRSM();
void DTRMV();
#else
double DNRM2_();
double DASUM_();
double DDOT_();
int IDAMAX_();
void DGEMM_();
void DGEMV_();
void DGER_();
void DTRSM_();
void DTRMV_();
#endif
#else
#ifdef NOUNDERBLAS
double dnrm2();
double dasum();
double ddot();
int idamax();
void dgemm();
void dgemv();
void dger();
void dtrsm();
void dtrmv();
#else
double dnrm2_();
double dasum_();
double ddot_();
int idamax_();
void dgemm_();
void dgemv_();
void dger_();
void dtrsm_();
void dtrmv_();
#endif
#endif

/*
  LAPACK next.
 */

#ifdef CAPSLAPACK
#ifdef NOUNDERLAPACK
void DPOTRF();
void DPOTRS();
void DPOTRI();
void DTRTRI();
#else
void DPOTRF_();
void DPOTRS_();
void DPOTRI_();
void DTRTRI_();
#endif
#else
#ifdef NOUNDERLAPACK
void dpotrf();
void dpotrs();
void dpotri();
void dtrtri();
#else
void dpotrf_();
void dpotrs_();
void dpotri_();
void dtrtri_();
#endif
#endif


