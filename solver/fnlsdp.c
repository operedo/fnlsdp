#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stddef.h>      /* definition of NULL */
#include <sys/time.h> 

#include "declarations.h"

double starttime;
double opotime;
double factortime;
double othertime;
double totaltime;
double othertime1;
double othertime2;
double othertime3;
struct timeval tp;
double t1;
double t2;
double endtime;

int timeval_substract (struct timeval *result,struct timeval *x,struct timeval *y);

int main(int argc, char *argv[]){
	gettimeofday(&tp, NULL);
  	starttime=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;

	int id, np, ret;
	/* dirs, files, tags ...*/
	const char* dir="../data/compleib_data/";
	const char* dirinit="../data/initial_points/";
	const char* code="REA1";
	const char* tag="01";
	char* solutionQP="solutionQP.dat-s";


	/*algorithm vars*/
	int max_iter=10000,step_1_fail=20;	
	double romax=10.0,beta=0.9,gamma=0.1,sigma=0.5,tolerancia_ro=10e-4;
	
	/* compleib data matrices */
	genmatrix A, B1, B,C1, C,D11,D12,D21;
	
	/* compleib initial point */
	genmatrix F0, Q0, V0;
	
	/* compleib matrix dimensions*/
	int nx,nw,nu,nz,ny;	
	int i,n,k;
	double rho=100.0;
	
	/* auxiliar matrices */
	genmatrix AF, CF, OUT, MF1, MF2;
	struct blockmatrix diagMF1, diagMF2;
	
	/* filter SDP */
	filter Fil;
	//double beta=0.9;
	//double gamma=0.1;

	/* sdp problem data */
	struct blockmatrix Csdp;
  	double *asdp;
  	struct constraintmatrix *constraintssdp;
  	
	/* sdp variables */
	struct blockmatrix X,Z;
  	double *y;
	
	/* sdp value of objectives functions */
  	double pobj,dobj;
	genmatrix *null=NULL;
	
	 /* Initialize the process grid */
	struct scalapackpar scapack;
  	struct paramstruc params;
  	int printlevel;

	/*load compleib data*/
	//int size=10974;
	//load_genmatrix("testeigenvalue/bcsstk17.dat",&A,size,size,0);
	//load_compleib(code, dir, &A, &B1, &B, &C1, &C, &D11, &D12, &D21, &nx, &nw, &nu, &nz, &ny, id);
	//load_initial_point(code,tag,dirinit,&F0,&Q0,&V0,nx,nu,ny, id);
	
	
	//initialize_filter(&Fil,500.0,500.0);
	
	
	MPI_Init(&argc,&argv);
  	MPI_Comm_rank (MPI_COMM_WORLD,&id);
  	MPI_Comm_size (MPI_COMM_WORLD,&np);
	scapack.id = id;
  	scapack.np = np;

  	switch (scapack.np)
        {
        	case 1: scapack.nprow=1; scapack.npcol=1;
          		break;
        	case 2: scapack.nprow=2; scapack.npcol=1;
          		break;
        	case 4: scapack.nprow=2; scapack.npcol=2;
          		break;
        	case 8: scapack.nprow=4; scapack.npcol=2;
          		break;
        	case 16: scapack.nprow=4; scapack.npcol=4;
          		break;
        	case 32: scapack.nprow=8; scapack.npcol=4;
          		break;
        	case 64: scapack.nprow=8; scapack.npcol=8;
          		break;
        	default:
          		if (scapack.id==0)
              			printf("Can not setup %d processors to a grid.\nPlease use 1,2,4,8,9,16,32 or 64 nodes to run or modify fnlsdp.c file. \n",scapack.np);
          		MPI_Finalize();
          		return(1);
	};
   	
	sl_init_(&scapack.ic,&scapack.nprow,&scapack.npcol);
	Cblacs_gridinfo(scapack.ic,&scapack.nprow,&scapack.npcol,&scapack.myrow,&scapack.mycol);
	
	/*
	if(id==0)print_filter(&Fil);
	if(acceptable(&Fil,10.0,3.0,beta,gamma)) {printf("acceptable!\n");add(&Fil,10.0,3.0);}
	if(id==0)print_filter(&Fil);
	if(acceptable(&Fil,8.0,3.1,beta,gamma)) {printf("acceptable!\n");add(&Fil,8.0,3.1);}
	if(id==0)print_filter(&Fil);
	if(acceptable(&Fil,6.0,3.2,beta,gamma)) {printf("acceptable!\n");add(&Fil,6.0,3.2);}
	if(id==0)print_filter(&Fil);
	extract(&Fil,10.0,3.0);
	if(id==0)print_filter(&Fil);
	extract(&Fil,8.0,3.1);
	if(id==0)print_filter(&Fil);
	*/	


//printf("scapack: %d,%d,%d,%d,%d\n",scapack.ic,scapack.npcol,scapack.nprow,scapack.mycol,scapack.myrow);

	
	//printf("f=%f\ntheta=%f\n",eval_f(&F0,&Q0, &V0,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id),eval_theta(&F0,&Q0, &V0,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id));
//printf("scapack: %d,%d,%d,%d,%d\n",scapack.ic,scapack.npcol,scapack.nprow,scapack.mycol,scapack.myrow);

	//test_nelmin(&F0,&Q0, &V0,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id);

	//fix_genmatrix(&Q0);
	//print_genmatrix(&V0);

	/*double ll=lambda1(&A,size,scapack,params,printlevel,id);	
	if(id==0)printf("lambda1=%f\n",ll);
	MPI_Barrier(MPI_COMM_WORLD);
	free_mat_gen(&A,0);*/

//printf("scapack: %d,%d,%d,%d,%d\n",scapack.ic,scapack.npcol,scapack.nprow,scapack.mycol,scapack.myrow);


	algorithm(code,tag,dir,dirinit,max_iter,romax,beta,gamma,sigma,tolerancia_ro,step_1_fail,scapack,params,printlevel,id);


	/*
	double *dx=(double *)calloc(nu*ny+nx*(nx+1),sizeof(double));	
	double *x_current=(double *)calloc(nu*ny+nx*(nx+1),sizeof(double));	

	printf("holaaaa\n");

	mats2vec(dx,&F0,&Q0,&V0,nu,ny,nx);
	
	for(i=0;i<nu*ny+nx*(nx+1);i++){
		printf("dx[%d]=%f\n",i,dx[i]);
	}

	mats2vec(x_current,&F0,&Q0,&V0,nu,ny,nx);
	
	for(i=0;i<nu*ny+nx*(nx+1);i++){
		printf("x_current[%d]=%f\n",i,x_current[i]);
	}
	
	printf("fobj=%f\n",eval_nabla_f_vec(dx,x_current, rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id));
	
	free(dx);
	free(x_current);
	*/

	/*solve qp*/
	//ret=0;
	//fix_genmatrix(&Q0);
	//ret=solve_qp(code,"01", &F0,&Q0, &V0,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id,solutionQP);
//printf("scapack: %d,%d,%d,%d,%d\n",scapack.ic,scapack.npcol,scapack.nprow,scapack.mycol,scapack.myrow);
	
	/*if(DEBUG_FNLSDP && id==0){
		printf("F1:\n");
		print_genmatrix(&F0);
		printf("Q1:\n");
		print_genmatrix(&Q0);
		printf("V1:\n");
		print_genmatrix(&V0);
	}*/
	
	
	//free_filter(&Fil);
	//free_initial_point(&F0,&Q0,&V0, np);
	//free_compleib(&A, &B1, &B, &C1, &C, &D11, &D12, &D21, np);
  	
	Cblacs_gridexit(scapack.ic);
	MPI_Finalize();
        gettimeofday(&tp, NULL);
  	endtime=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  	totaltime=endtime-starttime;
  	othertime=totaltime-opotime-factortime;
	if(id==0){
  		printf("Elements time: %f \n",opotime);
  		printf("Factor time: %f \n",factortime);
  		printf("Other time: %f \n",othertime);
  		printf("Total time: %f \n",totaltime);
	}
	return ret;
}

int timeval_substract (result, x, y)
	struct timeval *result, *x, *y;
{
	/* Perform the carry for the later subtraction by updating y. */
	if (x->tv_usec < y->tv_usec) {
		int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
		y->tv_usec -= 1000000 * nsec;
		y->tv_sec += nsec;
	}
	if (x->tv_usec - y->tv_usec > 1000000) {
		int nsec = (y->tv_usec - x->tv_usec) / 1000000;
		y->tv_usec += 1000000 * nsec;
		y->tv_sec -= nsec;
	}
	/* Compute the time remaining to wait.
	tv_usec  is certainly positive. */
	result->tv_sec = x->tv_sec - y->tv_sec;
	result->tv_usec = x->tv_usec - y->tv_usec;
	/* Return 1 if result is negative. */
	return x->tv_sec < y->tv_sec;
}

