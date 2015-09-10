#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "declarations.h"
//#include "fnlsdp.h"

int solve_qp(const char *code, const char *tag, double *step_x, genmatrix *F0,genmatrix *Q0, genmatrix *V0,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np, char* solutionQP){
	int id, pp;
	/* dirs, files, tags ...*/
	const char* dir="../data/compleib_data/";
	const char* dirinit="../data/initial_points/";
	int i,ii,jj,counter,n,k,ret;
	genmatrix AF, CF, OUT, MF1, MF2;
	struct blockmatrix diagMF1, diagMF2;
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

	if(DEBUG_FNLSDP && np==0){printf("A:\n");
	print_genmatrix(A);}
	if(DEBUG_FNLSDP && np==0){printf("F0:\n");
	print_genmatrix(F0);}
	/*build sdp in csdp format from linearized nlsdp in compleib format*/
	build_aux_mats1(A,&AF,-1,-1,np);
	build_aux_mats1(C1,&CF,-1,-1,np);
	build_aux_mats2(&OUT,F0->nrows,C1->ncols,np);
	double scale1=1.0;
	double scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,F0,C,&OUT,np);
	if(DEBUG_FNLSDP && np==0){
		printf("OUT:\n");
		print_genmatrix(&OUT);	
	}
	scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,B,&OUT,&AF,np);
	if(DEBUG_FNLSDP && np==0){
		printf("AF:\n");
		print_genmatrix(&AF);	
	}
	if(DEBUG_FNLSDP && np==0){
	printf("D12:\n");
	print_genmatrix(D12);
	printf("OUT:\n");
	print_genmatrix(&OUT);
	printf("CF:\n");
	print_genmatrix(&CF);
	}
	scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,D12,&OUT,&CF,np);
	if(DEBUG_FNLSDP && np==0){
		printf("CF:\n");
		print_genmatrix(&CF);	
	}
	if(DEBUG_FNLSDP && np==0)printf("build_C:\n");
	build_C(&Csdp, nx, ny, nu, Q0, V0, B1, &AF,rho,np);
	if(DEBUG_FNLSDP && np==0)print_blockmatrix(&Csdp);
	build_constraints(&constraintssdp, nx, ny, nz, nu, Q0, V0, B, C, D12, &AF, &CF,np);
	if(DEBUG_FNLSDP && np==0)printf("build constraints:\n");
	if(DEBUG_FNLSDP && np==0)print_constraintmatrix(nx,ny,nu,constraintssdp);
	build_a(&asdp,nx,ny,nu,np);
	if(DEBUG_FNLSDP && np==0)print_a(nx,ny,nu,asdp);
	n=2*(nu*ny+nx*(nx+1)+1)+4*nx*nx+nx;
	k=nu*ny+nx*(nx+1)+1;
	
	initparams(&params,&printlevel);
      	initsoln(n,k,Csdp,asdp,constraintssdp,&X,&y,&Z,&params.rank1);
	write_prob("prob.dat-s",n,k,Csdp,asdp,constraintssdp);	
	ret=easy_sdp(n,k,Csdp,asdp,constraintssdp,0.0,&X,&y,&Z,&pobj,&dobj,scapack,params,printlevel);
	
    	write_sol(solutionQP,n,k,X,y,Z);
	

	for(i=1;i<=nu*ny+nx*(nx+1);i++){
		step_x[i-1]=y[i];
	}


	/*for(i=1;i<=nu*ny;i++){
		F0->mat[i-1]=y[i];
	}
	counter=1;
	for(jj=1;jj<=nx;jj++){	
		for(ii=jj;ii<=nx;ii++){			
			Q0->mat[ijtok(ii,jj,nx)]=y[counter+nu*ny];
			counter++;
		}
	}
	for(ii=1;ii<=nx;ii++){	
		for(jj=ii+1;jj<=nx;jj++){			
			Q0->mat[ijtok(ii,jj,nx)]=Q0->mat[ijtok(jj,ii,nx)];
		}
	}	
	counter=1;
	for(jj=1;jj<=nx;jj++){	
		for(ii=jj;ii<=nx;ii++){			
			V0->mat[ijtok(ii,jj,nx)]=y[counter+nu*ny+(int)(nx*(nx+1)/2)];
			counter++;
		}
	}	
	for(ii=1;ii<=nx;ii++){	
		for(jj=ii+1;jj<=nx;jj++){			
			V0->mat[ijtok(ii,jj,nx)]=V0->mat[ijtok(jj,ii,nx)];
		}
	}*/


	/*
	if(DEBUG_FNLSDP && np==0){
	printf("Z:\n");
	print_blockmatrix(&Z);
	printf("y[N+1]:\n");
	printf("%f\n",y[k]);
	}*/
	free(asdp);
	free_mat(Csdp);
	free_constraints(nx,ny,nu,constraintssdp,np);
	free(y);
	free_mat(X);
	free_mat(Z);
	free_mat_gen(&AF,np);
	free_mat_gen(&CF,np);
	free_mat_gen(&OUT,np);
	return ret;
}

int solve_qp_vec(const char *code, const char *tag, double *step_dx, double *XX/* genmatrix *F0,genmatrix *Q0, genmatrix *V0*/,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np, char* solutionQP){
	int i,j;
	
	//for(i=0;i<nu*ny+nx*(nx+1);i++)
	//	printf("X[%d]=%f\n",i,X[i]);

	genmatrix F0,Q0,V0;
	build_aux_mats2(&Q0,nx,nx,np);
	build_aux_mats2(&F0,nu,ny,np);
	build_aux_mats2(&V0,nx,nx,np);

	for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			F0.mat[i+j*nu]=XX[i+j*nu];
	}

	int count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			Q0.mat[i+j*nx]=XX[count+nu*ny];
			if(i!=j)
				Q0.mat[j+i*nx]=XX[count+nu*ny];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			V0.mat[i+j*nx]=XX[count+nu*ny+nx*(nx+1)/2];
			if(i!=j)
				V0.mat[j+i*nx]=XX[count+nu*ny+nx*(nx+1)/2];
			count++;
		}
	}	
	/*printf("F0:\n");
	print_genmatrix(&F0);
	printf("Q0:\n");
	print_genmatrix(&Q0);
	printf("V0:\n");
	print_genmatrix(&V0);
	*/

	int id, pp;
	/* dirs, files, tags ...*/
	const char* dir="../data/compleib_data/";
	const char* dirinit="../data/initial_points/";
	int ii,jj,counter,n,k,ret;
	genmatrix AF, CF, OUT, MF1, MF2;
	struct blockmatrix diagMF1, diagMF2;
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

	if(DEBUG_FNLSDP && np==0){printf("A:\n");
	print_genmatrix(A);}
	if(DEBUG_FNLSDP && np==0){printf("F0:\n");
	print_genmatrix(&F0);}
	/*build sdp in csdp format from linearized nlsdp in compleib format*/
	build_aux_mats1(A,&AF,-1,-1,np);
	build_aux_mats1(C1,&CF,-1,-1,np);
	build_aux_mats2(&OUT,F0.nrows,C1->ncols,np);
	double scale1=1.0;
	double scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&F0,C,&OUT,np);
	if(DEBUG_FNLSDP && np==0){
		printf("OUT:\n");
		print_genmatrix(&OUT);	
	}
	scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,B,&OUT,&AF,np);
	if(DEBUG_FNLSDP && np==0){
		printf("AF:\n");
		print_genmatrix(&AF);	
	}
	if(DEBUG_FNLSDP && np==0){
	printf("D12:\n");
	print_genmatrix(D12);
	printf("OUT:\n");
	print_genmatrix(&OUT);
	printf("CF:\n");
	print_genmatrix(&CF);
	}
	scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,D12,&OUT,&CF,np);
	if(DEBUG_FNLSDP && np==0){
		printf("CF:\n");
		print_genmatrix(&CF);	
	}
	if(DEBUG_FNLSDP && np==0)printf("build_C:\n");
	build_C(&Csdp, nx, ny, nu, &Q0, &V0, B1, &AF,rho,np);
	if(DEBUG_FNLSDP && np==0)print_blockmatrix(&Csdp);
	build_constraints(&constraintssdp, nx, ny, nz, nu, &Q0, &V0, B, C, D12, &AF, &CF,np);
	if(DEBUG_FNLSDP && np==0)printf("build constraints:\n");
	if(DEBUG_FNLSDP && np==0)print_constraintmatrix(nx,ny,nu,constraintssdp);
	build_a(&asdp,nx,ny,nu,np);
	if(DEBUG_FNLSDP && np==0)print_a(nx,ny,nu,asdp);
	n=2*(nu*ny+nx*(nx+1)+1)+4*nx*nx+nx;
	k=nu*ny+nx*(nx+1)+1;
	
	initparams(&params,&printlevel);
      	initsoln(n,k,Csdp,asdp,constraintssdp,&X,&y,&Z,&params.rank1);
	write_prob("prob.dat-s",n,k,Csdp,asdp,constraintssdp);	
	ret=easy_sdp(n,k,Csdp,asdp,constraintssdp,0.0,&X,&y,&Z,&pobj,&dobj,scapack,params,printlevel);
	
    	write_sol(solutionQP,n,k,X,y,Z);
	

	for(i=1;i<=nu*ny+nx*(nx+1);i++){
		step_dx[i-1]=y[i];
	}


	/*for(i=1;i<=nu*ny;i++){
		F0->mat[i-1]=y[i];
	}
	counter=1;
	for(jj=1;jj<=nx;jj++){	
		for(ii=jj;ii<=nx;ii++){			
			Q0->mat[ijtok(ii,jj,nx)]=y[counter+nu*ny];
			counter++;
		}
	}
	for(ii=1;ii<=nx;ii++){	
		for(jj=ii+1;jj<=nx;jj++){			
			Q0->mat[ijtok(ii,jj,nx)]=Q0->mat[ijtok(jj,ii,nx)];
		}
	}	
	counter=1;
	for(jj=1;jj<=nx;jj++){	
		for(ii=jj;ii<=nx;ii++){			
			V0->mat[ijtok(ii,jj,nx)]=y[counter+nu*ny+(int)(nx*(nx+1)/2)];
			counter++;
		}
	}	
	for(ii=1;ii<=nx;ii++){	
		for(jj=ii+1;jj<=nx;jj++){			
			V0->mat[ijtok(ii,jj,nx)]=V0->mat[ijtok(jj,ii,nx)];
		}
	}*/


	/*
	if(DEBUG_FNLSDP && np==0){
	printf("Z:\n");
	print_blockmatrix(&Z);
	printf("y[N+1]:\n");
	printf("%f\n",y[k]);
	}*/
	free(asdp);
	free_mat(Csdp);
	free_constraints(nx,ny,nu,constraintssdp,np);
	free(y);
	free_mat(X);
	free_mat(Z);
	free_mat_gen(&AF,np);
	free_mat_gen(&CF,np);
	free_mat_gen(&OUT,np);
	free_mat_gen(&F0,np);
	free_mat_gen(&Q0,np);
	free_mat_gen(&V0,np);
	
	return ret;
}
