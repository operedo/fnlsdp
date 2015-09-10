#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "declarations.h"

#define max(a,b) ((a)<(b)?(b):(a))

double eval_f(genmatrix *F0,genmatrix *Q0, genmatrix *V0,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np){
	double tt;
	genmatrix CF, CFT, OUT, MF1, MF2, RES1, RES2;
	build_aux_mats1(C1,&CF,-1,-1,np);
	build_aux_mats2(&CFT,CF.ncols,CF.nrows,np);
	build_aux_mats2(&OUT,F0->nrows,C1->ncols,np);
	build_aux_mats2(&RES1,Q0->nrows,C1->nrows,np);
	build_aux_mats2(&RES2,C1->nrows,C1->nrows,np);
	double scale1=1.0;
	double scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,F0,C,&OUT,np);
	if(DEBUG_FNLSDP && np==0){
		printf("OUT:\n");
		print_genmatrix(&OUT);	
	}
	scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,D12,&OUT,&CF,np);
	if(DEBUG_FNLSDP && np==0){
		printf("CF:\n");
		print_genmatrix(&CF);	
	}
	
	copy_mat_gen_trans(&CF,&CFT,np);
	scale1=1.0;
	scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,Q0,&CFT,&RES1,np);
	mat_mult_raw_gen1(scale1,scale2,&CF,&RES1,&RES2,np);
	if(DEBUG_FNLSDP && np==0){
		printf("RES2:\n");
		print_genmatrix(&RES2);	
	}
	
	tt=trace(&RES2,np);

	free_mat_gen(&CF,np);
	free_mat_gen(&CFT,np);
	free_mat_gen(&OUT,np);
	free_mat_gen(&RES1,np);
	free_mat_gen(&RES2,np);
	return tt;
}

double eval_f_vec(double *X/*genmatrix *F0,genmatrix *Q0, genmatrix *V0*/,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np){
	int i,j;
	
	//for(i=0;i<nu*ny+nx*(nx+1);i++)
	//	printf("X[%d]=%f\n",i,X[i]);

	genmatrix F0,Q0,V0;
	build_aux_mats2(&Q0,nx,nx,np);
	build_aux_mats2(&F0,nu,ny,np);
	build_aux_mats2(&V0,nx,nx,np);

	for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			F0.mat[i+j*nu]=X[i+j*nu];
	}

	int count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			Q0.mat[i+j*nx]=X[count+nu*ny];
			if(i!=j)
				Q0.mat[j+i*nx]=X[count+nu*ny];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			V0.mat[i+j*nx]=X[count+nu*ny+nx*(nx+1)/2];
			if(i!=j)
				V0.mat[j+i*nx]=X[count+nu*ny+nx*(nx+1)/2];
			count++;
		}
	}	

	double tt;
	genmatrix CF, CFT, OUT, MF1, MF2, RES1, RES2;
	build_aux_mats1(C1,&CF,-1,-1,np);
	build_aux_mats2(&CFT,CF.ncols,CF.nrows,np);
	build_aux_mats2(&OUT,F0.nrows,C1->ncols,np);
	build_aux_mats2(&RES1,Q0.nrows,C1->nrows,np);
	build_aux_mats2(&RES2,C1->nrows,C1->nrows,np);
	double scale1=1.0;
	double scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&F0,C,&OUT,np);
	if(DEBUG_FNLSDP && np==0){
		printf("OUT:\n");
		print_genmatrix(&OUT);	
	}
	scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,D12,&OUT,&CF,np);
	if(DEBUG_FNLSDP && np==0){
		printf("CF:\n");
		print_genmatrix(&CF);	
	}
	
	copy_mat_gen_trans(&CF,&CFT,np);
	scale1=1.0;
	scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&Q0,&CFT,&RES1,np);
	mat_mult_raw_gen1(scale1,scale2,&CF,&RES1,&RES2,np);
	if(DEBUG_FNLSDP && np==0){
		printf("RES2:\n");
		print_genmatrix(&RES2);	
	}
	
	tt=trace(&RES2,np);

	outt->filter_cand_f=tt;

	free_mat_gen(&CF,np);
	free_mat_gen(&CFT,np);
	free_mat_gen(&OUT,np);
	free_mat_gen(&RES1,np);
	free_mat_gen(&RES2,np);
	
	free_mat_gen(&F0,np);
	free_mat_gen(&V0,np);
	free_mat_gen(&Q0,np);
	
	return tt;
}

double eval_nabla_f_vec(double *dX,double *X/*genmatrix *F0,genmatrix *Q0, genmatrix *V0*/,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np){
	int i,j;
	
	//for(i=0;i<nu*ny+nx*(nx+1);i++)
	//	printf("X[%d]=%f\n",i,X[i]);

	genmatrix F0,Q0,V0;
	build_aux_mats2(&Q0,nx,nx,np);
	build_aux_mats2(&F0,nu,ny,np);
	build_aux_mats2(&V0,nx,nx,np);

	for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			F0.mat[i+j*nu]=X[i+j*nu];
	}

	int count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			Q0.mat[i+j*nx]=X[count+nu*ny];
			if(i!=j)
				Q0.mat[j+i*nx]=X[count+nu*ny];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			V0.mat[i+j*nx]=X[count+nu*ny+nx*(nx+1)/2];
			if(i!=j)
				V0.mat[j+i*nx]=X[count+nu*ny+nx*(nx+1)/2];
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

	genmatrix dF,dQ,dV;
	build_aux_mats2(&dQ,nx,nx,np);
	build_aux_mats2(&dF,nu,ny,np);
	build_aux_mats2(&dV,nx,nx,np);

	for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			dF.mat[i+j*nu]=dX[i+j*nu];
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			dQ.mat[i+j*nx]=dX[count+nu*ny];
			if(i!=j)
				dQ.mat[j+i*nx]=dX[count+nu*ny];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			dV.mat[i+j*nx]=dX[count+nu*ny+nx*(nx+1)/2];
			if(i!=j)
				dV.mat[j+i*nx]=dX[count+nu*ny+nx*(nx+1)/2];
			count++;
		}
	}

	/*printf("dF:\n");
	print_genmatrix(&dF);
	printf("dQ:\n");
	print_genmatrix(&dQ);
	printf("dV:\n");
	print_genmatrix(&dV);
	*/


	double tt;
	genmatrix CF, CFT, dCF, dCFT, dFT, dQT, dVT, OUT, OUT2, OUT3, MF1, MF2, MF3, MF4, RES1, RES2, RES3, RES4;
	build_aux_mats1(C1,&CF,-1,-1,np);
	build_aux_mats2(&CFT,CF.ncols,CF.nrows,np);
	build_aux_mats2(&dCF,D12->nrows,C->ncols,np);
	build_aux_mats2(&dCFT,dCF.ncols,dCF.nrows,np);

	build_aux_mats2(&OUT,F0.nrows,C1->ncols,np);
	build_aux_mats2(&OUT2,Q0.nrows,C1->nrows,np);
	build_aux_mats2(&OUT3,Q0.nrows,D12->nrows,np);
	build_aux_mats2(&RES1,Q0.nrows,C1->nrows,np);
	build_aux_mats2(&RES2,C1->nrows,C1->nrows,np);
	build_aux_mats2(&RES3,D12->nrows,C1->nrows,np);	
	build_aux_mats2(&RES4,C1->nrows,D12->nrows,np);	

	build_aux_mats2(&MF1,dF.nrows,C->ncols,np);
	build_aux_mats2(&MF2,dF.ncols,dF.ncols,np);
	build_aux_mats2(&MF3,dQ.ncols,dQ.ncols,np);
	build_aux_mats2(&MF4,dV.ncols,dV.ncols,np);
	/*CF=C1+D12*F*C*/
	double scale1=1.0;
	double scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&F0,C,&OUT,np);
	if(DEBUG_FNLSDP && np==0){
		printf("OUT:\n");
		print_genmatrix(&OUT);	
	}
	scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,D12,&OUT,&CF,np);
	if(DEBUG_FNLSDP && np==0){
		printf("CF:\n");
		print_genmatrix(&CF);	
	}
	
	/*CFT=(C1+D12*F*C)'*/
	copy_mat_gen_trans(&CF,&CFT,np);
	if(DEBUG_FNLSDP && np==0){
		printf("CFT:\n");
		print_genmatrix(&CFT);	
	}
	
	/*CF*dQ*CF'*/
	scale1=1.0;
	scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&dQ,&CFT,&RES1,np);
	mat_mult_raw_gen1(scale1,scale2,&CF,&RES1,&RES2,np);
	if(DEBUG_FNLSDP && np==0){
		printf("CF*dQ*CF':\n");
		print_genmatrix(&RES2);	
	}
	

	/*dCF=D12*dF*C*/
	scale1=1.0;
	scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&dF,C,&MF1,np);
	mat_mult_raw_gen1(scale1,scale2,D12,&MF1,&dCF,np);
	if(DEBUG_FNLSDP && np==0){
		printf("dCF:\n");
		print_genmatrix(&dCF);	
	}
	
	/*dCFT=(D12*dF*C)'*/
	copy_mat_gen_trans(&dCF,&dCFT,np);
	if(DEBUG_FNLSDP && np==0){
		printf("dCFT:\n");
		print_genmatrix(&dCFT);	
	}
	
	/*dCF*Q*CF'*/
	scale1=1.0;
	scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&Q0,&CFT,&OUT2,np);
	mat_mult_raw_gen1(scale1,scale2,&dCF,&OUT2,&RES3,np);
	if(DEBUG_FNLSDP && np==0){
		printf("dCF*Q*CF':\n");
		print_genmatrix(&RES3);	
	}

	/*CF*Q*dCF'*/
	scale1=1.0;
	scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&Q0,&dCFT,&OUT3,np);
	mat_mult_raw_gen1(scale1,scale2,&CF,&OUT3,&RES4,np);
	if(DEBUG_FNLSDP && np==0){
		printf("CF*Q*dCF':\n");
		print_genmatrix(&RES4);	
	}

	copy_mat_gen_trans(&dF,&dFT,np);
	if(DEBUG_FNLSDP && np==0){
		printf("dFT:\n");
		print_genmatrix(&dFT);	
	}
	copy_mat_gen_trans(&dQ,&dQT,np);
	if(DEBUG_FNLSDP && np==0){
		printf("dQT:\n");
		print_genmatrix(&dQT);	
	}
	copy_mat_gen_trans(&dV,&dVT,np);
	if(DEBUG_FNLSDP && np==0){
		printf("dVT:\n");
		print_genmatrix(&dVT);	
	}
	/*dF'*DF*/
	scale1=1.0;
	scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&dFT, &dF,&MF2,np);	
	if(DEBUG_FNLSDP && np==0){
		printf("dF'*dF:\n");
		print_genmatrix(&MF2);	
	}

	/*dQ'*DQ*/
	scale1=1.0;
	scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&dQT, &dQ, &MF3,np);	
	if(DEBUG_FNLSDP && np==0){
		printf("dQ'*dQ:\n");
		print_genmatrix(&MF3);	
	}
	
	/*dV'*DV*/
	scale1=1.0;
	scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&dVT, &dV,&MF4,np);	
	if(DEBUG_FNLSDP && np==0){
		printf("dV'*dV:\n");
		print_genmatrix(&MF4);	
	}
	/*
	printf("trace(dCF*Q*CF')=%f\n",trace(&RES3,np));	
	printf("trace(CF*Q*dCF')=%f\n",trace(&RES4,np));	
	printf("trace(CF*dQ*CF')=%f\n",trace(&RES2,np));	
	printf("trace(dF'*dF)=%f\n",trace(&MF2,np));	
	printf("trace(dQ'*dQ)=%f\n",trace(&MF3,np));	
	printf("trace(dV'*dV)=%f\n",trace(&MF4,np));	
	*/
	tt=trace(&RES3,np) + trace(&RES4,np) + trace(&RES2,np) + trace(&MF2,np) + trace(&MF3,np) + trace(&MF4,np);

	//printf("tt=%f\n",tt);

	free_mat_gen(&CF,np);
	free_mat_gen(&CFT,np);
	free_mat_gen(&dCF,np);
	free_mat_gen(&dCFT,np);

	free_mat_gen(&dF,np);
	free_mat_gen(&dFT,np);
	free_mat_gen(&dQ,np);
	free_mat_gen(&dQT,np);
	free_mat_gen(&dV,np);
	free_mat_gen(&dVT,np);

	free_mat_gen(&OUT,np);
	free_mat_gen(&OUT2,np);
	free_mat_gen(&OUT3,np);

	free_mat_gen(&RES1,np);
	free_mat_gen(&RES2,np);
	free_mat_gen(&MF1,np);
	free_mat_gen(&MF2,np);
	free_mat_gen(&MF3,np);
	free_mat_gen(&MF4,np);
	
	free_mat_gen(&F0,np);
	free_mat_gen(&V0,np);
	free_mat_gen(&Q0,np);

	return tt;
}

double trace(genmatrix *A,int np){
	double tt=0.0;
	if(A->nrows!=A->ncols){
		if(DEBUG_FNLSDP && np==0)printf("fnlsdp_objective_theta: trace: different rows and cols in trace calculation.\n");
		exit(1);	
	}
	int i=0;
	for(i=0;i<A->nrows;i++){
		tt=tt+A->mat[i+i*A->nrows];
	}
	return tt;
}

double eval_theta(genmatrix *F0,genmatrix *Q0, genmatrix *V0,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np){
	double norma,lambda1ret;
	genmatrix AF, V0invadd, OUT, B1T, B1B1T, AFQ, AFQT, AFV, AFVT, I, M1, M2;
	double scale1=1.0;
	double scale2=0.0;
	build_aux_mats1(A,&AF,-1,-1,np);
	build_aux_mats1(V0,&V0invadd,-1,-1,np);
	build_aux_mats2(&OUT,F0->nrows,C->ncols,np);
	build_aux_mats2(&B1B1T,B1->nrows,B1->nrows,np);
	build_aux_mats2(&AFQ,AF.nrows,Q0->ncols,np);
	build_aux_mats2(&AFQT,Q0->ncols,AF.nrows,np);
	build_aux_mats2(&AFV,AF.nrows,V0->ncols,np);
	build_aux_mats2(&AFVT,V0->ncols,AF.nrows,np);
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

	scale2=0.0;	
	copy_mat_gen_trans(B1,&B1T,np);
	copy_mat_gen_inverseadd(V0,&V0invadd,np);
	mat_mult_raw_gen1(scale1,scale2,B1,&B1T,&B1B1T,np);
	mat_mult_raw_gen1(scale1,scale2,&AF,Q0,&AFQ,np);
	mat_mult_raw_gen1(scale1,scale2,&AF,V0,&AFV,np);
	make_i_gen(&I,nx,nx,np);
	copy_mat_gen_trans(&AFQ,&AFQT,np);
	copy_mat_gen_trans(&AFV,&AFVT,np);
	
	build_aux_mats1(&B1B1T,&M1,-1,-1,np);	
	build_aux_mats1(&I,&M2,-1,-1,np);	

	scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,&AFQT,&I,&M1,np);
	mat_mult_raw_gen1(scale1,scale2,&AFQ,&I,&M1,np);
	
	mat_mult_raw_gen1(scale1,scale2,&AFVT,&I,&M2,np);
	mat_mult_raw_gen1(scale1,scale2,&AFV,&I,&M2,np);

	lambda1ret=lambda1(&V0invadd,nx,scapack,params,printlevel,np); 
	lambda1ret=max(0.0,lambda1ret);

	norma=norm2(nx,M1.mat)+norm2(nx,M2.mat)+lambda1ret;

	free_mat_gen(&AF,np);
	free_mat_gen(&V0invadd,np);
	free_mat_gen(&OUT,np);
	free_mat_gen(&B1T,np);
	free_mat_gen(&B1B1T,np);
	free_mat_gen(&AFQ,np);
	free_mat_gen(&AFQT,np);
	free_mat_gen(&AFV,np);
	free_mat_gen(&AFVT,np);
	free_mat_gen(&M1,np);
	free_mat_gen(&M2,np);
	free_mat_gen(&I,np);
	
	return norma;
}

double eval_theta_vec(double *X/*genmatrix *F0,genmatrix *Q0, genmatrix *V0*/,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np){
	double norma=0.0,lambda1ret;
	int i,j;
	genmatrix F0,Q0,V0;
	build_aux_mats2(&Q0,nx,nx,np);
	build_aux_mats2(&F0,nu,ny,np);
	build_aux_mats2(&V0,nx,nx,np);

	for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			F0.mat[i+j*nu]=X[i+j*nu];
	}

	int count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			Q0.mat[i+j*nx]=X[count+nu*ny];
			if(i!=j)
				Q0.mat[j+i*nx]=X[count+nu*ny];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			V0.mat[i+j*nx]=X[count+nu*ny+nx*(nx+1)/2];
			if(i!=j)
				V0.mat[j+i*nx]=X[count+nu*ny+nx*(nx+1)/2];
			count++;
		}
	}

	genmatrix AF, V0invadd, OUT, B1T, B1B1T, AFQ, AFQT, AFV, AFVT, I, M1, M2;
	double scale1=1.0;
	double scale2=0.0;
	build_aux_mats1(A,&AF,-1,-1,np);
	build_aux_mats1(&V0,&V0invadd,-1,-1,np);
	build_aux_mats2(&OUT,F0.nrows,C->ncols,np);
	build_aux_mats2(&B1B1T,B1->nrows,B1->nrows,np);
	build_aux_mats2(&AFQ,AF.nrows,Q0.ncols,np);
	build_aux_mats2(&AFQT,Q0.ncols,AF.nrows,np);
	build_aux_mats2(&AFV,AF.nrows,V0.ncols,np);
	build_aux_mats2(&AFVT,V0.ncols,AF.nrows,np);
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

	scale2=0.0;	
	copy_mat_gen_trans(B1,&B1T,np);
	copy_mat_gen_inverseadd(&V0,&V0invadd,np);
	mat_mult_raw_gen1(scale1,scale2,B1,&B1T,&B1B1T,np);
	mat_mult_raw_gen1(scale1,scale2,&AF,&Q0,&AFQ,np);
	mat_mult_raw_gen1(scale1,scale2,&AF,&V0,&AFV,np);
	make_i_gen(&I,nx,nx,np);
	copy_mat_gen_trans(&AFQ,&AFQT,np);
	copy_mat_gen_trans(&AFV,&AFVT,np);
	
	build_aux_mats1(&B1B1T,&M1,-1,-1,np);	
	build_aux_mats1(&I,&M2,-1,-1,np);	

	scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,&AFQT,&I,&M1,np);
	mat_mult_raw_gen1(scale1,scale2,&AFQ,&I,&M1,np);
	
	mat_mult_raw_gen1(scale1,scale2,&AFVT,&I,&M2,np);
	mat_mult_raw_gen1(scale1,scale2,&AFV,&I,&M2,np);

	lambda1ret=lambda1(&V0invadd,nx,scapack,params,printlevel,np); 
	lambda1ret=max(0.0,lambda1ret);

	double restr1=norm2(nx,M1.mat);
	double restr2=norm2(nx,M2.mat);  

	norma=restr1+restr2+lambda1ret;


	outt->filter_cand_t=norma;
	outt->filter_cand_constraints_restr1=restr1;
	outt->filter_cand_constraints_restr2=restr2;
	outt->filter_cand_constraints_restr3=lambda1ret;


	free_mat_gen(&AF,np);
	free_mat_gen(&V0invadd,np);
	free_mat_gen(&OUT,np);
	free_mat_gen(&B1T,np);
	free_mat_gen(&B1B1T,np);
	free_mat_gen(&AFQ,np);
	free_mat_gen(&AFQT,np);
	free_mat_gen(&AFV,np);
	free_mat_gen(&AFVT,np);
	free_mat_gen(&M1,np);
	free_mat_gen(&M2,np);
	free_mat_gen(&I,np);

	free_mat_gen(&F0,np);
	free_mat_gen(&V0,np);
	free_mat_gen(&Q0,np);
	
	return norma;
}


double lambda1(genmatrix *a, int n, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np)
{
	double lambda1ret;
	//MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

	int ictxt, prow, pcol, myrow, mycol;
	int brow,bcol;
	int info, lwork, liwork;
	int *desca,*descz;
	double *a0, *z0, *w0, *querywork, *work, *gap;
	int *queryiwork, *iwork, *ifail, *iclustr;	

	int izero=0, ione=1;
	double mone=-1.0;
	double dzero=0.0;
	char jobz, uplo, range;
	int i,j;
	int LDA, LDB;
	double temp1, temp2;
	
	desca=(int *)malloc(9*sizeof(int));
	descz=(int *)malloc(9*sizeof(int));

	brow=2;//32;
	bcol=2;//32;
	prow=scapack.nprow;
	pcol=scapack.npcol;

	jobz='N';
	uplo='U';
	range='I';

	
	//Cblacs_get(0,0,&ictxt);
	//sl_init_(&ictxt,&prow,&pcol);
	//Cblacs_gridinit(&ictxt,"Row",prow, pcol);
	Cblacs_gridinfo(scapack.ic, &prow,&pcol, &myrow, &mycol);


	//Compute the size of the local matrices
	LDA=numroc_(&n, &brow, &myrow, &izero, &prow);
	LDB=numroc_(&n, &bcol, &mycol, &izero, &pcol);


	//allocat space for A, Z and W
	a0=(double *) malloc(LDA*LDB*sizeof(double));
	w0=(double *) malloc(n*sizeof(double));
	z0=(double *) malloc(LDA*LDB*sizeof(double));
	
	//initialize the array descriptor
	descinit_(desca, &n, &n, &brow, &bcol, &izero, &izero, /*&ictxt*/&scapack.ic, &LDA, &info);
	descinit_(descz, &n, &n, &brow, &bcol, &izero, &izero, /*&ictxt*/&scapack.ic, &LDA, &info);

	//distribute matrix to grid
	for(j=1;j<n+1;j++){
		for(i=1;i<n+1;i++){
			pdelset_(a0,&i,&j,desca,&(a->mat[i-1 + (j-1)*(a->nrows)]));
		}
	}

	char cmach = 'U';
	double abstol=-1.0;
	
	querywork=(double *)malloc(1*sizeof(double));
	lwork=-1;
	queryiwork=(int *)malloc(1*sizeof(int));
	liwork=-1;

	ifail=(int *)malloc(n*sizeof(int));
	iclustr=(int *)malloc(2*prow*pcol*sizeof(int));

	gap=(double *)malloc(prow*pcol*sizeof(double));

	int found=0;
	int nz=0;


	//Cblacs_barrier(/*ictxt*/scapack.ic,"All");
	MPI_Barrier(MPI_COMM_WORLD);
	
#ifdef NOUNDERLAPACK
  #ifdef CAPSLAPACK
	PDSYEVX(&jobz, &range, &uplo, &n, a0, &ione, &ione, desca, &dzero, &dzero, &n, &n, &mone, &found, &nz, w0, &mone, z0, &ione, &ione, descz, querywork, &lwork, queryiwork, &liwork, ifail, iclustr, gap, &info);
  #else
	pdsyevx(&jobz, &range, &uplo, &n, a0, &ione, &ione, desca, &dzero, &dzero, &n, &n, &mone, &found, &nz, w0, &mone, z0, &ione, &ione, descz, querywork, &lwork, queryiwork, &liwork, ifail, iclustr, gap, &info);
  #endif
#else
  #ifdef CAPSLAPACK
	PDSYEVX_(&jobz, &range, &uplo, &n, a0, &ione, &ione, desca, &dzero, &dzero, &n, &n, &mone, &found, &nz, w0, &mone, z0, &ione, &ione, descz, querywork, &lwork, queryiwork, &liwork, ifail, iclustr, gap, &info);
  #else
	pdsyevx_(&jobz, &range, &uplo, &n, a0, &ione, &ione, desca, &dzero, &dzero, &n, &n, &mone, &found, &nz, w0, &mone, z0, &ione, &ione, descz, querywork, &lwork, queryiwork, &liwork, ifail, iclustr, gap, &info);
  #endif
#endif
	MPI_Barrier(MPI_COMM_WORLD);

	lwork=((int)querywork[0]);//*sizeof(double);
	liwork=((int)queryiwork[0]);//*sizeof(int);
	

	iwork=(int *)malloc( liwork * sizeof(int));
	work=(double *)malloc( lwork*sizeof(double));//lwork*sizeof(double));

	Cblacs_barrier(/*ictxt*/scapack.ic,"All");
	
#ifdef NOUNDERLAPACK
  #ifdef CAPSLAPACK
	PDSYEVX(&jobz, &range, &uplo, &n, a0, &ione, &ione, desca, &dzero, &dzero, &n, &n, &mone, &found, &nz, w0, &mone, z0, &ione, &ione, descz, work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);
  #else
	pdsyevx(&jobz, &range, &uplo, &n, a0, &ione, &ione, desca, &dzero, &dzero, &n, &n, &mone, &found, &nz, w0, &mone, z0, &ione, &ione, descz, work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);
  #endif
#else
  #ifdef CAPSLAPACK
	PDSYEVX_(&jobz, &range, &uplo, &n, a0, &ione, &ione, desca, &dzero, &dzero, &n, &n, &mone, &found, &nz, w0, &mone, z0, &ione, &ione, descz, work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);
  #else
	pdsyevx_(&jobz, &range, &uplo, &n, a0, &ione, &ione, desca, &dzero, &dzero, &n, &n, &mone, &found, &nz, w0, &mone, z0, &ione, &ione, descz, work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);
  #endif
#endif
	/*if(myrow==0 && mycol==0){
		printf("e1=%f\n", w0[0]);
	}*/

	MPI_Barrier(MPI_COMM_WORLD);
	
	lambda1ret=w0[0];

	free(querywork);	free(queryiwork);
	free(gap);	free(iclustr);	free(ifail);
	free(iwork);	free(work);
	free(w0);	free(a0);	free(z0);
	free(desca);		free(descz);
      	
   	//if (myrow==0 && mycol ==0)
	//	Cblacs_gridexit(ictxt);
	//MPI_Finalize();
	return lambda1ret;
}

/*void get_theta_vec(double *X, double theta,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np){
	double norma=0.0,lambda1ret;
	int i,j;
	genmatrix F0,Q0,V0;
	build_aux_mats2(&Q0,nx,nx,np);
	build_aux_mats2(&F0,nu,ny,np);
	build_aux_mats2(&V0,nx,nx,np);

	for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			F0.mat[i+j*nu]=X[i+j*nu];
	}

	int count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			Q0.mat[i+j*nx]=X[count+nu*ny];
			if(i!=j)
				Q0.mat[j+i*nx]=X[count+nu*ny];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			V0.mat[i+j*nx]=X[count+nu*ny+nx*(nx+1)/2];
			if(i!=j)
				V0.mat[j+i*nx]=X[count+nu*ny+nx*(nx+1)/2];
			count++;
		}
	}

	genmatrix AF, V0invadd, OUT, B1T, B1B1T, AFQ, AFQT, AFV, AFVT, I, M1, M2;
	double scale1=1.0;
	double scale2=0.0;
	build_aux_mats1(A,&AF,-1,-1,np);
	build_aux_mats1(&V0,&V0invadd,-1,-1,np);
	build_aux_mats2(&OUT,F0.nrows,C->ncols,np);
	build_aux_mats2(&B1B1T,B1->nrows,B1->nrows,np);
	build_aux_mats2(&AFQ,AF.nrows,Q0.ncols,np);
	build_aux_mats2(&AFQT,Q0.ncols,AF.nrows,np);
	build_aux_mats2(&AFV,AF.nrows,V0.ncols,np);
	build_aux_mats2(&AFVT,V0.ncols,AF.nrows,np);
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

	scale2=0.0;	
	copy_mat_gen_trans(B1,&B1T,np);
	copy_mat_gen_inverseadd(&V0,&V0invadd,np);
	mat_mult_raw_gen1(scale1,scale2,B1,&B1T,&B1B1T,np);
	mat_mult_raw_gen1(scale1,scale2,&AF,&Q0,&AFQ,np);
	mat_mult_raw_gen1(scale1,scale2,&AF,&V0,&AFV,np);
	make_i_gen(&I,nx,nx,np);
	copy_mat_gen_trans(&AFQ,&AFQT,np);
	copy_mat_gen_trans(&AFV,&AFVT,np);
	
	build_aux_mats1(&B1B1T,&M1,-1,-1,np);	
	build_aux_mats1(&I,&M2,-1,-1,np);	

	scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,&AFQT,&I,&M1,np);
	mat_mult_raw_gen1(scale1,scale2,&AFQ,&I,&M1,np);
	
	mat_mult_raw_gen1(scale1,scale2,&AFVT,&I,&M2,np);
	mat_mult_raw_gen1(scale1,scale2,&AFV,&I,&M2,np);

	lambda1ret=lambda1(&V0invadd,nx,scapack,params,printlevel,np); 
	lambda1ret=max(0.0,lambda1ret);

	norma=norm2(nx,M1.mat)+norm2(nx,M2.mat)+lambda1ret;

	free_mat_gen(&AF,np);
	free_mat_gen(&V0invadd,np);
	free_mat_gen(&OUT,np);
	free_mat_gen(&B1T,np);
	free_mat_gen(&B1B1T,np);
	free_mat_gen(&AFQ,np);
	free_mat_gen(&AFQT,np);
	free_mat_gen(&AFV,np);
	free_mat_gen(&AFVT,np);
	free_mat_gen(&M1,np);
	free_mat_gen(&M2,np);
	free_mat_gen(&I,np);

	free_mat_gen(&F0,np);
	free_mat_gen(&V0,np);
	free_mat_gen(&Q0,np);
	
	return norma;
} */
