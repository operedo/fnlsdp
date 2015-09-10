/**
 * \file fnlsdp_build_mats.c 
 * \author Oscar Peredo
 * \date 09 Feb 2010
 *
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "declarations.h"

/** 
 * \brief Copy general matrix A to general matrix B (B memory must be alloc)
 * \param A
 * \param B 
 * \param np
 */
void copy_mat_gen(genmatrix *A,genmatrix *B,int np)
{
  	int i,j;
	//if(A.nrows!=B.nrows || A.ncols!=B.ncols){
	//	if(DEBUG_FNLSDP)printf("copy_mat_get: matrix dimensions must agree.");
	//	exit(-1);
	//}
  	B->mat=(double *)calloc((A->nrows)*(A->ncols),sizeof(double));
	for (i=0; i< (A->nrows)*(A->ncols); i++)
	{
		B->mat[i]=A->mat[i];
	}
	B->nrows=A->nrows;
	B->ncols=A->ncols;
}

/** 
 * \brief Copy general matrix A (transposed) to general matrix B (B memory must be alloc)
 * \param A
 * \param B
 * \param np
 *
 */
void copy_mat_gen_trans(genmatrix *A,genmatrix *B,int np)
{
  	int i,j;
	//if(A.nrows!=B.nrows || A.ncols!=B.ncols){
	//	if(DEBUG_FNLSDP)printf("copy_mat_get: matrix dimensions must agree.");
	//	exit(-1);
	//}
  	B->mat=(double *)calloc((A->nrows)*(A->ncols),sizeof(double));
	for (i=0; i< (A->nrows); i++)
		for (j=0; j< (A->ncols); j++)
			B->mat[j+i*A->ncols]=A->mat[i+j*A->nrows];
	B->nrows=A->ncols;
	B->ncols=A->nrows;
}

/** 
 * \brief Copy general matrix A (additive inverse) to general matrix B (B memory must be alloc)
 * \param A
 * \param B
 * \param np
 */
void copy_mat_gen_inverseadd(genmatrix *A,genmatrix *B,int np)
{
  	int i,j;
	//if(A.nrows!=B.nrows || A.ncols!=B.ncols){
	//	if(DEBUG_FNLSDP)printf("copy_mat_get: matrix dimensions must agree.");
	//	exit(-1);
	//}
  	B->mat=(double *)calloc((A->nrows)*(A->ncols),sizeof(double));
	for (i=0; i< (A->nrows); i++)
		for (j=0; j< (A->ncols); j++)
			B->mat[i+j*A->ncols]=-A->mat[i+j*A->nrows];
	B->nrows=A->nrows;
	B->ncols=A->ncols;
}


/** 
 * \brief Calculate alpha*A*B+beta*C with DGEMM
 * \param scale1 
 * \param scale2
 * \param ap
 * \param bp
 * \param cp
 * \param np
 */
void mat_mult_raw_gen1(scale1,scale2,ap,bp,cp,np)
     	//int m;
	//int n;
	//int k;
     	double scale1;
     	double scale2;
     	genmatrix *ap;
     	genmatrix *bp;
     	genmatrix *cp;
	int np;
{
	int n,m,k;
	m=ap->nrows;
	n=bp->ncols;
	k=ap->ncols;

	if(ap->nrows!=cp->nrows || bp->ncols!=cp->ncols || ap->ncols!=bp->nrows){
		if(DEBUG_FNLSDP && np==0)printf("mat_mult_raw_gen1: matrix dimensions must agree.");
		exit(-1);
	}
#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
	  DGEMM("N","N",&m,&n,&k,&scale1,ap->mat,&m,bp->mat,&k,&scale2,cp->mat,&m);
#else
	  dgemm("N","N",&m,&n,&k,&scale1,ap->mat,&m,bp->mat,&k,&scale2,cp->mat,&m);
#endif
#else
#ifdef CAPSBLAS
	  DGEMM_("N","N",&m,&n,&k,&scale1,ap->mat,&m,bp->mat,&k,&scale2,cp->mat,&m);
#else
	  dgemm_("N","N",&m,&n,&k,&scale1,ap->mat,&m,bp->mat,&k,&scale2,cp->mat,&m);
#endif
#endif

}

/** 
 * \brief Calculate alpha*A*B'+beta*C with DGEMM
 * \param scale1 
 * \param scale2
 * \param ap
 * \param bp
 * \param cp
 * \param np
 */
void mat_mult_raw_gen2(scale1,scale2,ap,bp,cp,np)
     	//int m;
	//int n;
	//int k;
     	double scale1;
     	double scale2;
     	genmatrix *ap;
     	genmatrix *bp;
     	genmatrix *cp;
	int np;
{
	int n,m,k;
	m=ap->nrows;
	n=bp->ncols;
	k=ap->ncols;

	//if(ap->nrows!=cp->nrows || bp->ncols!=cp->ncols || ap->ncols!=bp->nrows){
	if(ap->nrows!=cp->nrows || bp->nrows!=cp->ncols || ap->ncols!=bp->ncols){
		if(DEBUG_FNLSDP && np==0)printf("mat_mult_raw_gen2: matrix dimensions must agree.");
		exit(-1);
	}
#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
	  DGEMM("N","T",&m,&n,&k,&scale1,ap->mat,&m,bp->mat,&n,&scale2,cp->mat,&m);
#else
	  dgemm("N","T",&m,&n,&k,&scale1,ap->mat,&m,bp->mat,&n,&scale2,cp->mat,&m);
#endif
#else
#ifdef CAPSBLAS
	  DGEMM_("N","T",&m,&n,&k,&scale1,ap->mat,&m,bp->mat,&n,&scale2,cp->mat,&m);
#else
	  dgemm_("N","T",&m,&n,&k,&scale1,ap->mat,&m,bp->mat,&n,&scale2,cp->mat,&m);
#endif
#endif

}

/** 
 * \brief Calculate alpha*A'*B+beta*C with DGEMM
 * \param scale1 
 * \param scale2
 * \param ap
 * \param bp
 * \param cp
 * \param np
 */
void mat_mult_raw_gen3(scale1,scale2,ap,bp,cp,np)
     	//int m;
	//int n;
	//int k;
     	double scale1;
     	double scale2;
     	genmatrix *ap;
     	genmatrix *bp;
     	genmatrix *cp;
	int np;
{
	int n,m,k;
	m=ap->nrows;
	n=bp->ncols;
	k=ap->ncols;

	//if(ap->nrows!=cp->nrows || bp->ncols!=cp->ncols || ap->ncols!=bp->nrows){
	if(ap->ncols!=cp->nrows || bp->ncols!=cp->ncols || ap->nrows!=bp->nrows){
		if(DEBUG_FNLSDP && np==0)printf("mat_mult_raw_gen3: matrix dimensions must agree.");
		exit(-1);
	}
#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
	  DGEMM("T","N",&m,&n,&k,&scale1,ap->mat,&k,bp->mat,&k,&scale2,cp->mat,&m);
#else
	  dgemm("T","N",&m,&n,&k,&scale1,ap->mat,&k,bp->mat,&k,&scale2,cp->mat,&m);
#endif
#else
#ifdef CAPSBLAS
	  DGEMM_("T","N",&m,&n,&k,&scale1,ap->mat,&k,bp->mat,&k,&scale2,cp->mat,&m);
#else
	  dgemm_("T","N",&m,&n,&k,&scale1,ap->mat,&k,bp->mat,&k,&scale2,cp->mat,&m);
#endif
#endif

}

/**
 * \brief Allocate out matrix with dimensions data->nrows x data->ncols, or copy from data matrix with the same dimensions
 * \param data
 * \param out
 * \param nrows
 * \param ncols
 * \param np
 */
void build_aux_mats1(genmatrix *data, genmatrix *out, int nrows, int ncols, int np){
	if(nrows==-1 && ncols==-1)
		copy_mat_gen(data,out,np);
	else{ 
		if(nrows==0 && ncols==0){
			out->mat=(double *)calloc((data->nrows)*(data->ncols),sizeof(double));
		}
		else{
			if(DEBUG_FNLSDP && np==0)printf("build_aux_mats1: bad arguments.\n");
		}
	}	
}

/**
 * \brief Allocate out matrix with dimensions nrows x ncols
 * \param data
 * \param out
 * \param nrows
 * \param ncols
 * \param np
 */
void build_aux_mats2(genmatrix *out, int nrows, int ncols, int np){
	if(nrows>0 && ncols>0){
		out->mat=(double *)calloc(nrows*ncols,sizeof(double));
		out->nrows=nrows;
		out->ncols=ncols;
	}
	else{
		if(DEBUG_FNLSDP && np==0)printf("build_aux_mats2: bad arguments.\n");
	} 
}

/*void build_aux_mats(genmatrix *A,genmatrix *AF, genmatrix *C1, genmatrix *CF, genmatrix *F0, genmatrix *OUT){
	copy_mat_gen(A,AF);
	copy_mat_gen(C1,CF);
	OUT->mat=(double *)calloc((F0->nrows)*(C1->ncols),sizeof(double));
	OUT->nrows=F0->nrows;
	OUT->ncols=C1->ncols;
}*/

/**
 * \brief Allocate identity matrix
 * \param I
 * \param nrows
 * \param ncols
 * \param np
 */
void make_i_gen(genmatrix *I,int nrows, int ncols,int np)
{
  	int k,i,j;
	I->mat=(double *)calloc(nrows*ncols,sizeof(double));
	I->nrows=nrows;
	I->ncols=ncols;
	for (k=1; k<=nrows*ncols; k++)
    	{
		i=((k-1)%nrows)+1;
		j=floor((k-1)/nrows)+1;
		if(i==j)
			I->mat[k-1]=1.0;
    	}
}

void make_diag(genmatrix *M, struct blockmatrix *A,int np)
{
  	int blk,i;
	double *p;
  	A->nblocks=1;
  	A->blocks=(struct blockrec *)malloc(sizeof(struct blockrec)*(2));
  	if (A->blocks == NULL)
    	{
      		printf("make_diag: Storage allocation failed!\n");
      		exit(10);
    	};
	A->blocks[1].blockcategory=DIAG;
	A->blocks[1].blocksize=(M->nrows)*(M->nrows);
	A->blocks[1].data.vec=(double *)malloc(sizeof(double)*(M->nrows*M->nrows+1));
	if (A->blocks[1].data.vec == NULL)
	{
		printf("make_diag: Storage allocation failed!\n");
		exit(10);
	};
	p=A->blocks[1].data.vec;
	for (i=1; i<=A->blocks[1].blocksize; i++)
		p[i]=M->mat[i-1];
}


void build_C(struct blockmatrix *pC, int nx, int ny, int nu, genmatrix* Q0, genmatrix *V0, genmatrix *B1, genmatrix *AF,double rho,int np){
	genmatrix MF1, MF2;
	genmatrix B1T, AFT;
	struct blockmatrix diagMF1, diagMF2;
	double scale1, scale2;
	int i,j,blk,blksz,nblocks;

	nblocks=4;
	pC->nblocks=nblocks;
	pC->blocks=(struct blockrec *)malloc((nblocks+1)*sizeof(struct blockrec));
  	if (pC->blocks == NULL)
    	{
      		printf("build_C: Storage allocation failed!\n");
      		exit(10);
    	}

	build_aux_mats2(&MF1,nx,nx,np);
	scale1=1.0; scale2=0.0;
	copy_mat_gen_trans(B1,&B1T,np);
	mat_mult_raw_gen1(scale1,scale2,B1,&B1T,&MF1,np);
	scale1=1.0; scale2=1.0;
	copy_mat_gen_trans(AF,&AFT,np);
	mat_mult_raw_gen1(scale1,scale2,Q0,&AFT,&MF1,np);
	scale1=1.0; scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,AF,Q0,&MF1,np);
	if(DEBUG_FNLSDP && np==0){printf("MF1:\n");
	print_genmatrix(&MF1);}	

	make_i_gen(&MF2,nx,nx,np);
	scale1=1.0; scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,V0,&AFT,&MF2,np);
	scale1=1.0; scale2=1.0;
	mat_mult_raw_gen1(scale1,scale2,AF,V0,&MF2,np);
	if(DEBUG_FNLSDP && np==0){printf("MF2:\n");
	print_genmatrix(&MF2);}
	
	/*
	* First Block: 
	* blkdiag(
	* -->	[I,zero;zero',0],
	*	diag([MF1(:);-MF1(:);MF2(:);-MF2(:)]),
	*	V0,
	*	[I,zero;zero',rho]
	*	)
	*/


	blk=1;
	blksz=nu*ny+nx*(nx+1)+1;
	pC->blocks[blk].blocksize=blksz;
	pC->blocks[blk].blockcategory=/*DIAG;*/MATRIX;
	pC->blocks[blk].data.mat=(double *)calloc((/*1+abs(blksz)*/blksz*blksz),sizeof(double));
	if (pC->blocks[blk].data.mat == NULL)
	{
		printf("build_C:Storage allocation failed!\n");
		exit(10);
	};
	for (i=1; i<=blksz; i++){
		for(j=1;j<=blksz;j++){
			if(i<blksz && i==j)
			//pC->blocks[blk].data.vec[i]=-1.0;
				pC->blocks[blk].data.mat[ijtok(i,i,blksz)]=-1.0;
			else
				pC->blocks[blk].data.mat[ijtok(i,j,blksz)]=0.0;
		}
	}

	/*
	* Mid Block: 
	* blkdiag(
	*	[I,zero;zero',0],
	* -->	diag([MF1(:);-MF1(:);MF2(:);-MF2(:)]),
	*	      ------
	*	V0,
	* -->	[I,zero;zero',rho]
	*	)
	*/	
	blk=2;
	blksz=4*(MF1.nrows)*(MF1.nrows); 
	pC->blocks[blk].blocksize=blksz;
	pC->blocks[blk].blockcategory=DIAG;
	pC->blocks[blk].data.vec=(double *)calloc((1+abs(blksz)),sizeof(double));
	if (pC->blocks[blk].data.vec == NULL)
	{
		printf("build_C: Storage allocation failed!\n");
		exit(10);
	};
	for (i=1; i<=blksz; i++){
		if(i>=1 && i<=MF1.nrows*MF1.nrows){
			if(MF1.mat[i-1]!=0.0)
				pC->blocks[blk].data.vec[i]=-MF1.mat[i-1];
		}
		if(i>MF1.nrows*MF1.nrows && i<=2*MF1.nrows*MF1.nrows){
			if(MF1.mat[i-(1+MF1.nrows*MF1.nrows)]!=0.0)
				pC->blocks[blk].data.vec[i]=-(-MF1.mat[i-(1+MF1.nrows*MF1.nrows)]);
		}
		if(i>2*MF1.nrows*MF1.nrows && i<=3*MF1.nrows*MF1.nrows){
			if(MF2.mat[i-(1+2*MF1.nrows*MF1.nrows)]!=0.0)
				pC->blocks[blk].data.vec[i]=-MF2.mat[i-(1+2*MF1.nrows*MF1.nrows)];
		}
		if(i>3*MF1.nrows*MF1.nrows && i<=4*MF1.nrows*MF1.nrows){
			if(MF2.mat[i-(1+3*MF1.nrows*MF1.nrows)]!=0.0)
				pC->blocks[blk].data.vec[i]=-(-MF2.mat[i-(1+3*MF1.nrows*MF1.nrows)]);
		}
	}

	/*
	* Mid Block: 
	* blkdiag(
	*	[I,zero;zero',0],
	* 	diag([MF1(:);-MF1(:);MF2(:);-MF2(:)]),
	* -->	V0,
	* 	[I,zero;zero',rho]
	*	)
	*/	
	blk=3;
	blksz=V0->nrows;
	pC->blocks[blk].blocksize=blksz;
	pC->blocks[blk].blockcategory=MATRIX;
	pC->blocks[blk].data.mat=(double *)calloc((blksz*blksz),sizeof(double));
	if (pC->blocks[blk].data.mat == NULL)
	{
		printf("build_C:Storage allocation failed!\n");
		exit(10);
	};
	for (j=1; j<=blksz; j++)
		for (i=j; i<=blksz; i++)
			//if(V0->mat[ijtok(i,j,blksz)]!=0.0)
				pC->blocks[blk].data.mat[ijtok(i,j,blksz)]=-V0->mat[ijtok(i,j,blksz)];
	/*
	* Last Block: 
	* blkdiag(
	*	[I,zero;zero',0],
	*	diag([MF1(:);-MF1(:);MF2(:);-MF2(:)]),
	*	V0,
	* -->	[I,zero;zero',rho]
	*	)
	*/

	blk=4;
	blksz=nu*ny+nx*(nx+1)+1;
	pC->blocks[blk].blocksize=blksz;
	pC->blocks[blk].blockcategory=MATRIX;
	pC->blocks[blk].data.mat=(double *)calloc((blksz*blksz),sizeof(double));
	if (pC->blocks[blk].data.mat == NULL)
	{
		printf("build_C:Storage allocation failed!\n");
		exit(10);
	};
	for (i=1; i<=blksz; i++){
		if(i<blksz)
			pC->blocks[blk].data.mat[ijtok(i,i,blksz)]=-1.0;
		else
			pC->blocks[blk].data.mat[ijtok(i,i,blksz)]=-rho*rho*(nu*ny+nx*(nx+1));
	}	
	
	free_mat_gen(&B1T,np);
	free_mat_gen(&AFT,np);
	
	free_mat_gen(&MF1,np);
	free_mat_gen(&MF2,np);		
}

void build_a(double **a, int nx, int ny, int nu,int np){
	int i;
	*a=(double *)calloc((nu*ny+nx*(nx+1)+1+1),sizeof(double));
	for(i=1;i<=nu*ny+nx*(nx+1)+1;i++){
		if(i<nu*ny+nx*(nx+1)+1)
			(*a)[i]=0.0;
		if(i==nu*ny+nx*(nx+1)+1)
			(*a)[i]=1.0;
	}
}


void build_constraints(struct constraintmatrix **pconstraints, int nx, int ny, int nz,int nu, genmatrix* Q0, genmatrix *V0, genmatrix *B, genmatrix *C,genmatrix *D12, genmatrix *AF,genmatrix * CF,int np){
  	struct sparseblock *p;
  	struct sparseblock *oldp;

	genmatrix CT,D12T,AUX1, AUX2, AUX3, AUX4, AA;	
	genmatrix Ctilde,Ctildep,CtildeT,CtildepT,BT, CFT;
	double *gamma;
	double scale1, scale2;
	int i,ib,j,ii,iii,jj,k,kaux,kk,kkk,pk,blk,blksz,nblocks;

	pk=nu*ny+nx*(nx+1)+1;


	/*
	* Build Gamma=[2*vec(D12'*CF*Q0*C');vec(CF'*CF);0_nx_x_nx]
	*/
	copy_mat_gen_trans(C,&CT,np);
	copy_mat_gen_trans(D12,&D12T,np);
	build_aux_mats2(&AUX1,nx,ny,np);
	build_aux_mats2(&AUX2,nz,ny,np);
	build_aux_mats2(&AUX3,nu,ny,np);
	build_aux_mats2(&AUX4,nx,nx,np);
	
	/*printf("Q0:\n");
	print_genmatrix(Q0);
	printf("C:\n");
	print_genmatrix(C);*/
	scale1=1.0; scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,Q0,&CT,&AUX1,np);
	if(DEBUG_FNLSDP && np==0){printf("AUX1:\n");
	print_genmatrix(&AUX1);}
	scale1=1.0; scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,CF,&AUX1,&AUX2,np);
	if(DEBUG_FNLSDP && np==0){printf("AUX2:\n");
	print_genmatrix(&AUX2);}
	scale1=1.0; scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,&D12T,&AUX2,&AUX3,np);
	if(DEBUG_FNLSDP && np==0){printf("AUX3:\n");
	print_genmatrix(&AUX3);}

	scale1=1.0; scale2=0.0;
	copy_mat_gen_trans(CF,&CFT,np);
	mat_mult_raw_gen1(scale1,scale2,&CFT,CF,&AUX4,np);

	if(DEBUG_FNLSDP && np==0){printf("CF'*CF:\n");
	print_genmatrix(&AUX4);}

	double *AUX4_sym=(double *)calloc((int)(nx*(nx+1)/2),sizeof(double));
	int counter=0;
	for(j=1;j<=nx;j++){
		for(i=j;i<=nx;i++){
			if(i==j){
				if(DEBUG_FNLSDP && np==0)printf("ijtok(%d,%d,%d)=%d\n",i,j,nx,(int)ijtok(i,j,nx));
				AUX4_sym[counter]=AUX4.mat[(int)ijtok(i,j,nx)];
				counter++;	
			}
			if(i>j){
				if(DEBUG_FNLSDP && np==0)printf("ijtok(%d,%d,%d)=%d\n",i,j,nx,(int)ijtok(i,j,nx));
				AUX4_sym[counter]=2*AUX4.mat[(int)ijtok(i,j,nx)];
				counter++;	
			}
		}
	}

	gamma=(double *)calloc((nu*ny+nx*(nx+1)),sizeof(double));
	if(DEBUG_FNLSDP && np==0)printf("gamma:\n");
	for(i=0;i<nu*ny+nx*(nx+1);i++){
		if(i<nu*ny)
			gamma[i]=2*AUX3.mat[i];
		if(i<nu*ny+(int)(nx*(nx+1)/2) && i>=nu*ny)
			gamma[i]=AUX4_sym[i-nu*ny];
		if(i>=nu*ny+(int)(nx*(nx+1)/2) )
			gamma[i]=0.0;
		if(DEBUG_FNLSDP && np==0)printf("gamma[%d]=%f\n",i,gamma[i]);
	}

	free(AUX4_sym);

	/*
	* Build AA=[]
	*/

	build_aux_mats2(&Ctilde,ny,nx,np);
	scale1=1.0; scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,C,Q0,&Ctilde,np);
	if(DEBUG_FNLSDP && np==0){printf("Ctilde:\n");
	print_genmatrix(&Ctilde);}
	copy_mat_gen_trans(&Ctilde,&CtildeT,np);
	if(DEBUG_FNLSDP && np==0){printf("CtildeT:\n");
	print_genmatrix(&CtildeT);}
	build_aux_mats2(&Ctildep,ny,nx,np);
	scale1=1.0; scale2=0.0;
	mat_mult_raw_gen1(scale1,scale2,C,V0,&Ctildep,np);
	if(DEBUG_FNLSDP && np==0){printf("Ctildep:\n");
	print_genmatrix(&Ctildep);}
	copy_mat_gen_trans(B,&BT,np);
	if(DEBUG_FNLSDP && np==0){printf("BT:\n");
	print_genmatrix(&BT);}
	copy_mat_gen_trans(&Ctildep,&CtildepT,np);
	if(DEBUG_FNLSDP && np==0){printf("CtildepT:\n");
	print_genmatrix(&CtildepT);}


	

	/*
	* Build constraints...
	*/
	
  	*pconstraints=(struct constraintmatrix *)malloc((pk+1)*sizeof(struct constraintmatrix));

  	if (*pconstraints == NULL)
    	{
      		printf("Storage allocation failed!\n");
      		exit(10);
    	}
  
  	/*
   	* Null out all pointers in constraints.
   	*/
  	for (i=1; i<=pk; i++)
    	{
      		(*pconstraints)[i].blocks=NULL;
    	}

 	/*
   	* Now, go through each of the blks in each of the constraint matrices,
   	* and allocate space for the entries and indices.
   	*/
	int numentries;//=3+4*nx*nx*(nu*ny+2*nx*nx)+4+2;
	int totalsize=2*(nu*ny+nx*(nx+1)+1)+4*nx*nx+nx;
	int numblocks=1;

	int block;
	int endblocks;


	genmatrix index, index2;
	build_aux_mats2(&index,nx,nx,np);
	build_aux_mats2(&index2,(int)(nx*(nx+1)/2),2,np);
	
	counter=1;
	for(j=1;j<=nx;j++){
		for(i=j;i<=nx;i++){
			index.mat[ijtok(i,j,nx)]=counter;
			counter++;
		}
	}
	for(i=1;i<=nx;i++){
		for(j=i+1;j<=nx;j++){
			index.mat[ijtok(i,j,nx)]=index.mat[ijtok(j,i,nx)];
			counter++;
		}
	}
	

	counter=1;
	for(j=1;j<=nx;j++){
		for(i=j;i<=nx;i++){
			index2.mat[ijtok(counter,1,(int)(nx*(nx+1)/2))]=i;	
			index2.mat[ijtok(counter,2,(int)(nx*(nx+1)/2))]=j;	
			counter++;
		}
	}


	if(DEBUG_FNLSDP && np==0){
		printf("index:\n");
		print_genmatrix(&index);
		printf("index2:\n");
		print_genmatrix(&index2);
	}

	genmatrix TT3,TTT1,TTT2;
	build_aux_mats2(&TT3,nx*nx,(int)(nx*(nx+1)/2),np);
	build_aux_mats2(&TTT1,nx,(int)(nx*(nx+1)/2),np);
	build_aux_mats2(&TTT2,nx,(int)(nx*(nx+1)/2),np);

	int col=0;
	for(col=1;col<=nx;col++){
		for(i=1;i<=nx;i++){
			for(j=1;j<=(int)(nx*(nx+1)/2);j++){
				TTT1.mat[ijtok(i,j,nx)]=0.0;
				TTT2.mat[ijtok(i,j,nx)]=0.0;
			}
		} 
		for(j=1;j<=nx;j++){
			for(i=1;i<=nx;i++){
				TTT1.mat[ijtok(j,(int)index.mat[ijtok(i,j,nx)],nx)]=AF->mat[ijtok(col,i,nx)];
				if(DEBUG_FNLSDP && np==0)printf("TTT1.mat[%d,%d]=%f\n",j,(int)index.mat[ijtok(i,j,nx)],TTT1.mat[ijtok(j,(int)index.mat[ijtok(i,j,nx)],nx)]);
				TTT2.mat[ijtok(j,(int)index.mat[ijtok(col,i,nx)],nx)]=AF->mat[ijtok(j,i,nx)];
				if(DEBUG_FNLSDP && np==0)printf("TTT2.mat[%d,%d]=%f\n",j,(int)index.mat[ijtok(col,i,nx)],TTT2.mat[ijtok(j,(int)index.mat[ijtok(col,i,nx)],nx)]);
			}
		}
		for(i=1;i<=nx;i++){
			for(j=1;j<=(int)(nx*(nx+1)/2);j++){
				//printf("TT3.mat[ijtok(%d,%d,%d)]\n",i+nx*(col-1),j,nx*nx);
				TT3.mat[ijtok(i+nx*(col-1),j,nx*nx)]=TTT1.mat[ijtok(i,j,nx)] + TTT2.mat[ijtok(i,j,nx)];
			}
		}
	}
	

/*	for(col=1;col<=nx;col++){
		for(i=1;i<=nx;i++){
			for(j=1;j<=nx;j++){
				TT3.mat[ijtok(i+,j,nx*nx)]=TTT1.mat[ijtok(i,j)] + TTT2.mat[ijtok(i,j)];
				//TTT1.mat[ijtok(j,(int)index.mat[ijtok(i,j,nx)],nx)]=AF->mat[ijtok(col,i,nx)];
				//TTT2.mat[ijtok(j,(int)index.mat[ijtok(col,i,nx)],nx)]=AF->mat[ijtok(j,i,nx)];
			}
		}
 	}
*/	
	if(DEBUG_FNLSDP && np==0){
		printf("TT3:\n");
		print_genmatrix(&TT3);
	}
	



//////////////////////////

  	for (i=1; i<=pk; i++)
    	{
		if(i<pk)
			numblocks=4;
		else
			numblocks=1;
		
		block=1;
		if(DEBUG_FNLSDP && np==0)printf("numblocks=%d\n",numblocks);
	
		while (block<=numblocks)
		{
			if(DEBUG_FNLSDP && np==0)printf("allocating block %d from constraint %d\n",block,i);
			p=(struct sparseblock *)malloc(sizeof(struct sparseblock));
      			if (p==NULL)
			{
	  			if(DEBUG_FNLSDP && np==0)printf("Storage allocation failed!\n");
	  			exit(10);
			};
      			p->constraintnum=i;
      			p->blocknum=block;
      			p->numentries=1;
      			p->next=NULL;
      			p->entries=NULL;
      			p->iindices=NULL;
      			p->jindices=NULL;
      			p->blocksize=1;
      			if(block==1)
				(*pconstraints)[i].blocks=p;
			
			if(block>1)
				oldp->next=p;
			
	  		if(block==1){
				if(numblocks==4){
					p->blocksize=nu*ny+nx*(nx+1)+1;
					p->numentries=2;//3;
				}
				if(numblocks==1){
					p->blocksize=nu*ny+nx*(nx+1)+1;
					p->numentries=1;
				}
			}
			if(block==4){
				p->blocksize=nu*ny+nx*(nx+1)+1;
				p->numentries=1;//2;
			}
			if(block==3){
				p->blocksize=nx;
				p->numentries=1;//3;
			}
			if(block==2){
				p->blocksize=4*nx*nx;
				p->numentries=4*nx*nx;
			}
			p->entries=(double *)calloc((p->numentries+1),sizeof(double));
          		if (p->entries == NULL)
	    		{
	      			if(DEBUG_FNLSDP && np==0)printf("Storage allocation failed!\n");
	      			exit(10);
	    		}
			
#ifdef NOSHORTS
	  		p->iindices=(int *)calloc((p->numentries+1),sizeof(int));
#else
	 		p->iindices=(unsigned short *)calloc((p->numentries+1),sizeof(unsigned short));
#endif
          		if (p->iindices == NULL)
	    		{
	      			if(DEBUG_FNLSDP && np==0)printf("Storage allocation failed!\n");
	      			exit(10);
	    		}

#ifdef NOSHORTS
	  		p->jindices=(int *)calloc((p->numentries+1),sizeof(int));
#else
	  		p->jindices=(unsigned short *)calloc((p->numentries+1),sizeof(unsigned short));
#endif
          		if (p->jindices == NULL)
	    		{
	      			if(DEBUG_FNLSDP && np==0)printf("Storage allocation failed!\n");
	      			exit(10);
	    		}

			
			// Block 1 to nu*ny+nx*(nx+1)
			//
    			//for j=1:(nu*ny+nx*(nx+1))
        		//	II=eye(nu*ny+nx*(nx+1),nu*ny+nx*(nx+1));
        		//	ej=II(:,j);
        		//	Cero=zeros(nu*ny+nx*(nx+1),nu*ny+nx*(nx+1));
        		//	Aj=A(:,j);
        		//	Ej=E(:,:,j);
        		//	F(:,:,j)=blkdiag([Cero,ej;ej',gamma(j)],diag(Aj),Ej,[Cero,ej;ej',0]);
    			//end

				
			numentries=1;
			if(block==1){
				if(numblocks==4){
					p->entries[numentries]=gamma[i-1];
	      				p->iindices[numentries]=nu*ny+nx*(nx+1)+1;
	      				p->jindices[numentries]=nu*ny+nx*(nx+1)+1;	
					if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,p->iindices[numentries],p->jindices[numentries],p->entries[numentries]);
					numentries++;
					p->entries[numentries]=1.0;
	      				p->iindices[numentries]=i;
	      				p->jindices[numentries]=nu*ny+nx*(nx+1)+1;
					if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,p->iindices[numentries],p->jindices[numentries],p->entries[numentries]);
					numentries++;
					//p->entries[numentries]=1.0;
	      				//p->iindices[numentries]=nu*ny+2*nx*nx+1;
	      				//p->jindices[numentries]=i;
					//numentries++;
				}
				if(numblocks==1){
	      				p->entries[numentries]=1.0;
	      				p->iindices[numentries]=nu*ny+nx*(nx+1)+1;
	      				p->jindices[numentries]=nu*ny+nx*(nx+1)+1;
					if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,p->iindices[numentries],p->jindices[numentries],p->entries[numentries]);
					numentries++;
				}	
			}
			
			if(block==3){
				if(i>=nu*ny+(int)(nx*(nx+1)/2)+1){
					kk=i % (nu*ny+(int)(nx*(nx+1)/2));
					ii=(int)index2.mat[ijtok(kk,1,(int)(nx*(nx+1)/2))];//((kk-1) % nx)+1;
					jj=(int)index2.mat[ijtok(kk,2,(int)(nx*(nx+1)/2))];//(floor((kk-1)/(nx))+1);
					//printf("block 3: ii=%d, jj=%d, p->blocksize=%d\n",ii,jj,p->blocksize);
					if(ii==jj){
						p->entries[numentries]=1.0;
	      					p->iindices[numentries]=ii;
	      					p->jindices[numentries]=ii;
						if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,p->iindices[numentries],p->jindices[numentries],p->entries[numentries]);
						numentries++;
						/*p->entries[numentries]=0.0;
	      					p->iindices[numentries]=1;
	      					p->jindices[numentries]=1;
						if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,p->iindices[numentries],p->jindices[numentries],p->entries[numentries]);
						numentries++;*/
						//p->entries[numentries]=0.0;
	      					//p->iindices[numentries]=p->blocksize;
	      					//p->jindices[numentries]=p->blocksize;
						//numentries++;
					}
					else{
						p->entries[numentries]=1.0;
	      					p->iindices[numentries]=ii;
	      					p->jindices[numentries]=jj;					
						if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,p->iindices[numentries],p->jindices[numentries],p->entries[numentries]);
						numentries++;
						//p->entries[numentries]=0.5;
	      					//p->iindices[numentries]=jj;
	      					//p->jindices[numentries]=ii;
						//numentries++;
						/*p->entries[numentries]=0.0;
	      					p->iindices[numentries]=1;
	      					p->jindices[numentries]=1;
						if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,p->iindices[numentries],p->jindices[numentries],p->entries[numentries]);
						numentries++;*/
					}
				}
				else{
						if(p->blocksize<3){
      							if(DEBUG_FNLSDP && np==0)printf("build_constraints: Too small blocksize!\n");
      							exit(10);
						}
						p->entries[numentries]=0.0;
	      					p->iindices[numentries]=1;
	      					p->jindices[numentries]=1;
						if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,p->iindices[numentries],p->jindices[numentries],p->entries[numentries]);
						numentries++;
						/*p->entries[numentries]=0.0;
	      					p->iindices[numentries]=p->blocksize;
	      					p->jindices[numentries]=p->blocksize;
						if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,p->iindices[numentries],p->jindices[numentries],p->entries[numentries]);
						numentries++;*/
						//p->entries[numentries]=0.0;
	      					//p->iindices[numentries]=2;
	      					//p->jindices[numentries]=2;
						//numentries++;
				}
			}
			
			if(block==4){	
				p->entries[numentries]=1.0;
	      			p->iindices[numentries]=i;
	      			p->jindices[numentries]=nu*ny+nx*(nx+1)+1;
				if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,p->iindices[numentries],p->jindices[numentries],p->entries[numentries]);
				numentries++;
				//p->entries[numentries]=1.0;
	      			//p->iindices[numentries]=nu*ny+2*nx*nx+1;
	      			//p->jindices[numentries]=i;
				//numentries++;			
			}
			
			if(block==2){
					
				for(k=1;k<=4*nx*nx;k++){
					if(i>=1 && i<=nu*ny){
						if(k>=1 && k<=nx*nx){
							ii= ((i-1) % nu)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nu)); 
							kk = ((k-1) % nx)+1;
							kkk = (int)ceil(((double)k) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);
							//printf("i=%d, ii=%d, iii=%d\n",i,ii,iii);
							//printf("(k>=1 && k<=%d)\n",nx*nx);
							//printf("AA[%d]=%f (%d) (%d,%d)\n",k,Ctilde.mat[ijtok(kkk,iii,Ctilde.nrows)]*B->mat[ijtok(kk,ii,B->nrows)]+BT.mat[ijtok(ii,kkk,BT.nrows)],ijtok(kk,iii,CtildeT.nrows),CtildeT.nrows,CtildeT.ncols);
							p->entries[numentries]=Ctilde.mat[ijtok(iii,kkk,Ctilde.nrows)]*B->mat[ijtok(kk,ii,B->nrows)]+BT.mat[ijtok(ii,kkk,BT.nrows)]*CtildeT.mat[ijtok(kk,iii,CtildeT.nrows)];							
						}
						if(k>nx*nx && k<=2*nx*nx){
							kaux=k-nx*nx;
							ii= ((i-1) % nu)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nu)); 
							kk = ((kaux-1) % nx)+1;
							kkk = (int)ceil(((double)kaux) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);
							//printf("(k>=%d && k<=%d)\n",nx*nx,2*nx*nx);
							p->entries[numentries]=-(Ctilde.mat[ijtok(iii,kkk,Ctilde.nrows)]*B->mat[ijtok(kk,ii,B->nrows)]+BT.mat[ijtok(ii,kkk,BT.nrows)]*CtildeT.mat[ijtok(kk,iii,CtildeT.nrows)]);
						}
						if(k>2*nx*nx && k<=3*nx*nx){
							kaux=k-2*nx*nx;
							ii= ((i-1) % nu)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nu)); 
							kk = ((kaux-1) % nx)+1;
							kkk = (int)ceil(((double)kaux) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);
							//printf("(k>=%d && k<=%d)\n",2*nx*nx,3*nx*nx);
							p->entries[numentries]=Ctildep.mat[ijtok(iii,kkk,Ctildep.nrows)]*B->mat[ijtok(kk,ii,B->nrows)]+BT.mat[ijtok(ii,kkk,BT.nrows)]*CtildepT.mat[ijtok(kk,iii,CtildepT.nrows)];
						}
						if(k>3*nx*nx && k<=4*nx*nx){
							kaux=k-3*nx*nx;
							ii= ((i-1) % nu)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nu)); 
							kk = ((kaux-1) % nx)+1;
							kkk = (int)ceil(((double)kaux) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);							
							//printf("(k>=%d && k<=%d)\n",3*nx*nx,4*nx*nx);
							p->entries[numentries]=-(Ctildep.mat[ijtok(iii,kkk,Ctildep.nrows)]*B->mat[ijtok(kk,ii,B->nrows)]+BT.mat[ijtok(ii,kkk,BT.nrows)]*CtildepT.mat[ijtok(kk,iii,CtildepT.nrows)]);
						}	
					}
					if(i>nu*ny && i<=(int)(nx*(nx+1)/2)+nu*ny){
						if(k>=1 && k<=nx*nx){
							ii= ((i-1) % nx)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nx)); 
							kk = ((k-1) % nx)+1;
							kkk = (int)ceil(((double)k) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);
							//printf("i=%d, ii=%d, iii=%d\n",i,ii,iii);
							//printf("(k>=1 && k<=%d)\n",nx*nx);
							//printf("AA[%d]=%f (%d) (%d,%d)\n",k,Ctilde.mat[ijtok(kkk,iii,Ctilde.nrows)]*B->mat[ijtok(kk,ii,B->nrows)]+BT.mat[ijtok(ii,kkk,BT.nrows)],ijtok(kk,iii,CtildeT.nrows),CtildeT.nrows,CtildeT.ncols);
						
							if(DEBUG_FNLSDP && np==0)printf("TT3.mat[ijtok(%d,%d,nx*nx)]=%f\n",k,i-nu*ny,TT3.mat[ijtok(k,i-nu*ny,nx*nx)]);	
							p->entries[numentries]=TT3.mat[ijtok(k,i-nu*ny,nx*nx)];
							/*if(kk==ii)
								p->entries[numentries]=AF->mat[ijtok(kk,ii,AF->nrows)]+AF->mat[ijtok(iii,kkk,AF->nrows)];
							else
								p->entries[numentries]=AF->mat[ijtok(kk,ii,AF->nrows)];	
							*/
						}
						if(k>nx*nx && k<=2*nx*nx){
							kaux=k-nx*nx;
							ii= ((i-1) % nx)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nx)); 
							kk = ((kaux-1) % nx)+1;
							kkk = (int)ceil(((double)kaux) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);
							//printf("(k>=%d && k<=%d)\n",nx*nx,2*nx*nx);
							
							p->entries[numentries]=-TT3.mat[ijtok(k-nx*nx,i-nu*ny,nx*nx)];
							/*if(kk==ii)
								p->entries[numentries]=-(AF->mat[ijtok(kk,ii,AF->nrows)]+AF->mat[ijtok(iii,kkk,AF->nrows)]);
							else
								p->entries[numentries]=-AF->mat[ijtok(kk,ii,AF->nrows)];
							*/
						}
						if(k>2*nx*nx && k<=3*nx*nx){
							kaux=k-2*nx*nx;
							ii= ((i-1) % nx)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nx)); 
							kk = ((kaux-1) % nx)+1;
							kkk = (int)ceil(((double)kaux) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);
							//printf("(k>=%d && k<=%d)\n",2*nx*nx,3*nx*nx);
							p->entries[numentries]=0.0;
						}
						if(k>3*nx*nx && k<=4*nx*nx){
							kaux=k-3*nx*nx;
							ii= ((i-1) % nx)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nx)); 
							kk = ((kaux-1) % nx)+1;
							kkk = (int)ceil(((double)kaux) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);							
							//printf("(k>=%d && k<=%d)\n",3*nx*nx,4*nx*nx);
							p->entries[numentries]=0.0;
						}	
					}
					if(i>(int)(nx*(nx+1)/2)+nu*ny && i<=nx*(nx+1)+nu*ny){
						if(k>=1 && k<=nx*nx){
							ii= ((i-1) % nx)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nx)); 
							kk = ((k-1) % nx)+1;
							kkk = (int)ceil(((double)k) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);
							//printf("i=%d, ii=%d, iii=%d\n",i,ii,iii);
							//printf("(k>=1 && k<=%d)\n",nx*nx);
							//printf("AA[%d]=%f (%d) (%d,%d)\n",k,Ctilde.mat[ijtok(kkk,iii,Ctilde.nrows)]*B->mat[ijtok(kk,ii,B->nrows)]+BT.mat[ijtok(ii,kkk,BT.nrows)],ijtok(kk,iii,CtildeT.nrows),CtildeT.nrows,CtildeT.ncols);
							p->entries[numentries]=0.0;							
						}
						if(k>nx*nx && k<=2*nx*nx){
							kaux=k-nx*nx;
							ii= ((i-1) % nx)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nx)); 
							kk = ((kaux-1) % nx)+1;
							kkk = (int)ceil(((double)kaux) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);
							//printf("(k>=%d && k<=%d)\n",nx*nx,2*nx*nx);
							p->entries[numentries]=0.0;
						}
						if(k>2*nx*nx && k<=3*nx*nx){
							kaux=k-2*nx*nx;
							ii= ((i-1) % nx)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nx)); 
							kk = ((kaux-1) % nx)+1;
							kkk = (int)ceil(((double)kaux) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);
							//printf("(k>=%d && k<=%d)\n",2*nx*nx,3*nx*nx);
							
							p->entries[numentries]=TT3.mat[ijtok(k-2*nx*nx,i-(int)(nx*(nx+1)/2)-nu*ny,nx*nx)];

							/*if(kk==ii)
								p->entries[numentries]=AF->mat[ijtok(kk,ii,AF->nrows)]+AF->mat[ijtok(iii,kkk,AF->nrows)];
							else
								p->entries[numentries]=AF->mat[ijtok(kk,ii,AF->nrows)];*/
						}
						if(k>3*nx*nx && k<=4*nx*nx){
							kaux=k-3*nx*nx;
							ii= ((i-1) % nx)+1; // intern col
							iii= (int)ceil(((double)i) /((double)nx)); 
							kk = ((kaux-1) % nx)+1;
							kkk = (int)ceil(((double)kaux) /((double)nx));
							//printf("k=%d, kk=%d, kkk=%d\t",k,kk,kkk);							
							//printf("(k>=%d && k<=%d)\n",3*nx*nx,4*nx*nx);
							
							p->entries[numentries]=-TT3.mat[ijtok(k-3*nx*nx,i-(int)(nx*(nx+1)/2)-nu*ny,nx*nx)];

							/*if(kk==ii)
								p->entries[numentries]=-(AF->mat[ijtok(kk,ii,AF->nrows)]+AF->mat[ijtok(iii,kkk,AF->nrows)]);
							else
								p->entries[numentries]=-AF->mat[ijtok(kk,ii,AF->nrows)];*/
						}	
					}	      				
					p->iindices[numentries]=k;
	      				p->jindices[numentries]=k;
					if(DEBUG_FNLSDP && np==0)printf("%d %d %d %d %f\n",i,block,k,k,p->entries[numentries]);
					numentries++;				
				}
			}
		
			oldp=p;	
	  		p=p->next;
			block++;
		}
    	}
	free(gamma);
	free_mat_gen(&Ctilde,np);
	free_mat_gen(&Ctildep,np);
	free_mat_gen(&CtildeT,np);
	free_mat_gen(&CtildepT,np);
	free_mat_gen(&BT,np);
	free_mat_gen(&CT,np);
	free_mat_gen(&CFT,np);
	free_mat_gen(&AUX1,np);
	free_mat_gen(&AUX2,np);
	free_mat_gen(&AUX3,np);
	free_mat_gen(&AUX4,np);
	free_mat_gen(&index,np);
	free_mat_gen(&index2,np);
	free_mat_gen(&TTT1,np);
	free_mat_gen(&TTT2,np);
	free_mat_gen(&TT3,np);
}

void free_constraints(nx,ny,nu,constraints,np)
     int nx;
     int ny;
     int nu;
     struct constraintmatrix *constraints;
     int np;
{
  int i;
  struct sparseblock *ptr;
  struct sparseblock *oldptr;
  int k=nu*ny+nx*(nx+1)+1;
  if (constraints != NULL)
    {
      for (i=1; i<=k; i++)
	{
	  /*
	   * Get rid of constraint i.
	   */
	  
	  ptr=constraints[i].blocks;
	  while (ptr != NULL)
	    {
	      free(ptr->entries);
	      free(ptr->iindices);
	      free(ptr->jindices);
	      oldptr=ptr;
	      ptr=ptr->next;
	      free(oldptr);
	    };
	};
      /*
       * Finally, free the constraints array.
       */

      free(constraints);
    };

}


void free_mat_gen(genmatrix *A,int np){	
	free(A->mat);
}

void vec2mats(double *X, genmatrix *F0, genmatrix *Q0, genmatrix *V0,int nu, int ny, int nx){
	int i,j;

	F0->nrows=nu;
	F0->ncols=ny;
	Q0->nrows=nx;
	Q0->ncols=nx;
	V0->nrows=nx;
	V0->ncols=nx;

	for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			F0->mat[i+j*nu]=X[i+j*nu];
	}

	int count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			Q0->mat[i+j*nx]=X[count+nu*ny];
			if(i!=j)
				Q0->mat[j+i*nx]=X[count+nu*ny];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			V0->mat[i+j*nx]=X[count+nu*ny+nx*(nx+1)/2];
			if(i!=j)
				V0->mat[j+i*nx]=X[count+nu*ny+nx*(nx+1)/2];
			count++;
		}
	}	
}

void mats2vec(double *X, genmatrix *F0, genmatrix *Q0, genmatrix *V0,int nu, int ny, int nx){
	int i,j;
	
	//for(i=0;i<nu*ny+nx*(nx+1);i++)
	//	printf("X[%d]=%f\n",i,X[i]);

	for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			X[i+j*nu]=F0->mat[i+j*nu];
	}

	int count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			X[count+nu*ny]=Q0->mat[i+j*nx];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			X[count+nu*ny+nx*(nx+1)/2]=V0->mat[i+j*nx];
			count++;
		}
	}
}
