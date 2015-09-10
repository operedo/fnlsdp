/*
 *  This is an easy to call version of the sdp routine.  It takes as
 *  input a problem (n,k,C,a,constraints,constant_offset), and an 
 *  initial solution (X,y,Z), allocates working storage, and calls sdp() 
 *  to solve the problem.  The solution is returned in X,y,Z,pobj,dobj, and 
 *  the return code from sdp is returned as the return value from easy_sdp.
 *
 */

/* This file was modified compared to version 5.0.
 *  It is not meant to be compiled without ScaLAPACK and BLACS.
 *  
 *  Modified by Ivan Ivanov  31 Jan 2007
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "declarations.h"

int easy_sdp(n,k,C,a,constraints,constant_offset,pX,py,pZ,ppobj,pdobj,scapack,params,printlevel)
     int n;
     int k;
     struct blockmatrix C;
     double *a;
     struct constraintmatrix *constraints;
     double constant_offset;
     struct blockmatrix *pX;
     double **py;
     struct blockmatrix *pZ;
     double *ppobj;
     double *pdobj;
     struct scalapackpar scapack;
     struct paramstruc params;
     int printlevel;
{
  int ret;
  struct constraintmatrix fill;
  struct blockmatrix work1;
  struct blockmatrix work2;
  struct blockmatrix work3;
  struct blockmatrix bestx;
  struct blockmatrix bestz;
  struct blockmatrix Zi;
  struct blockmatrix dZ;
  struct blockmatrix dX;
  struct blockmatrix cholxinv;
  struct blockmatrix cholzinv;
  double *workvec1;
  double *workvec2;
  double *workvec3;
  double *workvec4;
  double *workvec5;
  double *workvec6;
  double *workvec7;
  double *workvec8;
  double *diagO; 
  double *Fp;
  double *O;
  double *dy;
  double *dy1;
  double *rhs;
  double *besty;
  int ldam;
  struct sparseblock **byblocks;
  struct sparseblock *ptr;
  struct sparseblock *oldptr;
  int i;
  int j;
  int blk;
  struct sparseblock *p;
  struct sparseblock *q;
  struct sparseblock *prev=NULL;
  double gap;
  int nnz;

  int nfree=0;
  int *freevarblocks=NULL;
  int *freevarind=NULL;
  int block;

  int descOd[9],izero=0,ione=1,itemp,nru,info=0,nrhs;
  double *O_d, *mind_tmp, *mind_gen, *RHS_d;

  
  /*
   *  Allocate working storage
   */

  alloc_mat(C,&work1);
  alloc_mat(C,&work2);
  alloc_mat(C,&work3);
  alloc_mat_packed(C,&bestx);
  alloc_mat_packed(C,&bestz);
  alloc_mat_packed(C,&cholxinv);
  alloc_mat_packed(C,&cholzinv);

  besty=(double *)malloc(sizeof(double)*(k+1));
   if (besty == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   if (n > k)
     {
       workvec1=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec1=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec1 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   if (n > k)
     {
       workvec2=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec2=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec2 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   if (n > k)
     {
       workvec3=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec3=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec3 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   if (n > k)
     {
       workvec4=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec4=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec4 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   if (n > k)
     {
       workvec5=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec5=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec5 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   if (n > k)
     {
       workvec6=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec6=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec6 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   if (n > k)
     {
       workvec7=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec7=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec7 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   if (n > k)
     {
       workvec8=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec8=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec8 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };


   if (n > k)
     {
       diagO=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       diagO=(double *)malloc(sizeof(double)*(k+1));
     };
   if (diagO == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };



   rhs=malloc(sizeof(double)*(k+1));
   if (rhs == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   dy=malloc(sizeof(double)*(k+1));
   if (dy == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   dy1=malloc(sizeof(double)*(k+1));
   if (dy1 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   Fp=malloc(sizeof(double)*(k+1));
   if (Fp == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   /*
    *  Work out the leading dimension for the array.  Note that we may not
    *  want to use k itself, for cache issues.
    */
   if ((k % 2) == 0)
     ldam=k+1;
   else
     ldam=k;


  /* initialize parameters needed for  scalapack lib */
   
   scapack.nb=64;
  

   if (scapack.nb>ldam){ 
       scapack.nb = (int) (ldam/scapack.nprow); 
    }
   else{
     i= (int) (ldam/scapack.nprow); 
     if (i < scapack.nb)
        scapack.nb = i;  
   }


   scapack.lr = numroc_(&ldam, &scapack.nb, &scapack.myrow,&izero,&scapack.nprow);
   scapack.lc = numroc_(&ldam, &scapack.nb, &scapack.mycol,&izero,&scapack.npcol);
   nrhs = numroc_(&ione, &scapack.nb, &scapack.mycol, &izero, &scapack.npcol);


   O = (double *)malloc(scapack.lr*scapack.lc*sizeof(double));
   if (O==NULL) {printf("Memory allocation error of distributed matrix O on proc %dx%d\n",scapack.myrow,scapack.mycol); exit(10);}

    RHS_d = (double *)calloc(scapack.lr*nrhs,sizeof(double));
   if (RHS_d==NULL){printf("Memory allocation error of distributed vector RHS_d on proc %dx%d\n",scapack.myrow,scapack.mycol); exit(10);}

   O_d = (double *)calloc(scapack.lr*scapack.lc,sizeof(double));
   if (O_d==NULL) {printf("Memory allocation error of distributed matrix O_d on proc %dx%d\n",scapack.myrow,scapack.mycol); exit(10);}

   mind_tmp = (double*)malloc(scapack.np*sizeof(double));
   if (mind_tmp==NULL) {printf("Memory allocation error of distributed vector min_tmp on proc %dx%d\n",scapack.myrow,scapack.mycol); exit(10);}

    mind_gen = (double*)malloc(scapack.np*sizeof(double));
   if (mind_gen==NULL) {printf("Memory allocation error of distributed vector min_tmp on proc %dx%d\n",scapack.myrow,scapack.mycol); exit(10);}

   nru = numroc_(&ldam, &scapack.nb, &scapack.myrow, &izero, &scapack.nprow);
   itemp = max(1,nru);
   descinit_(descOd, &ldam, &ldam, &scapack.nb, &scapack.nb, &izero, &izero, &scapack.ic,&itemp, &info);



   alloc_mat(C,&Zi);
   alloc_mat(C,&dZ);
   alloc_mat(C,&dX);

   /*
    *  Fill in lots of details in the constraints data structure that haven't
    *  necessarily been done before now.
    */

   /*
    * Set up the cross links used by op_o
    * While we're at it, determine which blocks are sparse and dense.
    */

   /*
    * Next, setup issparse and NULL out all nextbyblock pointers.
    */

   for (i=1; i<=k; i++)
     {
       p=constraints[i].blocks;
       while (p != NULL)
	 {
	   /*
	    * First, set issparse.
	    */
	   if (((p->numentries) > 0.25*(p->blocksize)) && ((p->numentries) > 15))
	     {
	       p->issparse=0;
	     }
	   else
	     {
	       p->issparse=1;
	     };
	   
	   if (C.blocks[p->blocknum].blockcategory == DIAG)
	     p->issparse=1;
	   
	   /*
	    * Setup the cross links.
	    */
	   
	   p->nextbyblock=NULL;
	   p=p->next;
	 };
     };
   
   /*
    * Now, cross link.
    */
   for (i=1; i<=k; i++)
     {
       p=constraints[i].blocks;
       while (p != NULL)
	 {
	   if (p->nextbyblock == NULL)
	     {
	       blk=p->blocknum;
	       
	       /*
		* link in the remaining blocks.
		*/
	       for (j=i+1; j<=k; j++)
		 {
		   q=constraints[j].blocks;
		   
		   while (q != NULL)
		     {
		       if (q->blocknum == p->blocknum)
			 {
			   if (p->nextbyblock == NULL)
			     {
			       p->nextbyblock=q;
			       q->nextbyblock=NULL;
			       prev=q;
			     }
			   else
			     {
			       prev->nextbyblock=q;
			       q->nextbyblock=NULL;
			       prev=q;
			     };
			   break;
			 };
		       q=q->next;
		     };
		 };
	     };
	   p=p->next;
	 };
     };

   /*
    * If necessary, print out information on sparsity of blocks.
    */
   
   if ((printlevel >= 4)&&(scapack.id==0))
     {
       for (i=1; i<=k; i++)
	 {
	   p=constraints[i].blocks;
	   while (p != NULL)
	     {
	       printf("%d,%d,%d,%d \n",i,p->blocknum,p->issparse,p->numentries);
	       p=p->next;
	     };
	 };
     };
   
   /*
    * Allocate space for byblocks pointers.
    */

   byblocks=(struct sparseblock **)malloc((C.nblocks+1)*sizeof(struct sparseblock *));
   if ((byblocks == NULL)&&(scapack.id==0))
     {
       printf("Storage Allocation Failed!\n");
       exit(10);
     };

   for (i=1; i<=C.nblocks; i++)
     byblocks[i]=NULL;

   /*
    * Fill in byblocks pointers.
    */
   for (i=1; i<=k; i++)
     {
       ptr=constraints[i].blocks;
       while (ptr != NULL)
	 {
	   if (byblocks[ptr->blocknum]==NULL)
	     {
	       byblocks[ptr->blocknum]=ptr;
	     };
	   ptr=ptr->next;
	 };
     };

   /*
    *  Compute "fill".  This data structure tells us which elements in the
    *  block diagonal matrix have corresponding elements in one of the
    *  constraints, and which constraint this element first appears in.
    *
    */

   makefill(k,C,constraints,&fill,work1,printlevel,&params.rank1, scapack.id);

   /*
    * Compute the nonzero structure of O.
    */

   nnz=structnnz(n,k,C,constraints);

  
   /*
    * Identify any split free variables in the problem formulation.
    */

   for (block=1; block<=C.nblocks; block++)
     {
       if (C.blocks[block].blockcategory==DIAG)
	 {
	   j=1;
	   while (j<= C.blocks[block].blocksize-1)
	     {
	       if (fabs(C.blocks[block].data.vec[j]+C.blocks[block].data.vec[j+1]) <= 1.0e-8*(1+fabs(C.blocks[block].data.vec[j])))
		 {
		   zero_mat(work1);
		   work1.blocks[block].data.vec[j]=1e6;
		   work1.blocks[block].data.vec[j+1]=1e6;
		   op_a(k,constraints,work1,workvec1,&params.rank1);
		   if (norm2(k,workvec1+1) < 0.01 )
		     {
		       /*
			* This is a split free variable pair.
			*/
		       if ((printlevel >= 3)&&(scapack.id==0))
			 {
			   printf("Detected free variable pair: %d, %d \n",block,j);
			 };
		       nfree++;
		       if (nfree==1)
			 {
			   freevarind=(int *)malloc((n+1)*sizeof(int));
			   freevarblocks=(int *)malloc((n+1)*sizeof(int));
			   if (freevarblocks==NULL)
			     {
			       if (printlevel >= 1)
				 printf("Couldn't allocate storage on node %d!\n", scapack.id);
			       exit(1000);
			     };
			 };
		       freevarind[nfree]=j;
		       freevarblocks[nfree]=block;
		       j=j+1;
		     };
		 };
	       j=j+1;
	     };
	 };

     };

   /*
    * Sort entries in diagonal blocks of constraints.
    */

   sort_entries(k,C,constraints);


   /*
    *  Now, call sdp().
    */
 
   ret=sdp(n,k,C,a,constant_offset,constraints,byblocks,fill,*pX,*py,*pZ,
	   cholxinv,cholzinv,ppobj,pdobj,work1,work2,work3,workvec1,
	   workvec2,workvec3,workvec4,workvec5,workvec6,workvec7,workvec8,
	   diagO,bestx,besty,bestz,Zi,O,O_d,RHS_d,descOd,mind_tmp,mind_gen,rhs,dZ,dX,
           dy,dy1,Fp,nfree,freevarblocks,freevarind,printlevel,params,scapack);


   if ((printlevel >= 1)&&(scapack.id==0))
     {
       if (ret==0)
	 printf(" %d Success: SDP solved\n",scapack.id);
       if (ret==1)
	 printf("Success: SDP is primal infeasible\n");
       if (ret==2)
	 printf("Success: SDP is dual infeasible\n");
       if (ret==3)
	 printf("Partial Success: SDP solved with reduced accuracy\n");
       if (ret >= 4)
	 printf("Failure: return code is %d \n",ret);

       if (ret==1)
	 {
	   op_at(k,*py,constraints,work1,&params.rank1);
	   addscaledmat(work1,-1.0,*pZ,work1);
	   printf("Certificate of primal infeasibility: a'*y=%.5e, ||A'(y)-Z||=%.5e\n",-1.0,Fnorm(work1));
	 };

       if (ret==2)
	 {
	   op_a(k,constraints,*pX,workvec1,&params.rank1);
	   printf("Certificate of dual infeasibility: tr(CX)=%.5e, ||A(X)||=%.5e\n",trace_prod(C,*pX),norm2(k,workvec1+1));
	 };

       if ((ret==0) || (ret>=3))
	 {
	   if (printlevel >= 3)
	     {
	       printf("XZ Gap: %.7e \n",trace_prod(*pZ,*pX));
	       gap=*pdobj-*ppobj;
	       printf("Real Gap: %.7e \n",gap);
	     };

	   if (printlevel >= 1)
	     {
               gap=*pdobj-*ppobj;

	       printf(" Primal objective value: %.7e \n",*ppobj);
	       printf(" Dual objective value: %.7e \n",*pdobj);
	       printf(" Relative primal infeasibility: %.2e \n",
		      pinfeas(k,constraints,*pX,a,workvec1,&params.rank1));
	       printf(" Relative dual infeasibility: %.2e \n",
		      dinfeas(k,C,constraints,*py,*pZ,work1,&params.rank1));
	       printf(" Real Relative Gap: %.2e \n",gap/(1+fabs(*pdobj)+fabs(*ppobj)));
	       printf(" XZ Relative Gap: %.2e \n",trace_prod(*pZ,*pX)/(1+fabs(*pdobj)+fabs(*ppobj)));

	       printf(" DIMACS error measures: %.2e %.2e %.2e %.2e %.2e %.2e\n",
		      pinfeas(k,constraints,*pX,a,workvec1,&params.rank1)*(1+norm2(k,a+1))/
		      (1+norminf(k,a+1)),
		      0.0,
		      dimacserr3(k,C,constraints,*py,*pZ,work1,&params.rank1),
		      0.0,
		      gap/(1+fabs(*pdobj)+fabs(*ppobj)),
		      trace_prod(*pZ,*pX)/(1+fabs(*pdobj)+fabs(*ppobj)));
             
	     };
	 };

     };



   /*
    *  Now, free up all of the storage.
    */

   free_mat(work1);
   free_mat(work2);
   free_mat(work3);
   free_mat_packed(bestx);
   free_mat_packed(bestz);
   free_mat_packed(cholxinv);
   free_mat_packed(cholzinv);

   free_mat(Zi);
   free_mat(dZ);
   free_mat(dX);

   free(besty);
   free(workvec1);
   free(workvec2);
   free(workvec3);
   free(workvec4);
   free(workvec5);
   free(workvec6);
   free(workvec7);
   free(workvec8);
   free(rhs);
   free(dy);
   free(dy1);
   free(Fp);
   free(O);
   free(O_d);
   free(RHS_d);
   free(mind_tmp);
   free(mind_gen);
   free(diagO);
   free(byblocks);
   if (nfree > 0)
     {
       free(freevarblocks);
       free(freevarind);
     };

 
   /*
    * Free up the fill data structure.
    */

   ptr=fill.blocks;
   while (ptr != NULL)
     {
       free(ptr->entries);
       free(ptr->iindices);
       free(ptr->jindices);
       oldptr=ptr;
       ptr=ptr->next;
       free(oldptr);
     };


  /*
   * Finally, free the constraints array.
   */
   
  
   return(ret);

}

int structnnz(n,k,C,constraints)
     int n;
     int k;
     struct blockmatrix C;
     struct constraintmatrix *constraints;
{
  int i,j;
  int ii,jj;
  int nnz;
  struct sparseblock *ptri;
  struct sparseblock *ptrj;
  

  nnz=0;
  for (i=1; i<=k; i++)
    for (j=1; j<=k; j++)
      {
	ptri=constraints[i].blocks;

	while (ptri != NULL)
	  {

	    ptrj=constraints[j].blocks;
	    while (ptrj != NULL)
	      {

		if (ptri->blocknum == ptrj->blocknum)
		  {
		    if (C.blocks[ptri->blocknum].blockcategory==MATRIX)
		      {
			nnz++;
			goto NEXTJ;
		      }
		    else
		      { /* DIAG block */
			for (ii=1; ii<=ptri->numentries; ii++)
			  for (jj=1; jj<=ptrj->numentries; jj++)
			    {
			      if (ptri->iindices[ii]==ptrj->iindices[jj])
				{
				  nnz++;
				  goto NEXTJ;
				};
			    };
		      };
		  };
		ptrj=ptrj->next;
	      };

	    ptri=ptri->next;
	  }; /* end while */

      NEXTJ:;
      }; /* end nested fors */
 
  return(nnz);
}

int actnnz(n,lda,A)
     int n;
     int lda;
     double A[];
{
  int i,j;
  int nnz;
 
  nnz=0;
  for (i=1; i<=n; i++)
    {
      if (A[ijtok(i,i,lda)] != 0.0)
	nnz++;
      for (j=i+1; j<=n; j++)
	{
	  if (A[ijtok(i,j,lda)] != 0.0)
	    {
	      nnz++;
	      nnz++;
	    };
	};
    };

  return(nnz);
}

int bandwidth(n,lda,A)
     int n;
     int lda;
     double A[];
{
  int i;
  int j;
  int bw;

  bw=0;
  for (j=2; j<=n; j++)
    {
      for (i=1; i<=j-1; i++)
	{
	  if (A[ijtok(i,j,lda)] != 0.0)
	    {
	      if ((j-i) > bw)
		bw=j-i;
	      break;
	    };
	};
    };

  return(bw);
}

