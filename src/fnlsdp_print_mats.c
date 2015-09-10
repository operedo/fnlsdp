/**
 * \file fnlsdp_print_mats.c 
 * \author Oscar Peredo
 * \date 09 Feb 2010
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "declarations.h"
//#include "fnlsdp.h"


/**
 * \brief Print general matrix
 * \param A
 */
void print_genmatrix(genmatrix *A){
		int i=0,j=0;
		for(i=0;i<A->nrows;i++){
			for(j=0;j<A->ncols;j++)
				printf("%f\t",A->mat[i+j*A->nrows]);
			printf("\n");
		}
}
/**
 * \brief Fix general matrix (set to 0 values less than 10e-8)
 * \param A
 */
void fix_genmatrix(genmatrix *A){
		int i=0,j=0;
		for(i=0;i<A->nrows;i++){
			for(j=0;j<A->ncols;j++){
				if(A->mat[i+j*A->nrows]<=1e-8 && A->mat[i+j*A->nrows]>=-1e-8){
					A->mat[i+j*A->nrows]=0.0;
				}
			}
		}
}


/**
 * \brief Print block matrix
 * \param A
 */
void print_blockmatrix(struct blockmatrix *A){
  	int blk,i,j;
  	double *p;

  	for (blk=1; blk<=A->nblocks; blk++)
    	{
      		printf("block %d:\n",blk);
		switch (A->blocks[blk].blockcategory) 
		{
			case DIAG:
	  			p=A->blocks[blk].data.vec;
	  			for (i=1; i<=A->blocks[blk].blocksize; i++)
	    				printf("%f\n",A->blocks[blk].data.vec[i]);
	  			break;
			case MATRIX:
	  			p=A->blocks[blk].data.mat;
	  			for (i=0; i<A->blocks[blk].blocksize; i++){
					for(j=0;j<A->blocks[blk].blocksize;j++)
						printf("%f\t",A->blocks[blk].data.mat[i+j*(A->blocks[blk].blocksize)]);
					printf("\n");
				}
	  			break;
			default:
	  			printf("print_blockmatrix: illegal block type\n");
	  			exit(12);
		};
    	}
}

/**
 * \brief  Print objective function of linear SDP
 * \param nx
 * \param ny
 * \param nu
 * \param a
 */
void print_a(int nx, int ny, int nu, double *a){
	int i;
	printf("Objective function:\n");
	for(i=0;i<nu*ny+nx*(nx+1)+1;i++){
		printf("a[%d]=%f\n",i,a[i]);
	}
}
/**
 * \brief Print constraint matrices of linear SDP
 * \param nx
 * \param ny
 * \param nu
 * \param constraints
 */
void print_constraintmatrix(int nx, int ny, int nu, struct constraintmatrix *constraints){
  	int i,j,k;
  	struct sparseblock *p;
  	struct sparseblock *oldptr;

	printf("Constraints:\n");
  	k=nu*ny+nx*(nx+1)+1;
	int block;
/*  	if (constraints != NULL)
    	{
      		for (i=1; i<=k; i++)
		{
			printf("Constraint %d:\n",i);
	  		ptr=constraints[i].blocks;
			//printf("num of blocks: %d\n",ptr->blocknum);
	  		block=1;
			while (ptr != NULL)
	    		{
	      			printf("block %d:\n",block);
				printf("blocksize: %d\n",ptr->blocksize);
				for(j=1;j<=ptr->numentries;j++){
					printf("entry[%d][%d]=%f\n",ptr->iindices[j],ptr->jindices[j],ptr->entries[j]);
				}
	      			ptr=ptr->next;
	    			block++;
			};
		};
    	};
*/
  	for (i=1; i<=k; i++)
    	{
      		p=constraints[i].blocks;
      		while (p != NULL)
		{
	  		for (j=1; j<=p->numentries; j++)
	    		{
	      			printf("%d %d %d %d %.18e \n",i,p->blocknum,
		      		p->iindices[j],
		      		p->jindices[j],
		      		p->entries[j]);
	    		};
	  		p=p->next;
		};
    	};
}
