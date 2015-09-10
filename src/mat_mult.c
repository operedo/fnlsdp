/*
 *  Compute C=scale1*A*B+scale2*C.  
 *  Note that C must consist of dense matrix and vector blocks- no sparse
 *  blocks or eye's or other special cases.
 *
 *  A and B can have blocks of all supported types.  Unsupported types
 *  generate exit(1).
 *
 *  It is assumed that all three matrices are of compatible block strucutre.
 *
 *  We use dgemm to do the actual work on normal dense blocks.
 *
 */

/* This file was modified compared to version 5.0.
*  Two functions are new:  vec_mult_mat_raw() and
*  vec_mult_vec(). They are part of op_o().
*  
*  Modified by Ivan Ivanov  31 Jan 2007
*/


#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"

void mat_mult(scale1,scale2,A,B,C)
     double scale1,scale2;
     struct blockmatrix A,B,C;
{
  int blk,i,n;
  double *ap;
  double *bp;
  double *cp;

  /*
   * In theory, the BLAS ensures that if scale2=0, then C will not be 
   * accessed before being written to.  In practice, this is not always 
   * true, so we initilize C to zeros for safety.
   */

  if (scale2 == 0.0)
    zero_mat(C);


  /*
   * Work through the blocks one at a time.
   */
  
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory) 
	{
	case DIAG:
          if (scale2 != 0.0)
	    {
	      for (i=1; i<=A.blocks[blk].blocksize; i++)
		{
		  C.blocks[blk].data.vec[i]=scale1*A.blocks[blk].data.vec[i]
		    *B.blocks[blk].data.vec[i]
		    +scale2*C.blocks[blk].data.vec[i];
		};
	    }
	  else
	    {
	      for (i=1; i<=A.blocks[blk].blocksize; i++)
		{
		  C.blocks[blk].data.vec[i]=scale1*A.blocks[blk].data.vec[i]
		    *B.blocks[blk].data.vec[i];
		};
	    };
	  break;
	case MATRIX:
	  /*
	   * Call dgemm to do the matrix multiplication.
	   */
	  
	  n=A.blocks[blk].blocksize;
	  ap=A.blocks[blk].data.mat;
	  bp=B.blocks[blk].data.mat;
	  cp=C.blocks[blk].data.mat;

	  mat_mult_raw(n,scale1,scale2,ap,bp,cp);
	  break;
	default:
	  printf("mat_mult illegal block type!\n");
	  exit(12);
	};

    };
}

void mat_mult_raw(n,scale1,scale2,ap,bp,cp)
     int n;
     double scale1;
     double scale2;
     double *ap;
     double *bp;
     double *cp;
{
#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
	  DGEMM("N","N",&n,&n,&n,&scale1,ap,&n,bp,&n,&scale2,cp,&n);
#else
	  dgemm("N","N",&n,&n,&n,&scale1,ap,&n,bp,&n,&scale2,cp,&n);
#endif
#else
#ifdef CAPSBLAS
	  DGEMM_("N","N",&n,&n,&n,&scale1,ap,&n,bp,&n,&scale2,cp,&n);
#else
	  dgemm_("N","N",&n,&n,&n,&scale1,ap,&n,bp,&n,&scale2,cp,&n);
#endif
#endif

}

void vec_mult_mat_raw(n,scale1,scale2,ap,bp,cp)
     int n;
     double scale1;
     double scale2;
     double *ap;
     double *bp;
     double *cp;
{
int ione=1;

#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
          DGEMV("T",&n,&n,&scale1,ap,&n,bp,&ione,&scale2,cp,&ione);
#else
          dgemv("T",&n,&n,&scale1,ap,&n,bp,&ione,&scale2,cp,&ione);
#endif
#else
#ifdef CAPSBLAS
          DGEMV_("T",&n,&n,&scale1,ap,&n,bp,&ione,&scale2,cp,&ione);
#else
          dgemv_("T",&n,&n,&scale1,ap,&n,bp,&ione,&scale2,cp,&ione);
#endif
#endif

}

double vec_mult_vec(n,bp,cp)
     int n;
     double *bp;
     double *cp;
{
int ione=1;
double res=0.0;

#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
          res = DDOT(&n,bp,&ione,cp,&ione);
#else
          res = ddot(&n,bp,&ione,cp,&ione);
#endif
#else
#ifdef CAPSBLAS
          res = DDOT_(&n,bp,&ione,cp,&ione);
#else
          res = ddot_(&n,bp,&ione,cp,&ione);
#endif
#endif

return(res);
}


