/*
  Make an identity matrix.
 */

#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"

void make_i(A)
     struct blockmatrix A;
{
  int blk,i;
  double *p;

  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory) 
	{
	case DIAG:
	  p=A.blocks[blk].data.vec;
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    p[i]=1.0;
	  break;
	case MATRIX:
	  p=A.blocks[blk].data.mat;
	  for (i=0; i<=A.blocks[blk].blocksize*A.blocks[blk].blocksize-1; i++)
	    p[i]=0.0;
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    p[ijtok(i,i,A.blocks[blk].blocksize)]=1.0;
	  break;
	default:
	  printf("make_i illegal block type\n");
	  exit(12);
	};
    }
}






