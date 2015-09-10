/*
  Add a matrix plus a multiple of a second matrix and put the result in a
  third matrix.

    C=A+scale*B

  */

#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"

void addscaledmat(A,scale,B,C)
     struct blockmatrix A;
     double scale;
     struct blockmatrix B;
     struct blockmatrix C;
{
  int blk;
  int i;

  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory)
	{
	case DIAG:

	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    C.blocks[blk].data.vec[i] = A.blocks[blk].data.vec[i] + scale*B.blocks[blk].data.vec[i];
	  break;
	case MATRIX:
	  for (i=0; i<=A.blocks[blk].blocksize*A.blocks[blk].blocksize-1; i++)
	    C.blocks[blk].data.vec[i] = A.blocks[blk].data.vec[i] + scale*B.blocks[blk].data.vec[i];
	  break;
	case PACKEDMATRIX:
	default:
	  printf("addscaledmat illegal block type \n");
	  exit(12);
	};
    };


}


