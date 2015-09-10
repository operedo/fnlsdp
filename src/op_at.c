/*
  Compute A'(y).
  */


/* This file was modified compared to version 5.0.
*  Changes conrern the use of parameter rank1 with 
*  nonzero value.
*  
*  Modified by Ivan Ivanov  31 Jan 2007
*/


#include <stdlib.h>
#include "declarations.h"


void op_at(k,y,constraints,result,prank1)
     int k;
     double *y;
     struct constraintmatrix *constraints;
     struct blockmatrix result;
     int *prank1;
{
  int i,j,jj;
#ifndef BIT64
  int p,q;
#else
  long int p,q;
#endif
  int blk;
  double ent,entex;
  struct sparseblock *ptr;

  zero_mat(result);

  for (i=1; i<=k; i++)
    {
      if (y[i] == 0.0)
	{
	  continue;
	};

      ptr=constraints[i].blocks;

      while (ptr != NULL)
	{
	  blk=ptr->blocknum;

	  if (result.blocks[blk].blockcategory  == DIAG)
	    {
	      for (j=1; j<=ptr->numentries; j++)
		{
		  ent=ptr->entries[j];
		  p=ptr->iindices[j];
                  /*modified due to rank1 variable*/
                  if ((*prank1>0)&&(*prank1>=blk)){
                      ent = ent*ent;
                  }
		  result.blocks[blk].data.vec[p] += y[i]*ent;
		};
	    }
	  else
	    {
	      for (j=1; j<=ptr->numentries; j++)
		{
		  ent=ptr->entries[j];
                   /*modified due to rank1 variable*/
                  if ((*prank1>0)&&(*prank1>=blk)){
                      for (jj=j; jj <=ptr->numentries ;jj++){
                           p=ijtok(ptr->iindices[j],ptr->jindices[jj],ptr->blocksize);
                           q=ijtok(ptr->jindices[jj],ptr->iindices[j],ptr->blocksize);

                          entex=ptr->entries[jj];

                           result.blocks[blk].data.mat[p] += y[i]*ent*entex;

                            if (p != q){
                               result.blocks[blk].data.mat[q] += y[i]*ent*entex;
                            };
                      }
                  }
                  else{  /* the standard way of computing it*/
                       p=ijtok(ptr->iindices[j],ptr->jindices[j],ptr->blocksize);
                       q=ijtok(ptr->jindices[j],ptr->iindices[j],ptr->blocksize);

                       result.blocks[blk].data.mat[p] += y[i]*ent;
                       if (p != q)
                           result.blocks[blk].data.mat[q] += y[i]*ent;
                  };
                };

            };

	  ptr=ptr->next;
	};
    };

}

