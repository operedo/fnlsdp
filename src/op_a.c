/*
  Compute A(X).
  */


/* This file was modified compared to version 5.0.
*  Changes conrern the use of parameter rank1 with 
*  nonzero value.
*  
*  Modified by Ivan Ivanov  31 Jan 2007
*/


#include <stdlib.h>
#include "declarations.h"


void op_a(k,constraints,X,result,prank1)
     int k;
     struct constraintmatrix *constraints;
     struct blockmatrix X;
     double *result;
     int *prank1;
{
  int i,j,jj;
  int p,q;
  int blk;
  double ent,entex;
  double *mat;
  double *vec;
  int nume;
  struct sparseblock *ptr;
  double contrib;

  for (i=1; i<=k; i++)
    {
      result[i]=0.0;

      contrib=0.0;
      ptr=constraints[i].blocks;
      while (ptr != NULL)
	{
	  blk=ptr->blocknum;
	  nume=ptr->numentries;

	  if (X.blocks[blk].blockcategory == DIAG)
	    {
	      vec=X.blocks[blk].data.vec;
	      for (j=1; j<=nume; j++)
		{
		  ent=ptr->entries[j];
		  p=ptr->iindices[j];
                  /*modified due to rank1 variable*/
                  if ((*prank1>0) && (*prank1 >= blk)){
                      ent = ent*ent;
                  }
		  contrib += ent*vec[p];
		};
	    }
	  else
	    {
	      mat=X.blocks[blk].data.mat;
	      for (j=1; j<=nume; j++)
		{
		  ent=ptr->entries[j];
                   /*modified due to rank1 variable*/
                  if ((*prank1>0) && (*prank1 >= blk)){
                      for (jj=j; jj <=nume ;jj++){
                           p=ijtok(ptr->iindices[j],ptr->jindices[jj],ptr->blocksize);
                           q=ijtok(ptr->jindices[jj],ptr->iindices[j],ptr->blocksize);
                           if (p == q)
                           {
                               contrib += ent*ent*mat[p];
                            }
                            else
                            {
                               entex=ptr->entries[jj];
                               contrib += ent*entex*(mat[p]+mat[q]);
                             };

                      }
                  }
                  else{  /* the standard way of computing it*/
                        p=ijtok(ptr->iindices[j],ptr->jindices[j],ptr->blocksize);
                        q=ijtok(ptr->jindices[j],ptr->iindices[j],ptr->blocksize);

                       if (p == q)
                       {
                          contrib += ent*mat[p];
                        }
                       else
                       {
                         contrib += ent*(mat[p]+mat[q]);
                        };
                   };
                };

            };


	  ptr=ptr->next;
	};
     
      result[i] += contrib;
    };
}



