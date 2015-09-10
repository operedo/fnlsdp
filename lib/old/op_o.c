/*
 * Compute a matrix representation of the operator 
 *
 * O(.)=A(inv(Z)*A'(.)*X) 
 *
 * The ith colum of the result is O(e_i), where e_i is the vector with a 1 in
 * position i and 0's elsewhere.  
 *
 */


/* This file was modified compared to version 5.0.
*  Changes conrern several issues - the use of parameter rank1 with 
*  nonzero value, distributed way of computing O, etc.
*  
*  Modified by Ivan Ivanov  31 Jan 2007
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "declarations.h"


void
op_o(k, constraints, byblocks, Zi, X, O, O_d, descOd, work1, work2, work3, scapack, prank1)
  int             k;
  struct constraintmatrix *constraints;
  struct sparseblock **byblocks;
  struct blockmatrix Zi;
  struct blockmatrix X;
  double         *O;
  double         *O_d;
  int           *descOd;
  struct blockmatrix work1;
  struct blockmatrix work2;
  struct blockmatrix work3;
  struct scalapackpar scapack;
  int *prank1;
{
  int             i, j;
  int             ii, jj;
  int             ldam;
#ifndef BIT64
  int p,q,r,s;
#else
  long int p,q,r,s;
#endif

  struct sparseblock *ptri;
  struct sparseblock *ptrj;

  int             blocknum;
  int             blocksize;
  double          contrib, contrib_X, contrib_Zi;
  double         *Ziblk;
  double         *Xblk;
  double         *workblk;
  double         *work2blk;
  double         *work3blk;
  double          enti, entj, scale1, scale2;
  int             gi,gj;
 
  int ictxt,nprow, npcol, myrow, mycol,indx_gc;
  int lr, lc, nb, np, ione=1, izero=0;
  
  /* init grid for 2D block cyclic row distribution of O  */
 
  myrow = scapack.myrow;
  mycol = scapack.mycol;
  nprow = scapack.nprow;
  npcol = scapack.npcol;
  nb = scapack.nb;
  ictxt = scapack.ic;
  np = scapack.np;
  
    /*
   * Work out the leading dimension for the array.  Note that we may not want
   * to use k itself, for cache issues. 
   */
  if ((k % 2) == 0)
    ldam = k + 1;
  else
    ldam = k;

  /* determine local dimensions and create the distributed matrix O_d */

  lr = scapack.lr;
  lc = scapack.lc; 

  indx_gc = indxl2g_(&lc,&nb,&mycol,&izero,&npcol);
  if (indx_gc > k){ 
    i=lc-1;
    indx_gc = indxl2g_(&i,&nb,&mycol,&izero,&npcol);
  }
 /*
   * First, zero out the O_d  matrix.
   */
   
  
  for (j = 1; j <= lc; j++){
    for (i = 1; i <= lr; i++)      /* i < j*/
       O_d[ijtok(i, j, lr)] = 0.0;
   
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Cblacs_barrier( ictxt,"All");

  
  /* Loop over i, then the blocks, then j */
  for (i = 1; i <= lr; i++) {

  /* i is transformed to global index gi */
    
    gi=indxl2g_(&i,&nb,&myrow,&izero,&nprow);

    if (gi > k){  
        continue;
    }
    ptri = constraints[gi].blocks; 

    while (ptri != NULL) {
      blocknum = ptri->blocknum;
      blocksize = ptri->blocksize;
    
       /* Diagonal blocks */
      if (ptri->issparse == 1 &&
	  X.blocks[blocknum].blockcategory == DIAG ) {
	Ziblk = Zi.blocks[blocknum].data.vec;
	Xblk = X.blocks[blocknum].data.vec;
        
        ptrj = ptri;
	while (ptrj!=NULL) {
            gj = ptrj->constraintnum;

            if (mycol==indxg2p_(&gj,&nb,&mycol,&izero,&npcol)){

                j = indxg2l_(&gj,&nb,&mycol,&izero,&npcol);
  
	  /*
	   * Do the contribution from constraints i and j of this block. 
	   */
	  contrib = 0.0;
	  p = 1;
	  q = 1;
	  while ((p <= ptri->numentries) &&
		 (q <= ptrj->numentries)) {
	    if (ptri->iindices[p] < ptrj->iindices[q]) {
	      p = p + 1;
	    } else {
	      if (ptri->iindices[p] > ptrj->iindices[q]) {
		q = q + 1;
	      } else {
		/*
		 * Match! 
		 */
                 /*modified due to rank1 variable*/
                  enti = ptri->entries[p];
                  entj = ptrj->entries[q];
                  if ((*prank1>0)&&(*prank1>=blocknum)){
                     enti = enti*enti;
                     entj = entj*entj;
                  }

		contrib += ptri->entries[p] * ptrj->entries[q] *
		  Ziblk[ptri->iindices[p]] *
		  Xblk[ptri->iindices[p]];
		p = p + 1;
		q = q + 1;
	      };
	    };
	  };
	  
          O_d[ijtok(i, j, lr)] += contrib; 

          } /* if p_of_i*/
          else if (gj > indx_gc){
              break;
          }

	  /* Next j */
          ptrj = ptrj->nextbyblock;
        }
      }
      /* Sparse matrix block */
      else if (ptri->issparse == 1 &&
	       X.blocks[blocknum].blockcategory == MATRIX) {
	Ziblk = Zi.blocks[blocknum].data.mat;
	Xblk = X.blocks[blocknum].data.mat;

        ptrj = ptri;
        while (ptrj != NULL) {
            gj = ptrj->constraintnum; 

/*
             * The following prefetch seems to give about a 5% performance
             * improvement on certain problems. e.g. truss8 on a 1200 Mhz
             * Athlon.
             *
             * It won't be compiled unless we're using gcc.
             */

#ifdef __GNUC__
#if (((__GNUC__ == 3) && (__GNUC_MINOR__ > 1)) || (__GNUC__ > 3))
            __builtin_prefetch(ptrj->nextbyblock, 0, 3);
#endif
#endif
          
             if ( mycol==indxg2p_(&gj,&nb,&mycol,&izero,&npcol)){            

                j = indxg2l_(&gj,&nb,&mycol,&izero,&npcol);
          
         /* Only process sparse-sparse pairs. */
       	
	    /*
	     * Do the contribution from constraints i and j of this block. 
	     */
	    contrib = 0.0; contrib_X = 0.0; contrib_Zi = 0.0;
	    for (ii = 1; ii <= ptri->numentries; ii++) {

	      enti = ptri->entries[ii];
	      p = ptri->iindices[ii];
	      q = ptri->jindices[ii];

	      /*
	       * We'll keep the p==q test outside of the inner loop. 
	       *
	       * In what follows, we've made use of the symmetry of Ziblk and
	       * Xblk, by permuting all array indices so that p and q are
	       * last. This means that we stay in columns p and q of Ziblk
	       * and Xblk as much as possible, improving locality. 
	       *
	       */

	      if (p == q) {
		for (jj = 1; jj <= ptrj->numentries; jj++) {
		  entj = ptrj->entries[jj];
		  r = ptrj->iindices[jj];
		  s = ptrj->jindices[jj];
                  /*modified due to rank1 variable*/

                  if ((*prank1>0)&&(*prank1>=blocknum)){

                     if (p == r) {
                    /* here p==r, q==s, when rank1==1 always p==q && r==s */
                         contrib_X += enti * entj * Xblk[ijtok(p, r, blocksize)];
                         contrib_Zi += enti * entj * Ziblk[ijtok(p, r, blocksize)];
                     } else {
                     /* here p!=r, q!=s */

                       contrib_X += enti * entj * Xblk[ijtok(r, p, blocksize)];
                        contrib_Zi += enti * entj * Ziblk[ijtok(r, p, blocksize)];
                     }

                  }
                  else{

  		     if (r == s) {
		       /* here p==q, r==s */
		        contrib += enti * entj *
		        Ziblk[ijtok(r, q, blocksize)] *
		        Xblk[ijtok(s, p, blocksize)];
		      } else {
		       /* here p=q, r!=s */
		       contrib += enti * entj *
		          (Ziblk[ijtok(r, q, blocksize)] *
		           Xblk[ijtok(s, p, blocksize)]
		           +
		           Ziblk[ijtok(s, q, blocksize)] *
		           Xblk[ijtok(r, p, blocksize)]);
                   } 
		  }
		}
	      } else {		/* p!= q */
		for (jj = 1; jj <= ptrj->numentries; jj++) {
		  entj = ptrj->entries[jj];
		  r = ptrj->iindices[jj];
		  s = ptrj->jindices[jj];

		  if (r == s) {
		    /* p!=q, but r=s */
		    contrib += enti * entj *
		      (Ziblk[ijtok(r, q, blocksize)] *
		       Xblk[ijtok(s, p, blocksize)]
		       +
		       Ziblk[ijtok(r, p, blocksize)] *
		       Xblk[ijtok(s, q, blocksize)]);
		  } else {
		    /* here, p!=q and r!=s */
		    contrib += enti * entj *
		      (Ziblk[ijtok(r, q, blocksize)] *
		       Xblk[ijtok(s, p, blocksize)] +
		       Ziblk[ijtok(r, p, blocksize)] *
		       Xblk[ijtok(s, q, blocksize)] +
		       Ziblk[ijtok(s, q, blocksize)] *
		       Xblk[ijtok(r, p, blocksize)] +
		       Ziblk[ijtok(s, p, blocksize)] *
		       Xblk[ijtok(r, q, blocksize)]);
		  }
		}
	      }
	    }
          /*modified due to rank1 variable*/
           if ((*prank1>0)&&(*prank1>=blocknum))
               O_d[ijtok(i, j, lr)] += contrib_X*contrib_Zi;
            else
               O_d[ijtok(i, j, lr)] += contrib;
 
             
          } /* end is mycol*/
          else if (gj > indx_gc){
             break;
          } 
           ptrj = ptrj->nextbyblock;   
	} /* end for*/

      }  /* end else if*/
      /* Dense matrix block */
      else {
	/*
              * put this block into a work matrix.
              */

              blocksize = work1.blocks[blocknum].blocksize;
              workblk = work1.blocks[blocknum].data.mat;
              work2blk = work2.blocks[blocknum].data.mat;
              work3blk = work3.blocks[blocknum].data.mat;
              Xblk = X.blocks[blocknum].data.mat;
              Ziblk = Zi.blocks[blocknum].data.mat;
              for (ii = 0; ii <= blocksize * blocksize - 1; ii++) {
                  workblk[ii] = 0.0; work2blk[ii] = 0.0; work3blk[ii] = 0.0;
              }
              for (ii = 1; ii <= ptri->numentries; ii++) {
                  enti = ptri->entries[ii];
                  p = ptri->iindices[ii];
                  q = ptri->jindices[ii];
                  /*modified due to rank1 */
                  if ((*prank1>0)&&(*prank1>=blocknum)){
                       workblk[ijtok(p,1,blocksize)] = enti;
                   }
                   else{
                        workblk[ijtok(p, q, blocksize)] = enti;
                        if (p != q)
                            workblk[ijtok(q, p, blocksize)] = enti;
                   }
              }
         /*
               * Now, multiply out Zi*work*X.
               */
               scale1 = 1.0;
               scale2 = 0.0;

           /*modified due to rank1 */
              if ((*prank1>0)&&(*prank1>=blocknum)){
                   vec_mult_mat_raw(blocksize,scale1,scale2,Xblk,workblk,work2blk);
                   vec_mult_mat_raw(blocksize,scale1,scale2,Ziblk,workblk,work3blk);
              }
              else{
                   mat_mult_raw(blocksize, scale1, scale2, Ziblk, workblk, work2blk);
                   mat_mult_raw(blocksize, scale1, scale2, work2blk, Xblk, workblk);
              }


       ptrj = ptri; 
       while (ptrj != NULL) {
	  gj = ptrj->constraintnum;

/*
             * The following prefetch seems to give about a 5% performance
             * improvement on certain problems. e.g. truss8 on a 1200 Mhz
             * Athlon.
             *
             * It won't be compiled unless we're using gcc.
             */


#ifdef __GNUC__
#if (((__GNUC__ == 3) && (__GNUC_MINOR__ > 1)) || (__GNUC__ > 3))
            __builtin_prefetch(ptrj->nextbyblock, 0, 3);
#endif
#endif

 
          if ( mycol==indxg2p_(&gj,&nb,&mycol,&izero,&npcol)){

                j = indxg2l_(&gj,&nb,&mycol,&izero,&npcol);
  
  	  /*
	   * Compute if block j is sparse; block j is dense, but block i has
	   * more elements; or, if both blocks have the same number of
	   * elements, and the constraint number of block i is less than or
	   * equal to j + one case extra not in original
	   */
	  if ((ptrj->issparse == 1 ||
	      ptri->numentries > ptrj->numentries ||
	      ptri->numentries == ptrj->numentries ||
		ptri->numentries < ptrj->numentries)&& gi<=gj
	    ) {
            
	    contrib = 0.0; contrib_X = 0.0; contrib_Zi = 0.0;

	    /*modified due to rank1, not checked */
             if ((*prank1>0)&&(*prank1>=blocknum))
               for (ii = 0; ii <= blocksize * blocksize - 1; ii++) {
                     workblk[ii] = 0.0;
                }

	    for (ii = 1; ii <= ptrj->numentries; ii++) {
	      entj = ptrj->entries[ii];
              /*modified due to rank1, not checked */
              if ((*prank1>0)&&(*prank1>=blocknum)){
                  p=ptrj->iindices[ii];
                  workblk[ijtok(p,1,blocksize)] = entj;
              }
              else{
                  p = ijtok(ptrj->iindices[ii], ptrj->jindices[ii], blocksize);
	          q = ijtok(ptrj->jindices[ii], ptrj->iindices[ii], blocksize);
	          contrib += entj * workblk[p];
	          if (p != q) {
		      contrib += entj * workblk[q];
	          };
	        };
              };

              /* modified due to rank1 */
            if ((*prank1>0)&&(*prank1>=blocknum)){
                contrib_X = vec_mult_vec(blocksize,work2blk,workblk);
                contrib_Zi = vec_mult_vec(blocksize,work3blk,workblk);
                contrib += contrib_X*contrib_Zi;
            }

          
              O_d[ijtok(i, j, lr)] += contrib;

           }  
                       	   
          } 	   
          else if (gj > indx_gc){
              break;
          }

           ptrj = ptrj->nextbyblock;
         }  
      }  
      ptri = ptri->next; 

     } 
  }  

    
  MPI_Barrier(MPI_COMM_WORLD);
 

/* make a copy of O_d to O */

   pdlacpy_("All", &ldam, &ldam, O_d, &ione, &ione, descOd, O, &ione, &ione, descOd);
  
}
