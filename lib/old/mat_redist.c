/* 
 * Redistributes a matrix from the root node to all the others.  
 * Makes a copy of a distributed matrix to the root node.
 */

/* This file is completely new and not part of version 5.0.
*  
*  Modified by Ivan Ivanov  31 Jan 2007
*/ 


#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "declarations.h"

#define TINY 1.0e-10

int max (int a, int b){
if (a>b) return(a); else return(b);
}


void 
ddistr(int ictxt, int n, int numc, int nb, double *A , double *A_d, int *descAd ){

int RootNodeic, ione=1, izero=0, isRootNode=0, nru, info;
int nprow, npcol, myrow, mycol,descA[9], itemp;
int i,k;

/*
#ifdef NOUNDERLAPACK
 sl_init__(&RootNodeic,&ione, &ione);
#else
  sl_init__(&RootNodeic,&ione, &ione);
#endif
*/



   
  sl_init_(&RootNodeic,&ione, &ione);
  
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);  

  if (myrow==0 && mycol==0) { isRootNode=1;}

  if (isRootNode){
     nru = numroc_(&n, &n, &myrow,&izero, &nprow );
     itemp = max(1, nru);
     descinit_(descA, &n, &numc, &n, &n, &izero, &izero, &RootNodeic, &itemp, &info);
  } 
  else{
     k=0;
     for(i=0;i<9;i++){ 
       descA[k]=0;
       k++;
     }
     descA[1]=-1;
  }

  nru = numroc_(&n, &nb, &myrow, &izero, &nprow);
  itemp = max(1,nru);
  descinit_(descAd, &n, &numc, &nb, &nb, &izero, &izero, &ictxt,&itemp, &info);  

  pdgemr2d_( &n, &numc, A, &ione, &ione, descA, A_d, &ione, &ione, descAd, &ictxt);
  
 if (isRootNode){ 
     Cblacs_gridexit(RootNodeic);
  }
}


void 
dgather(int ictxt, int n, int numc, int nb, double *A, double *A_d, int *descAd){

int RootNodeic, ione=1, izero=0, isRootNode=0, nru, info;
int nprow, npcol, myrow, mycol, descA[9], itemp;
int i,k;

   sl_init_(&RootNodeic, &ione, &ione);

   Cblacs_gridinfo(ictxt, &nprow,&npcol, &myrow, &mycol);

   if (myrow==0 && mycol ==0){ isRootNode = 1;}

   if(isRootNode){
      nru = numroc_(&n, &n, &myrow, &izero, &nprow);
      itemp = max(1,nru);
      descinit_(descA, &n, &numc, &n, &n, &izero, &izero, &RootNodeic, &itemp, &info );
   }
   else{
      k=0;
      for(i=0;i<9;i++){
          descA[k]=0;
          k++;
      }
      descA[1]=-1;
    }

   pdgemr2d_(&n,&numc,A_d,&ione, &ione, descAd, A, &ione, &ione, descA, &ictxt );

   if (isRootNode){
       Cblacs_gridexit(RootNodeic);
   }
  
}


