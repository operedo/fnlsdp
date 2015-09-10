/*
  Solve a system of equations using the Cholesky factorization of A.  
  Note that we assume that A is positive definite and that A has already
  been factored.
*/

/* This file was modified compared to version 5.0.
*  It is completely changed for a distributed computations only.
*  
*  Modified by Ivan Ivanov  31 Jan 2007
*/


#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"


int solvesys(m,ldam,descOd,O_d,RHS_d,rhs,scapack)
     int m;
     int ldam;
     int *descOd;
     double *O_d; 
     double *RHS_d;
     double *rhs;
     struct scalapackpar scapack;
{
  int incx;
  int info;
 
  int ione=1;
  int descRHSd[9];
  

  incx=1;


  if ((scapack.myrow>-1)&(scapack.mycol>-1)&(scapack.myrow<scapack.nprow)&(scapack.mycol<scapack.npcol)){

      MPI_Barrier(MPI_COMM_WORLD);
      ddistr(scapack.ic, m, ione, scapack.nb, rhs+1, RHS_d, descRHSd);
   
      Cblacs_barrier(scapack.ic,"All");

#ifdef NOUNDERLAPACK
  #ifdef CAPSLAPACK
         PDPOTRS("U",&m,&ione,O_d,&ione,&ione,descOd, RHS_d,&ione,&ione,descRHSd, &info);
  #else
         pdpotrs("U",&m,&ione,O_d,&ione,&ione,descOd, RHS_d,&ione,&ione,descRHSd, &info);
  #endif
#else
  #ifdef CAPSLAPACK
         PDPOTRS_("U",&m,&ione,O_d,&ione,&ione,descOd, RHS_d,&ione,&ione,descRHSd, &info);
  #else
         pdpotrs_("U",&m,&ione,O_d,&ione,&ione,descOd, RHS_d,&ione,&ione,descRHSd, &info);
  #endif
#endif
       MPI_Barrier(MPI_COMM_WORLD);
       
       dgather(scapack.ic, m, ione, scapack.nb, rhs+1, RHS_d, descRHSd);
     }
           
     MPI_Bcast(rhs,m+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
   	  
           if (info != 0)
	     {
	       return(6);
	     };

	   return(0);
}


