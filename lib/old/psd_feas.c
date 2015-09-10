/*
  Check a primal dual solution for strict feasibility.  Return 1 if
  feasible, and 0 if not.  
*/

/* This file was modified compared to version 5.0.
*  Changes conrern declaration of the functions.
*  
*  Modified by Ivan Ivanov  31 Jan 2007
*/

#include <stdio.h>
#include <math.h>
#include "declarations.h"

/*
 * This routine computes the relative primal infeasibility.
 */
double pinfeas(k,constraints,X,a,workvec,prank1)
     int k;
     struct constraintmatrix *constraints;
     struct blockmatrix X;
     double *a;
     double *workvec;
     int *prank1;
{
  double nrme;
  double nrma;
  int i;

  /*
   *  First, check that A(X)=a.
   */
 
  op_a(k,constraints,X,workvec,prank1);
  nrma=norm2(k,a+1);
  for (i=1; i<=k; i++)
    workvec[i]=workvec[i]-a[i];
  nrme=norm2(k,workvec+1);

  return(nrme/(1.0+nrma));

}


double dinfeas(k,C,constraints,y,Z,work1,prank1)
     int k;
     struct blockmatrix C;
     struct constraintmatrix *constraints;
     double *y;
     struct blockmatrix Z;
     struct blockmatrix work1;
     int *prank1;
{
  double nrme;
  double nrmC;

  /*
   * Next, check that A'(y)-C=Z
   */

  zero_mat(work1);
  
  op_at(k,y,constraints,work1,prank1);

  addscaledmat(work1,-1.0,C,work1);
  addscaledmat(work1,-1.0,Z,work1);

  /*
    Now, we've got the error in workn1.  We'll compute the F norm of this
    error and compare it to the F norm of C.
    */

  nrme=Fnorm(work1);

  nrmC=Fnorm(C);

  return(nrme/(1+nrmC));

}

double dimacserr3(k,C,constraints,y,Z,work1,prank1)
     int k;
     struct blockmatrix C;
     struct constraintmatrix *constraints;
     double *y;
     struct blockmatrix Z;
     struct blockmatrix work1;
     int *prank1;
{
  double nrme;

  /*
   * Next, check that A'(y)-C=Z
   */

  zero_mat(work1);

  op_at(k,y,constraints,work1,prank1);

  addscaledmat(work1,-1.0,C,work1);
  addscaledmat(work1,-1.0,Z,work1);

  /*
    Now, we've got the error in workn1.  We'll compute the F norm of this
    error and compare it to the F norm of C.
    */

  nrme=Knorm(work1);

  return(nrme/(1+matinfnorm(C)));

}


