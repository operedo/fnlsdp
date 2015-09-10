# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>

#include "declarations.h"

void test_nelmin (genmatrix *F0, genmatrix *Q0, genmatrix *V0,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np)
{
  int i,j;
  int icount;
  int ifault;
  int kcount;
  int konvge;
  int n;
  int numres;
  double reqmin;
  double *start;
  double *step;
  double *xmin;
  double ynewlo;

  n = F0->nrows*F0->ncols + Q0->nrows*(Q0->nrows+1)/2 + V0->nrows*(V0->nrows+1)/2;

  start = (double *)malloc(sizeof(double)*n);
  step = (double *)malloc(sizeof(double)*n);
  xmin = (double *)malloc(sizeof(double)*n);

  printf("\n");
  printf("TEST01\n");
  printf("  Apply NELMIN to ROSENBROCK function.\n");

	for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			start[i+nu*j]=F0->mat[i+j*nu];
	}

	int count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			start[count+nu*ny]=Q0->mat[i+j*nx];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			start[count+nu*ny+nx*(nx+1)/2]=V0->mat[i+j*nx];
			count++;
		}
	}
	

  reqmin = 1.0E-08;

	for(i=0;i<nu*ny+nx*(nx+1);i++)
  		step[i] = 1.0;

  konvge = 10;
  kcount = 500;

  printf("\n");
  printf("  Starting point X:\n");
  printf("\n");
  for ( i = 0; i < n; i++ )
  {
    printf("  %f\n", start[i]);
  }

printf("entro?\n");
  ynewlo = eval_theta_vec ( start,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );

  printf("\n");
  printf("  F(X) = %f\n", ynewlo);

  nelmin ( eval_theta_vec, n, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );

  printf("\n");
  printf("  Return code IFAULT = %f\n", ifault);
  printf("\n");
  printf("  Estimate of minimizing value X*:\n");
  printf("\n");
  for ( i = 0; i < n; i++ )
  {
    printf("  %f\n", xmin[i]);
  }

  printf("\n");
  printf("  F(X*) = %f\n", ynewlo);

  printf("\n");
  printf("  Number of iterations = %d\n", icount);
  printf("  Number of restarts =   %d\n", numres);

  free(start);
  free(step);
  free(xmin);

  return;
}

//****************************************************************************80

void nelmin ( double fn ( double *x,struct outputinfo *outt ,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np), int n, double *start, double *xmin, double *ynewlo, double reqmin, double *step, int konvge, int kcount, int *icount, int *numres, int *ifault ,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np)

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    FORTRAN77 version by R ONeill
//    C++ version by John Burkardt
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  p = (double *)malloc(sizeof(double)*n*(n+1));
  pstar = (double *)malloc(sizeof(double)*n);
  p2star = (double *)malloc(sizeof(double)*n);
  pbar =(double *)malloc(sizeof(double)*n); 
  y =(double *)malloc(sizeof(double)*(n+1));

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }
    y[n] = fn ( start,outt ,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np);
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = fn ( start,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );
      *icount = *icount + 1;
      start[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      ystar = fn ( pstar,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
        y2star = fn ( p2star,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
          y2star = fn ( p2star,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
              y[j] = fn ( xmin,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
          y2star = fn ( p2star,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
      z = fn ( xmin,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
      z = fn ( xmin,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  free(p);
  free(pstar);
  free(p2star);
  free(pbar);
  free(y);

  return;
}
//****************************************************************************80


void fminsearch (double *step_x,genmatrix *F0, genmatrix *Q0, genmatrix *V0,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np)
{
  	int i,j;
  	int icount;
  	int ifault;
  	int kcount;
  	int konvge;
  	int n;
  	int numres;
  	double reqmin;
  	double *start;
  	double *step;
  	double *xmin;
  	double ynewlo;

  	n = F0->nrows*F0->ncols + Q0->nrows*(Q0->nrows+1)/2 + V0->nrows*(V0->nrows+1)/2;

  	start = (double *)malloc(sizeof(double)*n);
  	step = (double *)malloc(sizeof(double)*n);
  	xmin = (double *)malloc(sizeof(double)*n);

	for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			start[i+nu*j]=F0->mat[i+j*nu];
	}

	int count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			start[count+nu*ny]=Q0->mat[i+j*nx];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			start[count+nu*ny+nx*(nx+1)/2]=V0->mat[i+j*nx];
			count++;
		}
	}
	

  	reqmin = 1.0E-08;

	for(i=0;i<nu*ny+nx*(nx+1);i++)
  		step[i] = 1.0;

  	konvge = 10;
  	kcount = 1000;
 
  	ynewlo = eval_theta_vec ( start,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );

  	nelmin ( eval_theta_vec, n, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );


	for(i=0;i<nu*ny+nx*(nx+1);i++)
		step_x[i]=xmin[i];


    	free(start);
  	free(step);
  	free(xmin);

  	return;
}

void fminsearch_vec (double *step_x,double *x_inicial/*genmatrix *F0, genmatrix *Q0, genmatrix *V0*/,struct outputinfo *outt,double rho,genmatrix *A,genmatrix *B1,genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int nx, int nw, int nu, int ny, int nz, struct scalapackpar scapack,struct paramstruc params, int printlevel,int np)
{
  	int i,j;
  	int icount;
  	int ifault;
  	int kcount;
  	int konvge;
  	int n;
  	int numres;
  	double reqmin;
  	double *start;
  	double *step;
  	double *xmin;
  	double ynewlo;

	//genmatrix F0, Q0, V0;
	//vec2mats(x_inicial,&F0,&Q0,&V0,nu,ny,nx);

  	n = nu*ny+nx*(nx+1);//F0.nrows*F0.ncols + Q0.nrows*(Q0.nrows+1)/2 + V0.nrows*(V0.nrows+1)/2;

  	start = (double *)malloc(sizeof(double)*n);
  	step = (double *)malloc(sizeof(double)*n);
  	xmin = (double *)malloc(sizeof(double)*n);

	for(i=0;i<n;i++)
		start[i]=x_inicial[i];

/*printf("\n");
  printf("  Starting point X:\n");
  printf("\n");
  for ( i = 0; i < n; i++ )
  {
    printf("  %f\n", start[i]);
  }*/
	/*for(i=0;i<nu;i++){
		for(j=0;j<ny;j++)
			start[i+nu*j]=F0.mat[i+j*nu];
	}

	int count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			start[count+nu*ny]=Q0.mat[i+j*nx];
			count++;
		}
	}

	count=0;
	for(j=0;j<nx;j++){
		for(i=j;i<nx;i++){
			start[count+nu*ny+nx*(nx+1)/2]=V0.mat[i+j*nx];
			count++;
		}
	}*/
	

  	reqmin = 1.0E-08;

	for(i=0;i<nu*ny+nx*(nx+1);i++)
  		step[i] = 1.0;

  	konvge = 10;
  	kcount = 1000;
 
  	ynewlo = eval_theta_vec ( start,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );

  	nelmin ( eval_theta_vec, n, start, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault,outt,rho,A,B1,B,C1,C,D11,D12,D21, nx, nw, nu, ny, nz, scapack, params, printlevel, np );


	for(i=0;i<nu*ny+nx*(nx+1);i++)
		step_x[i]=xmin[i];

  /*printf("\n");
  printf("  Return code IFAULT = %f\n", ifault);
  printf("\n");
  printf("  Estimate of minimizing value X*:\n");
  printf("\n");
  for ( i = 0; i < n; i++ )
  {
    printf("  %f\n", xmin[i]);
  }*/

  //printf("\n");
  //printf("  F(X*) = %f\n", ynewlo);

  //printf("\n");
  //printf("  Number of iterations = %d\n", icount);
  //printf("  Number of restarts =   %d\n", numres);


    	free(start);
  	free(step);
  	free(xmin);

  	return;
}

