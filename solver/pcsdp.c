/*
 * Modified with additional calls to MPI and ScaLAPACK for distributed parallel 
 * computations.
 *
 * Modified by Ivan Ivanov  31 Jan 2007
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "declarations.h"

int main(argc,argv)
     int argc;
     char *argv[];
{
  int ret;
  int n;
  int k;
  struct blockmatrix C;
  double *a;
  struct constraintmatrix *constraints;
  struct blockmatrix X,Z;
  double *y;
  double pobj,dobj;

  struct paramstruc params;
  int printlevel;

  struct scalapackpar scapack;

  int id, np;
    

 /* Initialize the process grid */

  MPI_Init(&argc,&argv);
  MPI_Comm_rank (MPI_COMM_WORLD,&id);
  MPI_Comm_size (MPI_COMM_WORLD,&np);

  scapack.id = id;
  scapack.np = np;

  switch (scapack.np)
        {
        case 1: scapack.nprow=1; scapack.npcol=1;
          break;
        case 2: scapack.nprow=2; scapack.npcol=1;
          break;
        case 4: scapack.nprow=2; scapack.npcol=2;
          break;
        case 8: scapack.nprow=4; scapack.npcol=2;
          break;
        case 16: scapack.nprow=4; scapack.npcol=4;
          break;
        case 32: scapack.nprow=8; scapack.npcol=4;
          break;
        case 64: scapack.nprow=8; scapack.npcol=8;
          break;
        default:
          if (scapack.id==0)
              printf("Can not setup %d processors to a grid.\nPlease use 1,2,4,8,9,16,32 or 64 nodes to run or modify csdp.c file. \n",scapack.np);
          MPI_Finalize();
          return(1);
        };

   sl_init_(&scapack.ic,&scapack.nprow,&scapack.npcol);
   Cblacs_gridinfo(scapack.ic,&scapack.nprow,&scapack.npcol,&scapack.myrow,&scapack.mycol);


  /*
   * Check that we have enough arguments.
   */

  if ((argc < 2) || (argc > 4)) 
    {
      printf("PCSDP 1.0\n");
      printf(" \n");
      printf("Usage: \n");
      printf("\n");
      printf("prun -v -np 2 pcsdp [-o <output file>] <input problem> [<final solution>] [<initial solution>]\n");
      exit(1);
    };
  
  
 /*
   * First, initialize the parameters.
   */
   initparams(&params,&printlevel);

  /*
   * Second, read the problem in.
   */

   ret=read_prob(argv[1],&n,&k,&C,&a,&constraints,1,&params.rank1,scapack.id);

 

  if (ret != 0)
    {
      printf("Giving up.\n");
      exit(10);
    };

  if (argc == 4)
    {
      ret=read_sol(argv[3],n,k,C,&X,&y,&Z);
      if (ret != 0)
	{
	  printf("Giving up.\n");
	  exit(10);
	};
    }
  else
    {
      initsoln(n,k,C,a,constraints,&X,&y,&Z,&params.rank1);
    };

  /*
   * Call the solver.
   */
  
 
  ret=easy_sdp(n,k,C,a,constraints,0.0,&X,&y,&Z,&pobj,&dobj,scapack,params,printlevel);


  /*
   * Write out the solution if necessary.
   */
  if (argc >= 3)
    write_sol(argv[2],n,k,X,y,Z);

/*
   * Free up memory.
   */

  free_prob(n,k,C,a,constraints,X,y,Z);

  Cblacs_gridexit(scapack.ic);

  MPI_Finalize();

  return(ret);
}

