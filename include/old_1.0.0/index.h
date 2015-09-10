/*
  Declarations needed to handle indexing into Fortran arrays and packed
  arrays.
*/
/* This file was modified compared to version 5.0
*  There are a number of new declarations added.
*  
*  Modified by Ivan Ivanov  31 Jan 2007
*/

#ifndef BIT64

/*
  First, to convert fortran i,j indices into a C vector index.
*/

#define ijtok(iiii,jjjj,lda) ((jjjj-1)*lda+iiii-1)

/* 
   Packed indexing.
 */

#define ijtokp(iii,jjj,lda) ((iii+jjj*(jjj-1)/2)-1)

/*
  Next, to convert C vector index into Fortran i,j indices.
*/

#define ktoi(k,lda) ((k % lda)+1)
#define ktoj(k,lda) ((k/lda)+1)

#else

/*
  First, to convert fortran i,j indices into a C vector index.
*/

#define ijtok(iiii,jjjj,lda) ((jjjj-1L)*lda+iiii-1L)

/* 
   Packed indexing.
 */

#define ijtokp(iii,jjj,lda) (((long int)iii+(long int)jjj*(jjj-1L)/2-1L))

/*
  Next, to convert C vector index into Fortran i,j indices.
*/

#define ktoi(k,lda) (((long int)k % lda)+1L)
#define ktoj(k,lda) (((long int)k/lda)+1L)


#endif

/*  Scalapack global-to-local index transformation */

#define  p_of_i(i,bs,p)  (((int)((i-1)/bs))%p)

#define  l_of_i(i,bs,p)  (int)((i-1)/(p*bs))

#define  x_of_i(i,bs,p)  ((i-1)%bs+1)

