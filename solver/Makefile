#
# Use this line to specify the C and Fortran compilers.
#
CC=mpicc
#gcc4
F77=mpif77
#
# Use this line to specify options for the C compiler.  You'll probably
# want to turn on optimizations. You may also have to use some of the 
# following flags:
#
#  -DCAPSBLAS         if BLAS routine names are capitalized.
#  -DCAPSLAPACK       if LAPACK routine names are capitalized.
#  -DNOUNDERBLAS      if BLAS routine names have no underscore.
#  -DNOUNDERLAPACK    if LAPACK routine names have no underscore.
#  -DBIT64            For I32LP64 systems.
#
CFLAGS= -g -O0 -I../include  -DBIT64 
#CFLAGS=-O3 -static -fopenmp -DNOSHORTS -DBIT64 -DUSEOPENMP -I../include -I/opt/intel/mkl/10.0.3.020/include
#
# Use this line to specify options for the Fortran compiler. In fact
# Fortran is used only to link the code with the external libraries.
# In case  you expirience problems as second_underscore in the filenames 
# of you BLAS, LAPACK or ScaLAPACk libraties it is a good idea to
# use the option -fno-second-underscore.
# 
FFLAGS=  #-fno-second-underscore
#
# Use this line to specify where the SDP and linear algebra libraries are
# to be found. 
#
# -L../lib            look in the ../lib directory
# -lsdp               get libsdp.a
# -llapack            get liblapack.a
# -lblas              get libblas.a
# -lm                 Get C math library.
#
# It's extremely likely that you'll have to change the LIBS= line for
# your particular system.  
# BLACS library might have slightly different names than the specified here,
# so use the proper names available on your cluster. Please, pay attention 
# on the order of ScaLAPACK and BLACS libs.
#
#MKLLIB=/opt/intel/mkl/10.0.3.020/lib/64
#LIBS=-static -L../lib -lsdp -lfnlsdp
#SCALAPACKLIB= -L$(MKLLIB) 
##BLAS= -lmkl_intel_lp64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack -lmkl_ipf -lguide -lpthread  
#BLAS= -lmkl_lapack -lgomp -lrt -lmkl_ipf -lguide -lpthread -lgfortran 
LIBLM= -lm


MKLLIB=/opt/intel/mkl/10.1.0.015/lib/em64t
LIBS= -L../lib -lfnlsdp -L$(MKLLIB) 
#SCALAPACKLIB= -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack -lmkl -lguide -lpthread
SCALAPACKLIB= -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack -lmkl -lguide -lpthread

#
#  This builds the pcsdp program.  
#
fnlsdp: fnlsdp.o 
	$(CC) $(CFLAGS) fnlsdp.o $(LIBS) $(SCALAPACKLIB) $(LIBLM)  -o fnlsdp
	rm -f *.o
#
# To clean out the directory:
#
clean:
	rm -f *.o*
	rm -f *.e*
	rm -f *.po*
	rm -f *.pe*
	rm -f fnlsdp
