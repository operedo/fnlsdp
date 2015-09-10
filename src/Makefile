#
# Use this line to specify your C compiler.
#
CC=mpicc
#CC=gcc4
#
# Specify flags to use in compiling the code.  

#CFLAGS=-O3 -static -fopenmp -DNOSHORTS -DBIT64 -DUSEOPENMP -I../include -I/opt/intel/mkl/10.0.3.020/include
CFLAGS= -g -O0 -I../include -DBIT64
#
#  This is used to build the library.  
#
libfnlsdp.a: fnlsdp_read_data.o fnlsdp_build_mats.o fnlsdp_print_mats.o fnlsdp_solve_qp.o fnlsdp_filter.o fnlsdp_objective_theta.o fnlsdp_restoration.o fnlsdp_algorithm.o fnlsdp_outputer.o readprob.o sdp.o op_o.o psd_feas.o op_a.o op_at.o  Fnorm.o calc_pobj.o calc_dobj.o trace_prod.o zero_mat.o mat_mult.o sym_mat.o copy_mat.o addscaledmat.o  make_i.o allocmat.o initsoln.o initparams.o add_mat.o writesol.o readsol.o easysdp.o writeprob.o solvesys.o makefill.o mat_multsp.o norms.o linesearch.o matvec.o chol.o  qreig.o tweakgap.o freeprob.o packed.o  sortentries.o mat_redist.o
	ar cr libfnlsdp.a  fnlsdp_read_data.o fnlsdp_build_mats.o fnlsdp_print_mats.o fnlsdp_solve_qp.o fnlsdp_filter.o fnlsdp_objective_theta.o fnlsdp_restoration.o fnlsdp_algorithm.o fnlsdp_outputer.o readprob.o sdp.o op_o.o psd_feas.o op_a.o op_at.o  Fnorm.o calc_pobj.o calc_dobj.o trace_prod.o zero_mat.o mat_mult.o sym_mat.o copy_mat.o addscaledmat.o make_i.o allocmat.o initsoln.o initparams.o add_mat.o writesol.o readsol.o easysdp.o writeprob.o solvesys.o makefill.o mat_multsp.o norms.o linesearch.o matvec.o chol.o qreig.o tweakgap.o freeprob.o packed.o sortentries.o mat_redist.o 
	cp libfnlsdp.a ../lib
##
# On some systems, you might need to add after "ar cr libsdp.a ..."
#        ranlib libsdp.a
#
#
# To clean things up.
#
clean:
	rm -f *.o
	rm -f libfnlsdp.a
	rm -f ../lib/libfnlsdp.a
