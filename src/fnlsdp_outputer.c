#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "declarations.h"

void initialize_outputer(struct outputinfo *outt,int nu, int ny, int nx, int id){
	//outt=(struct outputinfo *)malloc(sizeof(struct outputinfo));
	outt->nx=nx;
	outt->ny=nx;
	outt->nu=nu;
	//outt->iterate=(double *)calloc(nu*ny+nx*(nx+1),sizeof(double));
	//outt->iterate_deriv=(double *)calloc(nu*ny+nx*(nx+1),sizeof(double));
}

/*void free_outputer(struct outputinfo *outt, int id){
	free(outt->iterate);
	free(outt->iterate_deriv);
	free(outt);
}*/

void print_outputer(int flag, struct outputinfo *outt, int id){
	int nu=outt->nu;
	int ny=outt->ny;
	int nx=outt->nx;	

	int size=nu*ny+nx*(nx+1);

	if(flag==0){
		if(id==0)printf("\tit\tfr\t||dx||\tf\tth1\tth2\ttV\tro\ttipo_iter\n");
	}
	if(flag==1){
		if(id==0)printf("RESTORATION STEP\n");
	}
	if(flag==2){
		if(id==0)printf("FAILS OF STEP 1 (RESTORATION)\n");
	}
	if(flag==3){

		int num_iter=outt->alg_data_numerical_num_iter;
		int fr=outt->alg_data_numerical_fr;
		double ro=outt->alg_data_numerical_rho;
		double f=outt->filter_cand_f;
		double th1=outt->filter_cand_constraints_restr1;    
		double th2=outt->filter_cand_constraints_restr2; 
		double tV=outt->filter_cand_constraints_restr3;    
		char tipo_iter=outt->alg_data_char;
		if(id==0)printf("\t%2d\t%3.0d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%5.2f\t%c\n",num_iter,fr,outt->norm_iterate_deriv,f,th1, th2, tV,ro,tipo_iter);
	}
	if(flag==4){
		double norma_sol_qp=outt->norm_iterate_deriv;
		if(id==0)printf("\n Fin del algoritmo!\n Norma de la solucion de QP(x,ro) %11.8f \n\n",norma_sol_qp);
	}
	if(flag==5){
		if(id==0)printf("\n ro<tlerancia_ro\n");
	}
	if(flag>5)
		if(id==0)printf("ERROR: outputer");
}
