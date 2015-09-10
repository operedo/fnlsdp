#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "declarations.h"

#define min(a,b) ((a)>=(b)?(b):(a))
#define max(a,b) ((a)<(b)?(b):(a))

void algorithm(const char *code, const char *tag, const char *dir_compleib,const char *dir_points , int max_iter, double romax, double beta, double gamma, double sigma, double tolerancia_ro, int step_1_fail,  struct scalapackpar scapack,struct paramstruc params, int printlevel,int id){
	
	struct outputinfo outt;

	/* compleib data matrices */
	genmatrix A, B1, B,C1, C,D11,D12,D21;
	/* compleib initial point */
	genmatrix F0, Q0, V0;
	/* compleib matrix dimensions*/
	int nx,nw,nu,nz,ny;

	char tipo_iter;
	
	char* solutionQP="solutionQP.dat-s";

	int i;
	int stop=0;
	int step=0;
	int step_busca=0;
	int very_first_step=0;
	int num_iter=0;
	int cont=1;
	int paso=1;
	int ac=0;
	int fr=100;
	double rho=10.0;
	double rrho=1.0;
	double fobj;
	double t,f;
	double t1,f1;
	
	filter Fil,Fil_aux;
	double *x_inicial,*x_current,*dx,*suma_x_current_dx;

	load_compleib(code, dir_compleib, &A, &B1, &B, &C1, &C, &D11, &D12, &D21, &nx, &nw, &nu, &nz, &ny, id);
	printf("loading %s %s %s\n",code,tag,dir_points);
	load_initial_point(code,tag,dir_points,&F0,&Q0,&V0,nx,nu,ny, id);

	
	romax=romax*norm2(nx,A.mat)*norm2(nx,A.mat);


	initialize_outputer(&outt,nu,ny,nx,id);
	outt.alg_data_numerical_num_iter=1;
	outt.alg_data_numerical_fr=1;
	outt.alg_data_numerical_rho=1;
	outt.filter_cand_t=500;
	outt.filter_cand_f=500;
	outt.alg_data_char='w';

	print_outputer(0,&outt,id);

	x_inicial=(double *)calloc(nu*ny+nx*(nx+1),sizeof(double));
	x_current=(double *)calloc(nu*ny+nx*(nx+1),sizeof(double));
	dx=(double *)calloc(nu*ny+nx*(nx+1),sizeof(double));
	suma_x_current_dx=(double *)calloc(nu*ny+nx*(nx+1),sizeof(double)); 

	//printf("Begin main iteration\n");
	while(stop!=1 && num_iter<=max_iter){
/**********STEP 0*************/
		//printf("ITER: %d\n",num_iter);
		if(step==0){
			//printf("STEP 0:\n");
			num_iter=1;
			initialize_filter(&Fil,500.0,-500.0);
			initialize_filter(&Fil_aux,500.0,-500.0);
			/*printf("F0:\n");	
			print_genmatrix(&F0);			
			printf("Q0:\n");	
			print_genmatrix(&Q0);			
			printf("V0:\n");	
			print_genmatrix(&V0);*/			

			step=1;
		}
/**********STEP 1*************/
		if(step==1){
			//printf("STEP 1:\n");
			double *step_x=(double *)calloc(nu*ny+nx*(nx+1),sizeof(double));
			double *step_dx=(double *)calloc(nu*ny+nx*(nx+1),sizeof(double));
			ac=0;
			fr=100;
			step_busca=0;
			while(((fr!=0 && fr!=3 && fr!=7)  || ac==0) && step_busca<=step_1_fail){
				if(very_first_step==1){
					very_first_step=0;
					mats2vec(step_x,&F0,&Q0,&V0,nu,ny,nx);
				}
				else{
					if(step_busca==0){
						for(i=0;i<nu*ny+nx*(nx+1);i++)
							x_inicial[i]=x_current[i];	
					}
					else{
						for(i=0;i<nu*ny+nx*(nx+1);i++)
							x_inicial[i]=step_x[i];	
					}
					vec2mats(x_inicial,&F0,&Q0,&V0,nu,ny,nx);
					/*printf("x_inicial->F0:\n");	
					print_genmatrix(&F0);			
					printf("x_inicial->Q0:\n");	
					print_genmatrix(&Q0);			
					printf("x_inicial->V0:\n");	
					print_genmatrix(&V0);*/	
					fminsearch_vec(step_x,x_inicial,&outt,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id);
				}

				f=eval_f_vec(step_x,&outt,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id);
				t=eval_theta_vec(step_x,&outt,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id);
				//printf("t=%f, f=%f\n",t,f);
				ac=acceptable(&Fil,t,f,beta,gamma);
				//printf("ac=%d\n",ac);
				rho=1.0;
				if(ac==1){
					while((fr!=0 && fr!=3 && fr!=7) && rho<romax){
						rho=2*rho;
						fr=solve_qp_vec(code,"01",step_dx, step_x /*&F0,&Q0, &V0*/,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id,solutionQP);	
						//printf("fr=%d\n",fr);
					}
				}
				step_busca++;
			}
		
			if((fr!=0 && fr!=3 && fr!=7) || ac==0){
				print_outputer(2,&outt,id);
				//printf("FAILS OF STEP 1 (RESTORATION)\n");
				stop=1;
			}
			else{
				for(i=0;i<nu*ny+nx*(nx+1);i++){
					x_current[i]=step_x[i];	
					//printf("x_current[%d]=%f\n",i,x_current[i]);
				}
				for(i=0;i<nu*ny+nx*(nx+1);i++){
					dx[i]=step_dx[i];	
					//printf("dx[%d]=%f\n",i,dx[i]);
				}
			}
			if(stop==0){
				cont++;
				paso=1;
				step=3;

				outt.alg_data_numerical_num_iter=num_iter;
				outt.alg_data_numerical_fr=fr;
				outt.alg_data_numerical_rho=rho;
				outt.filter_cand_f=f;
				outt.filter_cand_t=t;
				outt.norm_iterate_deriv=norm2(nu*ny+nx*(nx+1),dx);
				printf("norm_iterate_deriv=%f\n",outt.norm_iterate_deriv);
				print_outputer(3,&outt,id);	
			}
			free(step_x);
			free(step_dx);
		}
/**********STEP 2*************/
		if(step==2){
			//printf("STEP 2:\n");
			if(paso!=1){
				f=eval_f_vec(x_current,&outt,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id);
				t=eval_theta_vec(x_current,&outt,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id);	
				fr=100;
				fr=solve_qp_vec(code,"01",dx, x_current /*&F0,&Q0, &V0*/,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id,solutionQP);	
				cont++;

				outt.alg_data_numerical_num_iter=num_iter;
				outt.alg_data_numerical_fr=fr;
				outt.alg_data_numerical_rho=rho;
				outt.filter_cand_f=f;
				outt.filter_cand_t=t;
				outt.norm_iterate_deriv=norm2(nu*ny+nx*(nx+1),dx);
				print_outputer(3,&outt,id);
			}
			if((fr!=0 && fr!=3 && fr!=7)){
				int v=acceptable(&Fil,t,f,beta,gamma);
				if(v==1){
					add(&Fil,t,f);
					add(&Fil_aux,t,f);
					outt.alg_data_char='h';
				}
				num_iter++;
				step=1;
			}
			else
				step=3;
		}
/**********STEP 3*************/
		if(step==3){
			//printf("STEP 3:\n");
			double norma_sol_qp=norm2(nu*ny+nx*(nx+1),dx);
			if(norma_sol_qp<=tolerancia_ro && (fr==0 || fr==3 || fr==7)){
				stop=1;
				outt.norm_iterate_deriv=norm2(nu*ny+nx*(nx+1),dx);
				print_outputer(4,&outt,id);
			}
			else
				step=4;
		}
/**********STEP 4*************/
		if(step==4){
			//printf("STEP 4:\n");
			for(i=0;i<nu*ny+nx*(nx+1);i++){
				suma_x_current_dx[i]=x_current[i]+dx[i];
			}	
			f1=eval_f_vec(suma_x_current_dx,&outt,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id);
			t1=eval_theta_vec(suma_x_current_dx,&outt,rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id);	
			add(&Fil_aux,t,f);
			int ac=acceptable(&Fil_aux,t1,f1,beta,gamma);
			extract(&Fil_aux,t,f);
			if(ac==0){
				rho=rho/2;
				if(rho<tolerancia_ro & (fr==0 || fr==3 || fr==7)){
					stop=1;
					outt.alg_data_char='e';
					outt.alg_data_numerical_num_iter=num_iter;
					outt.alg_data_numerical_fr=fr;
					outt.alg_data_numerical_rho=rho;
					outt.filter_cand_f=f;
					outt.filter_cand_t=t;
					outt.norm_iterate_deriv=norm2(nu*ny+nx*(nx+1),dx);
					print_outputer(3,&outt,id);
					print_outputer(5,&outt,id);
				}
				paso=4;
				tipo_iter='i';
				step=2;
			}
			else
				step=5;
		}
/**********STEP 5*************/
		if(step==5){
			//printf("STEP 5:\n");
			stop=1;
			fobj=eval_nabla_f_vec(dx,x_current, rho,&A,&B1,&B,&C1,&C,&D11,&D12,&D21,nx,nw,nu,ny,nz,scapack,params,printlevel,id);
			double b14= f + sigma*fobj;
			if(fobj<0 && b14<f1){
				rho=rho/2;
				step=2;
				tipo_iter='i';
				paso=5;
			}	
			else{
				step=6;
			}
		}
/**********STEP 6*************/
		if(step==6){
			//printf("STEP 6:\n");
			tipo_iter='f';
			if(fobj>=0){
				tipo_iter='h';
				int v=acceptable(&Fil,t,f,beta,gamma);
				if(v==1){
					add(&Fil,t,f);
					add(&Fil_aux,t,f);
				}
			}
			step=7;
		}
/**********STEP 7*************/
		if(step==7){
			//printf("STEP 7:\n");
			for(i=0;i<nu*ny+nx*(nx+1);i++){
				x_current[i]=x_current[i]+dx[i];
			}	
			num_iter++;
			double maxtmp=max(2*rho,rrho);
			rho=min(romax,maxtmp);
			step=2;
			paso=7;
		}
	}

	vec2mats(x_current,&F0,&Q0,&V0,nu,ny,nx);

	printf("F:\n");
	print_genmatrix(&F0);
	printf("Q:\n");
	print_genmatrix(&Q0);
	printf("V:\n");
	print_genmatrix(&V0);

	free(suma_x_current_dx);
	free(dx);
	free(x_inicial);
	free(x_current);
	//free_outputer(&outt,id);
	free_filter(&Fil);
	free_initial_point(&F0,&Q0,&V0, id);
	free_compleib(&A, &B1, &B, &C1, &C, &D11, &D12, &D21, id);
}
