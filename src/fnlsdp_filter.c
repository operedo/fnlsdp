#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "declarations.h"


void initialize_filter(filter *Fil, double ut, double uf){
	Fil->theta=(double *)malloc(sizeof(double)*1);
	Fil->f=(double *)malloc(sizeof(double)*1);
	Fil->theta[0]=ut;
	Fil->f[0]=uf;
	Fil->size=1;
}

void free_filter(filter *Fil){
	free(Fil->theta);
	free(Fil->f);
}

void print_filter(filter *Fil){
	int j=0;
	printf("Filter: size=%d\n",Fil->size);
	for(j=0;j<Fil->size;j++){
		printf("\t%f\t%f\n",Fil->theta[j],Fil->f[j]);
	}
}

int acceptable(filter *Fil,double t, double f, double beta, double gamma)
{
	int nn, verdad=1, j=0;
	nn=Fil->size;
	while(j<nn && verdad==1){
		if( ( t<=beta*Fil->theta[j]  ) || (f+gamma*t <= Fil->f[j]) )
			j++;
		else{
			j=nn+1;
			verdad=0;
		}
	}
	return verdad;	
}
void extract(filter *Fil, double t, double f){
	int nn,cont=0,j=0,i=0;
	nn=Fil->size;
	
	int flag[nn];
	int newSize=0;
	for(j=0;j<nn;j++){
		if( t==Fil->theta[j] && f==Fil->f[j] ){
			flag[j]=1;
			newSize++;
		}
		else
			flag[j]=0;
	}
	newSize=nn-newSize;
	
	double *newTheta=(double *)malloc(sizeof(double)*newSize);
	double *newF=(double *)malloc(sizeof(double)*newSize);

	j=0;
	for(i=0;i<nn;i++){
		if(flag[i]==0){
			newTheta[j]=Fil->theta[i];
			newF[j]=Fil->f[i];
			j++;
		}
	}
	//newTheta[newSize-1]=t;
	//newF[newSize-1]=f;

	free(Fil->theta);
	free(Fil->f);
	Fil->theta=newTheta;
	Fil->f=newF;
	Fil->size=newSize;
}

void add(filter *Fil, double t, double f){
	int nn,cont=0,j=0,i=0;
	nn=Fil->size;
	
	int flag[nn];
	int newSize=0;
	for(j=0;j<nn;j++){
		if( t<=Fil->theta[j] && f<=Fil->f[j] ){
			flag[j]=1;
			newSize++;
		}
		else
			flag[j]=0;
	}
	newSize=nn-newSize+1;
	
	double *newTheta=(double *)malloc(sizeof(double)*newSize);
	double *newF=(double *)malloc(sizeof(double)*newSize);

	j=0;
	for(i=0;i<nn;i++){
		if(flag[i]==0){
			newTheta[j]=Fil->theta[i];
			newF[j]=Fil->f[i];
			j++;
		}
	}
	newTheta[newSize-1]=t;
	newF[newSize-1]=f;

	free(Fil->theta);
	free(Fil->f);
	Fil->theta=newTheta;
	Fil->f=newF;
	Fil->size=newSize;
}

