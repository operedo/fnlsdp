#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "declarations.h"
//#include "fnlsdp.h"

/*typedef struct _genmatrix {
  int nrows;
  int ncols;
  double *mat;
} genmatrix;
*/

void load_genmatrix(const char* filenameA, genmatrix *A, int nrows, int ncols, int pp){
	FILE *fp;
	char buffer[128];
	int bufferint,linecounter;

	/* read A file*/
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameA);
	if ((fp = fopen(filenameA, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameA);
		exit(1);
	}
	A->mat=(double *)calloc((nrows)*(ncols),sizeof(double));
	A->nrows=nrows;
	A->ncols=ncols;
	linecounter=0;
	const char *dot = ".";
	char *ret;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			A->mat[linecounter]=atof(buffer);
		else
			A->mat[linecounter]=atoi(buffer);
		if(DEBUG_FNLSDP && pp==0)printf("A->mat[%d]=%f\n",linecounter, A->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);
}

void load_compleib(const char *code, const char *dir, genmatrix *A, genmatrix *B1, genmatrix *B, genmatrix *C1, genmatrix *C, genmatrix *D11, genmatrix *D12, genmatrix *D21, int *nx, int *nw, int *nu, int *nz, int *ny, int pp)
{
	FILE *fp;
	char buffer[128];
	int bufferint;

	const char *Dimsstr="_dims.dat";
	const char *Astr="_A.dat";
	const char *B1str="_B1.dat";
	const char *Bstr="_B.dat";
	const char *C1str="_C1.dat";
	const char *Cstr="_C.dat";
	const char *D11str="_D11.dat";
	const char *D12str="_D12.dat";
	const char *D21str="_D21.dat";
	
	char filenameDims[80]="";
	char filenameA[80]="";
	char filenameB1[80]="";
	char filenameB[80]="";
	char filenameC1[80]="";
	char filenameC[80]="";
	char filenameD11[80]="";
	char filenameD12[80]="";
	char filenameD21[80]="";
	
	/* read dims file*/
	strcat(filenameDims,dir);
	strcat(filenameDims,code);
	strcat(filenameDims,Dimsstr);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameDims);
	if ((fp = fopen(filenameDims, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameDims);
		exit(1);
	}

	int linecounter=0;
	while(fscanf(fp,"%d\n", &bufferint)!=EOF){
		if(linecounter==0) *nx=bufferint;
		if(linecounter==1) *nw=bufferint;
		if(linecounter==2) *nu=bufferint;
		if(linecounter==3) *nz=bufferint;
		if(linecounter==4) *ny=bufferint;
		linecounter++;	
	}
	if(DEBUG_FNLSDP && pp==0)printf("%d %d %d %d %d\n",*nx, *nw, *nu, *nz, *ny);
	fclose(fp);

	/* read A file*/
	strcat(filenameA,dir);
	strcat(filenameA,code);
	strcat(filenameA,Astr);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameA);
	if ((fp = fopen(filenameA, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameA);
		exit(1);
	}
	A->mat=(double *)malloc(sizeof(double)*(*nx)*(*nx));
	A->nrows=*nx;
	A->ncols=*nx;
	linecounter=0;
	const char *dot = ".";
	char *ret;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			A->mat[linecounter]=atof(buffer);
		else
			A->mat[linecounter]=atoi(buffer);
		if(DEBUG_FNLSDP && pp==0)printf("A->mat[%d]=%f\n",linecounter, A->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);
	
	/* read B1 file*/
	strcat(filenameB1,dir);
	strcat(filenameB1,code);
	strcat(filenameB1,B1str);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameB1);
	if ((fp = fopen(filenameB1, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameB1);
		exit(1);
	}
	B1->mat=(double *)malloc(sizeof(double)*(*nx)*(*nw));
	B1->nrows=*nx;
	B1->ncols=*nw;
	linecounter=0;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			B1->mat[linecounter]=atof(buffer);
		else
			B1->mat[linecounter]=atoi(buffer);		
		if(DEBUG_FNLSDP && pp==0)printf("B1->mat[%d]=%f\n", linecounter,B1->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);

	/* read B file */
	strcat(filenameB,dir);
	strcat(filenameB,code);
	strcat(filenameB,Bstr);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameB);
	if ((fp = fopen(filenameB, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameB);
		exit(1);
	}
	B->mat=(double *)malloc(sizeof(double)*(*nx)*(*nu));
	B->nrows=*nx;
	B->ncols=*nu;	
	linecounter=0;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			B->mat[linecounter]=atof(buffer);
		else
			B->mat[linecounter]=atoi(buffer);		
		if(DEBUG_FNLSDP && pp==0)printf("B->mat[%d]=%f\n", linecounter,B->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);

	/* read C1 file */
	strcat(filenameC1,dir);
	strcat(filenameC1,code);
	strcat(filenameC1,C1str);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameC1);
	if ((fp = fopen(filenameC1, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameC1);
		exit(1);
	}
	C1->mat=(double *)malloc(sizeof(double)*(*nz)*(*nx));
	C1->nrows=*nz;
	C1->ncols=*nx;	
	linecounter=0;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			C1->mat[linecounter]=atof(buffer);
		else
			C1->mat[linecounter]=atoi(buffer);		
		if(DEBUG_FNLSDP && pp==0)printf("C1->mat[%d]=%f\n", linecounter,C1->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);

	/* read C file*/
	strcat(filenameC,dir);
	strcat(filenameC,code);
	strcat(filenameC,Cstr);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameC);
	if ((fp = fopen(filenameC, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameC);
		exit(1);
	}
	C->mat=(double *)malloc(sizeof(double)*(*ny)*(*nx));
	C->nrows=*ny;
	C->ncols=*nx;	
	linecounter=0;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			C->mat[linecounter]=atof(buffer);
		else
			C->mat[linecounter]=atoi(buffer);		
		if(DEBUG_FNLSDP && pp==0)printf("C->mat[%d]=%f\n", linecounter,C->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);

	/* read D11 file*/
	strcat(filenameD11,dir);
	strcat(filenameD11,code);
	strcat(filenameD11,D11str);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameD11);
	if ((fp = fopen(filenameD11, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameD11);
		exit(1);
	}
	D11->mat=(double *)malloc(sizeof(double)*(*nz)*(*nw));
	D11->nrows=*nz;
	D11->ncols=*nw;	
	linecounter=0;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			D11->mat[linecounter]=atof(buffer);
		else
			D11->mat[linecounter]=atoi(buffer);		
		if(DEBUG_FNLSDP && pp==0)printf("D11->mat[%d]=%f\n", linecounter,D11->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);

	/* read D12 file*/
	strcat(filenameD12,dir);
	strcat(filenameD12,code);
	strcat(filenameD12,D12str);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameD12);
	if ((fp = fopen(filenameD12, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameD12);
		exit(1);
	}
	D12->mat=(double *)malloc(sizeof(double)*(*nz)*(*nu));
	D12->nrows=*nz;
	D12->ncols=*nu;	
	linecounter=0;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			D12->mat[linecounter]=atof(buffer);
		else
			D12->mat[linecounter]=atoi(buffer);		
		if(DEBUG_FNLSDP && pp==0)printf("D12->mat[%d]=%f\n", linecounter,D12->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);

	/* read D21 file*/
	strcat(filenameD21,dir);
	strcat(filenameD21,code);
	strcat(filenameD21,D21str);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameD21);
	if ((fp = fopen(filenameD21, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameD21);
		exit(1);
	}
	D21->mat=(double *)malloc(sizeof(double)*(*ny)*(*nw));
	D21->nrows=*ny;
	D21->ncols=*nw;	
	linecounter=0;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			D21->mat[linecounter]=atof(buffer);
		else
			D21->mat[linecounter]=atoi(buffer);		
		if(DEBUG_FNLSDP && pp==0)printf("D21->mat[%d]=%f\n", linecounter,D21->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);

}

void load_initial_point(const char *code, const char *id, const char *dir, genmatrix *F0, genmatrix *Q0, genmatrix *V0, int nx, int nu, int ny, int pp){
	FILE *fp;
	char buffer[128];

	const char *F0str="_F0.dat";
	const char *Q0str="_Q0.dat";
	const char *V0str="_V0.dat";
	
	char filenameF0[80]="";
	char filenameQ0[80]="";
	char filenameV0[80]="";
	
	/* read F0 file*/
	strcat(filenameF0,dir);
	strcat(filenameF0,code);
	strcat(filenameF0,"_");
	strcat(filenameF0,id);
	strcat(filenameF0,F0str);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameF0);
	if ((fp = fopen(filenameF0, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameF0);
		exit(1);
	}
	F0->mat=(double *)malloc(sizeof(double)*(nu)*(ny));
	F0->nrows=nu;
	F0->ncols=ny;
	int linecounter=0;
	const char *dot = ".";
	char *ret;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			F0->mat[linecounter]=atof(buffer);
		else
			F0->mat[linecounter]=atoi(buffer);
		if(DEBUG_FNLSDP && pp==0)printf("F0->mat[%d]=%f\n",linecounter, F0->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);

	/* read Q0 file*/
	strcat(filenameQ0,dir);
	strcat(filenameQ0,code);
	strcat(filenameQ0,"_");
	strcat(filenameQ0,id);
	strcat(filenameQ0,Q0str);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameQ0);
	if ((fp = fopen(filenameQ0, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameQ0);
		exit(1);
	}
	Q0->mat=(double *)malloc(sizeof(double)*(nx)*(nx));
	Q0->nrows=nx;
	Q0->ncols=nx;
	linecounter=0;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			Q0->mat[linecounter]=atof(buffer);
		else
			Q0->mat[linecounter]=atoi(buffer);
		if(DEBUG_FNLSDP && pp==0)printf("Q0->mat[%d]=%f\n",linecounter, Q0->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);

	/* read V0 file*/
	strcat(filenameV0,dir);
	strcat(filenameV0,code);
	strcat(filenameV0,"_");
	strcat(filenameV0,id);
	strcat(filenameV0,V0str);
	if(DEBUG_FNLSDP && pp==0)printf("%s\n",filenameV0);
	if ((fp = fopen(filenameV0, "r"))==NULL) {
		if(DEBUG_FNLSDP && pp==0)printf("Cannot open file %s.\n",filenameV0);
		exit(1);
	}
	V0->mat=(double *)malloc(sizeof(double)*(nx)*(nx));
	V0->nrows=nx;
	V0->ncols=nx;
	linecounter=0;
	while(fscanf(fp,"%s\n", buffer)!=EOF){
		if(DEBUG_FNLSDP && pp==0)printf("%s\n", buffer);
		if((ret=strstr(buffer,dot))!=NULL)
			V0->mat[linecounter]=atof(buffer);
		else
			V0->mat[linecounter]=atoi(buffer);
		if(DEBUG_FNLSDP && pp==0)printf("V0->mat[%d]=%f\n",linecounter, V0->mat[linecounter]);
		linecounter++;	
	}
	fclose(fp);	
}

void free_initial_point(genmatrix *F0, genmatrix *Q0, genmatrix *V0,int pp){
	free(F0->mat);
	free(Q0->mat);
	free(V0->mat);
}

void free_compleib(genmatrix *A, genmatrix *B1, genmatrix *B,genmatrix *C1,genmatrix *C,genmatrix *D11,genmatrix *D12,genmatrix *D21, int pp){
	free(A->mat);
	free(B1->mat);
	free(B->mat);
	free(C1->mat);
	free(C->mat);
	free(D11->mat);
	free(D12->mat);
	free(D21->mat);
}

