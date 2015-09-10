/* 
        Description: The ANSI C code of the PRVNS approach
        Programmer:  Wesklei Migliorini
        E-mail:      wesklei.m@gmail.com
        Date:	     01/09/2015
        Lisence:     Free
        Note:        The system was developed using Linux.
        To compile:  Type: make
        To run: ./algorithm input.in
 */

//Includes and defines/*{{{*/
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

#ifndef abss
#define abss(a)     (a<0 ? (-a) : a)
#endif

typedef int bool;
#define true 1
#define false 0
#define FAIL 0

//Testes Diversidade PRVNS
#define DEBUG  0 //1=> info, 2=> all
#define GRAFICO  0//1=> convergencia pop, 2=> k

#include <pthread.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mersenne.h"
#include "functions.c"
#include <fcntl.h>
#include <unistd.h>
#include <pthread.h>
#include <unistd.h> /* sleep() */
#include <string.h>/*}}}*/

//Solution structs/*{{{*/

typedef struct _VNS_ISLAND
{
	int id;
	int neighborhood_id;
	double **sol_migration; //sol_migration[MIGRATION_PERCENT]
	pVNS **pop; //pop[POP_SIZE + MIGRATION_PERCENT]	
	double *best;//[DIM]
	int best_index;
	double bestfo;
	
	pthread_mutex_t islands_mutex;
	pVNS *vns_params;
}pVNS_ISLAND;

typedef struct _VNS_SOLUTION
{
	double bestfo;         //best fo value
	int c_aval;		//avaliation count total
	int c_aval_best;	//alaviation count when found best			
	time_t stime;//start time
	time_t etime;//end time
	double t_total;//total time

}pVNS_SOLUTION;

//global parameters
typedef struct _VNS
{
	int P; /*The metric used, l1,l2,l_inf */
	int FUNCTION;
	//Problem definitions
	int DIM;    //number of problem variables
	int RUN;  /*Algorithm can run many times in order to see its robustness*/
	double *best; //[DIM];           //best solution found
	double bestfo;         //best fo value
	double LB; //lower bound of the variables
	double UB; //upper bound of the variables
	int KMAX,TMAX,AVAL_MAX;
	int RADII_T_FLAG;
	int know;
	double Q,RADII,RHO,EPSILON;
	double *r;
	int METHOD;
	int ECO_STEP;
	int EVO_STEP;
	int VNS_POP;
	int POP_INDEX;
	double DELTA_MUTATION;
	int G_MAX;	
	int MIGRATION_PERCENT;
	int MIGRATION_SIZE;//save the correct %
	int MIGRATION_INTERVAL;
	int ISLANDS;
	double PC; //probabilidade de crossover
	
	pVNS_SOLUTION solv;
}pVNS;/*}}}*/

//global variables for graphics only/*{{{*/
double **fo_geracoes; //fo_geracoes[RUN][Ger]; Ger < G_MAX
double *fo_mediaGeracoes; //fo_mediaGeracoes[Ger]; Ger < G_MAX
double **bestfo_geracoes; //bestfo_geracoes[Ger]; Ger < G_MAX
double *bestfo_mediaGeracoes; //bestfo_mediaGeracoes[Ger]; Ger < G_MAX
int execucao = -1;
int geracao=-1;/*}}}*/

//Functions declarations/*{{{*/
double randon( double inferior, double superior);
//VNS
bool radiiBetween(double *max, double *r, int know);
bool radiiLess(double *max, double *r, int know);
bool checkisLbUb(double *x, const int DIM,double lb, double ub);
void getRange(double lb, double ub, double *lower, double *upper, double new_value);
void shake(double *x, double *y, double *r, int know, const int DIM,double lb, double ub, int *cont, int RADII_T_FLAG, int p, int max_changes, double delta_percent);
bool lp(double* x, double* y, double* r, int know, const int DIM, int RADII_T_FLAG, int p);
bool lpInf(double* x, double* y, double* r, int know, const int DIM,int RADII_T_FLAG);
void neighborhoodChange(double* x, double* y, double *fx, const int DIM, int *k, const int FUNCTION, int *best_aval, int *aval, double *fy, int lb, int ub);
void *PRVNS(void *arg);

void *PPRVNS_Island(void *arg);//Parallel Populational Reduced VNS
void *PPRVNS_Master(void *arg);//MASTER Paralel Populational Reduced VNS
void evolucaoVizinhancaPRNS(int *r1, int *r2, int *p, pVNS *vns, int i, double *y, pVNS **x);
void randomIndexPRVNS(int *r1, int *r2, int *p, pVNS *vns, int i);
void trocaVizinhancaPRNS(pVNS *vns, int i, int *t, int *best_aval, int *best_index, double *bestfo, double *y, double *fy, pVNS **x);

//END VNS 
//aux functions
void printvector(double* x, const int DIM, char* seq);
void calculateStatistic(pVNS **sol, int maxr);
void freeArrays(int* POP_SIZE, double** pop, double* fo, double* best);
void AllocArrays(int* POP_SIZE, int* DIM, double*** pop, double** fo, double** best, double** U);
void grafico_linhas_x_y(char *data, char *xtitle, char *ytitle, char *title, char *legend, char *filename);
void grafico_duas_linhas(char *data1,char *data2, char *xtitle, char *ytitle, char *title,  char *legend1, char *legend2, char *filename);
void printProgress(double fx,int now,int total);
int getParametersVNS(FILE *file, pVNS *vns);
void showParametersVNS(pVNS *vns);
int f(double x1, int FUNCTION);/*}}}*/

double randon( double inferior, double superior)/*{{{*/
{
//  double aux = (float)inferior + ((superior - inferior)*rand()/(RAND_MAX+1.0));
   /* double aux = (double)inferior + ((superior - inferior)*((double)MT_randInt(RAND_MAX)/((double)RAND_MAX+1.0))); */
  double aux = (float)inferior + ((superior - inferior)*rand()/(RAND_MAX+1.0));

  /* printf("min=%f max=%f rand=%f\n",inferior,superior,aux); */
  return aux;
}/*}}}*/

/*Main program of the search algorithm*//*{{{*/
int main(int argc, char **argv)
{
	srand(time(NULL));
	MT_seed();

	//when teste flag is used
	if(argc > 0 && strcmp(argv[1],"teste") == 0)
	{
		//put teste functions her
		return 0;
	}
	/* Control Parameters of the search algorithm*/
	int FUNCTION;
	//Problem definitions
	int RUN;  /*Algorithm can run many times in order to see its robustness*/
	int i, j;
	double lb; //lower bound of the variables
	double ub; //upper bound of the variables

        pVNS *vns_run = (pVNS *) malloc (sizeof (pVNS));
	//read input file
	FILE *file = fopen( argv[1], "r" );
	if (getParametersVNS(file,vns_run) == -1 ){	//read input file
		return 0;
	}
	fclose(file);

	showParametersVNS(vns_run);

	RUN = vns_run->RUN;
	FUNCTION = vns_run->FUNCTION;
	prepararObjFunc(&FUNCTION,&lb,&ub);//initialize variables needed by function
	vns_run->LB = lb;
	vns_run->UB = ub;

	int maxr=RUN;

        pVNS **sol = (pVNS **) malloc (sizeof (pVNS*)*maxr);
	for(i=0;i<maxr;i++){
		sol[i] = (pVNS*) malloc (sizeof (pVNS));
		memcpy(sol[i],vns_run,sizeof(pVNS));
		sol[i]->RUN=i;

		sol[i]->best = malloc(sizeof(double) * vns_run->DIM);

		for (j=0; j<vns_run->DIM;j++) //each dimension
		{
			sol[i]->best[j] = randon(vns_run->LB,vns_run->UB);
		}
	}

	if(GRAFICO){//aloca para graficos/*{{{*/
		fo_geracoes = (double**) malloc(sizeof(double*) * RUN + 2);
		bestfo_geracoes = (double**) malloc(sizeof(double*) * RUN + 2);
		for(i=0;i<RUN;i++){
			if(vns_run->METHOD == 6) { //PRVNS
				fo_geracoes[i] = (double*) malloc(sizeof(double) * vns_run->G_MAX + 2 );
				bestfo_geracoes[i] = (double*) malloc(sizeof(double) * vns_run->G_MAX + 2);
			}else{
				fo_geracoes[i] = (double*) malloc(sizeof(double) * vns_run->TMAX + 2 );
				bestfo_geracoes[i] = (double*) malloc(sizeof(double) * vns_run->TMAX + 2);
			}
		}
			if(vns_run->METHOD == 6) { //PRVNS
				fo_mediaGeracoes = (double*) malloc(sizeof(double) * vns_run->G_MAX + 2);
				bestfo_mediaGeracoes = (double*) malloc(sizeof(double) * vns_run->G_MAX + 2);
			}else{
				fo_mediaGeracoes = (double*) malloc(sizeof(double) * vns_run->TMAX + 2);
				bestfo_mediaGeracoes = (double*) malloc(sizeof(double) * vns_run->TMAX + 2);
			}
	}	/*}}}*/


	//executa o metodo escolhido pelo arquivo de entrada
	double *r = (double*)malloc ( (vns_run->KMAX+1) *sizeof(double));
	double radii = vns_run->RADII;
	switch(vns_run->METHOD)
	{
		case 6:
		case 8:
			printf("Valores de raio usando PG e razao %f:\n",vns_run->Q);
			r[0] = 0.1f;
			r[1] = 0.3f;
			r[2] = 0.5f;
			r[3] = 0.7f;
			r[4] = 0.9f;
			for(j=0;j<vns_run->KMAX;j++){
				printf("k_%d => radii=%f \n",j+1,r[j]);
			}
			break;
	}

	for(i=0;i<maxr;i++){

		if(GRAFICO){/*{{{*/
			execucao++;
		}/*}}}*/

		switch (vns_run->METHOD) {

		case 6: //PRVNS
			sol[i] = PRVNS((void*)sol[i]);
			break;
		case 8: //PPRVNS

			sol[i] = PPRVNS_Master((void*)sol[i]);
			break;
			
		default:
			printf("Info: Invalid method.\n") ;
			exit(0);
		}

	}

	//GRAFICOS/*{{{*/
	//INICIO Grafico Convergencia inicializacao
#ifdef GRAFICO
	char *vns_plot_best;
	char *vns_plot_mean;
	char file_name[100];	

	if(GRAFICO == 1){
		vns_plot_best = (char *) malloc(sizeof(char) * vns_run->TMAX * 10);
		vns_plot_mean = (char *) malloc(sizeof(char) * vns_run->TMAX * 10);
		sprintf(file_name,"Convergencia/VNS_CONVERGENCIA_%d_%d_%d_%d",vns_run->METHOD,vns_run->DIM,vns_run->RUN,vns_run->FUNCTION);
	}
	//FIM Grafico Convergencia inicializacao

	if(GRAFICO == 1 && vns_run->METHOD == 4){//para RVNS
		for(i=0;i<=geracao;i++){

			double med_best = 0.0f;

			for(j=0;j<RUN;j++){
				med_best += bestfo_geracoes[j][i];
			}

			bestfo_mediaGeracoes[i] = med_best/RUN;

			//salva em string para gerar pelo gnuplot
			sprintf(vns_plot_best, "%s %d %lf\n ",vns_plot_best,i,bestfo_mediaGeracoes[i]);

			if(DEBUG){
				/* printf("\nMedia melhor execucoes Geracao[%d] =>%f\n",i,bestfo_mediaGeracoes[i]); */
			}
		}
		
		//escreve em arquivo
		grafico_linhas_x_y(vns_plot_best,"FO", "Geracoes", "RVNS", "Melhor",file_name);
	}else if(GRAFICO == 1 && vns_run->METHOD != 4){

		for(i=0;i<=geracao;i++){

			double med_med = 0.0f;
			double med_best = 0.0f;

			for(j=0;j<RUN;j++){
				med_med += fo_geracoes[j][i];
				med_best += bestfo_geracoes[j][i];
			}

			fo_mediaGeracoes[i] = med_med/RUN;
			bestfo_mediaGeracoes[i] = med_best/RUN;

			//salva em string para gerar pelo gnuplot
			sprintf(vns_plot_mean, "%s %d %lf\n ",vns_plot_mean,i,fo_mediaGeracoes[i]);
			sprintf(vns_plot_best, "%s %d %lf\n ",vns_plot_best,i,bestfo_mediaGeracoes[i]);

			if(DEBUG){
				printf("\nMedia execucoes Geracao[%d] =>%f\n",i,fo_mediaGeracoes[i]);
				printf("\nMedia melhor execucoes Geracao[%d] =>%f\n",i,bestfo_mediaGeracoes[i]);
			}
		}
		
		//escreve em arquivo
		grafico_duas_linhas(vns_plot_best,vns_plot_mean,"Geracoes","FO","PRVNS","Melhor","Media",file_name);
	}
#endif/*}}}*/

	//realiza os calculos necessarios
	calculateStatistic(sol,maxr);

	//liberar memoria
	for(i=0;i<maxr;i++){
		free(sol[i]);
	}
	free(sol);
	free(vns_run);

	pthread_exit(NULL);
	return 0;
}//end MAIN/*}}}*/

void calculateStatistic(pVNS **sol, int maxr){/*{{{*/
//faz os calculos necessarios, como media, desvio padrao, etc

	int i;
	double fo_best = sol[0]->solv.bestfo;
	double t_best = sol[0]->solv.t_total;
	int run_best = sol[0]->solv.c_aval;
	int run_when_best = sol[0]->solv.c_aval_best;
	double fo_sum = 0.0;
	double t_sum = 0.0;
	double run_sum = 0.0;
	double fo_mean = 0.0;
	double t_mean = 0.0;
	double run_mean = 0.0;
	double fo_vari = 0.0;
	double t_vari = 0.0;
	double run_vari = 0.0;
	double fo_devi = 0.0;
	double t_devi = 0.0;
	double run_devi = 0.0;
	double fo_sdevi = 0.0;
	double t_sdevi = 0.0;
	double run_sdevi = 0.0;

	for(i=0;i<maxr;i++){

		fo_sum += sol[i]->solv.bestfo;
		if(fo_best > sol[i]->solv.bestfo){
			fo_best = sol[i]->solv.bestfo;
		}

		t_sum += sol[i]->solv.t_total;
		if(t_best > sol[i]->solv.t_total){
			t_best = sol[i]->solv.t_total;
		}

		run_sum += sol[i]->solv.c_aval;
		if(run_best > sol[i]->solv.c_aval){
			run_best = sol[i]->solv.c_aval;
			run_when_best = sol[i]->solv.c_aval_best;
		}
	}

	fo_mean = fo_sum/maxr;
	t_mean = t_sum/maxr;
	run_mean = run_sum/maxr;

	for(i = 0; i < maxr; ++i){
		fo_devi += pow((sol[i]->solv.bestfo - fo_mean),2);
		t_devi += pow((sol[i]->solv.t_total - t_mean),2);
		run_devi += pow((sol[i]->solv.c_aval - run_mean),2);
	}

	fo_vari = fo_devi / (maxr);
	t_vari = t_devi / (maxr);
	run_vari = run_devi / (maxr);

	fo_sdevi = sqrt(fo_vari);
	t_sdevi = sqrt(t_vari);
	run_sdevi = sqrt(run_vari);

	printf("FO Standard Deviation: %6.10f\n", fo_sdevi);
	printf("Time Standard Deviation: %f\n", t_sdevi);
	printf("Iteration Standard Deviation: %6.3f\n", run_sdevi);
	printf("FO mean: %6.10f\n", fo_mean);
	printf("Time mean: %f\n", t_mean);
	printf("Iteration mean: %6.3f\n", run_mean);
	printf("FO best: %6.10f\n", fo_best);
	printf("Time best: %f\n", t_best);
	printf("Iteration count when found best: %d\n", run_when_best);
	printf("Iteration best: %d\n", run_best);

	printf("End running\n");
}/*}}}*/

void AllocArrays(int* POP_SIZE, int* DIM, double*** pop, double** fo, double** best, double** U)/* Dynamic array allocation *//*{{{*/
/* mol and sites */
{
	int j;
	double** aux = malloc ((*POP_SIZE)*sizeof(double*));

        for (j = 0;j < *POP_SIZE; j++){
		aux[j] = (double*)malloc (*DIM * sizeof(double));
	}

	*pop = aux;
	*fo = (double*)malloc (*POP_SIZE*sizeof(double));
	*best = (double*)malloc (*DIM*sizeof(double));
	*U = (double*)malloc(*DIM * sizeof(double));
}/*}}}*/

void freeArrays(int* POP_SIZE, double** pop, double* fo, double* best)/* Free arrays *//*{{{*/
{
	int i;
	free(fo);
	free(best);
	
	for (i = 0; i < *POP_SIZE; i++){
		free(pop[i]);
	}

	free(pop);
}/*}}}*/

int getParametersVNS(FILE *file, pVNS *vns)/*{{{*/
{
	int *P,*ECO_STEP,*EVO_STEP,*VNS_POP,*G_MAX,*FUNCTION,*DIM,*RUN;
	int *KMAX,*TMAX,*AVAL_MAX,*METHOD;
	int *RADII_T_FLAG;
	double *Q,*RADII,*RHO,*EPSILON,*PC;

	FUNCTION = (int*) malloc(sizeof(int));
	DIM = (int*) malloc(sizeof(int));
	RUN = (int*) malloc(sizeof(int));
	ECO_STEP = (int*) malloc(sizeof(int));
	EVO_STEP = (int*) malloc(sizeof(int));
	VNS_POP = (int*) malloc(sizeof(int));
	KMAX = (int*) malloc(sizeof(int));
	TMAX = (int*) malloc(sizeof(int));
	METHOD = (int*) malloc(sizeof(int));
	AVAL_MAX = (int*) malloc(sizeof(int));
	RADII_T_FLAG = (int*) malloc(sizeof(int));
	Q = (double*) malloc(sizeof(double));
	RADII = (double*) malloc(sizeof(double));
	RHO = (double*) malloc(sizeof(double));
	EPSILON = (double*) malloc(sizeof(double));
	P = (int*) malloc(sizeof(int));
	G_MAX = (int*) malloc(sizeof(int));
	PC = (double*) malloc(sizeof(double));

	if (file == 0)
	{
		printf( "Could not open ini file! Usage ./<exec> <file.in>\n" );
		return -1;
	}
	else 
	{
		ffscanf("RUN", file, "%d", &RUN);
		ffscanf("DIM", file, "%d", &DIM); 
		ffscanf("FUNCTION",file, "%d", &FUNCTION);
		ffscanf("KMAX", file, "%d", &KMAX);
		ffscanf("TMAX", file, "%d", &TMAX);
		ffscanf("Q", file, "%lf", &Q);
		ffscanf("RADII", file, "%lf", &RADII); 
		ffscanf("RADII_T_FLAG", file, "%d", &RADII_T_FLAG); 
		ffscanf("AVAL_MAX", file, "%d", &AVAL_MAX); 
		ffscanf("RHO", file, "%lf", &RHO); 
		ffscanf("EPSILON", file, "%lf", &EPSILON); 
		ffscanf("METHOD", file, "%d", &METHOD); 
		ffscanf("P", file, "%d", &P); 
		ffscanf("VNS_POP", file, "%d", &VNS_POP); 
		ffscanf("ECO_STEP", file, "%d", &ECO_STEP); 
		ffscanf("EVO_STEP", file, "%d", &EVO_STEP); 
		ffscanf("G_MAX", file, "%d", &G_MAX); 
		ffscanf("PC", file, "%lf", &PC); 

		vns->RUN = *RUN;
		vns->DIM = *DIM;
		vns->FUNCTION = *FUNCTION;
		vns->KMAX = *KMAX;
		vns->TMAX = *TMAX;
		vns->Q = *Q;
		vns->RADII = *RADII;
		vns->RADII_T_FLAG = *RADII_T_FLAG;
		vns->AVAL_MAX = *AVAL_MAX;
		vns->RHO = *RHO;
		vns->EPSILON = *EPSILON;
		vns->METHOD = *METHOD;
		vns->P = *P;
		vns->VNS_POP = *VNS_POP;
		vns->ECO_STEP = *ECO_STEP;
		vns->EVO_STEP = *EVO_STEP;
		vns->G_MAX = *G_MAX;
		vns->PC = *PC;

		return 1;
	}
}/*}}}*/

void showParametersVNS(pVNS *vns)/*{{{*/
{
	printf("***VNS PARAMETERS***\n");
	printf("RUNS = %d\n", vns->RUN);
	printf("DIM = %d\n", vns->DIM);
	printf("FUNCTION = %s\n", getFunctionName(vns->FUNCTION));
	printf("METHOD = %d\n", vns->METHOD);
	printf("****************\n");
	printf("KMAX = %d\n",vns->KMAX);
	printf("TMAX = %d\n", vns->TMAX);
	printf("Q = %f\n", vns->Q);
	printf("RADII = %f\n", vns->RADII);
	printf("RADII_T_FLAG = %d\n", vns->RADII_T_FLAG);
	printf("METRIC P = %d\n", vns->P);
	printf("VNS_POP = %d\n", vns->VNS_POP);
	printf("ECO_STEP = %d\n", vns->ECO_STEP);
	printf("EVO_STEP = %d\n", vns->EVO_STEP);
	printf("****************\n");
	printf("***HOOKE AND JEEVES PARAMETERS***\n");
	printf("AVAL_MAX = %d\n", vns->AVAL_MAX);
	printf("RHO = %f\n", vns->RHO);
	printf("EPSILON = %f\n", vns->EPSILON);
	printf("****************\n");
}/*}}}*/

int ffscanf(char *fieldname, FILE *fp, char *format, void** inbuffer)/*{{{*/
{
	char buffer[100];
	int len;
	int commentflag = 0;
	char *pch;
	char *pch2;
	do
	{
		if(fgets(buffer, 99, fp) == 0){
		       return FAIL;
		}

		buffer[99] = '\0';
		len = strlen(buffer);

		if (buffer[len - 1] == '\n'){
		       buffer[len - 1] = '\0';
		}

		switch (commentflag) 
		{
			case 0:
				
				if (strstr(buffer, "/*") != 0) 
				{
					commentflag = 1;
					if (strstr(buffer, "*/") != 0){
						commentflag = 2;
					}
				}

				break;
			case 1:

				if (strstr(buffer, "*/") != 0){
					commentflag = 2;
				}

				break;
			    case 2:

				if (strstr(buffer, "/*") != 0) 
				{
					commentflag = 1;
					if (strstr(buffer, "*/") != 0){
						commentflag = 2;
					}
				}
				else{
					commentflag = 0;
				}

				break;
		}	

	}while(commentflag != 0);	

	//separate field name: token = "="
	if (strstr (buffer, fieldname) != 0)
	{
		pch = strtok (buffer,"=");

		pch2 = pch;
		while (pch != NULL)
		{
			pch2 = pch;
			pch = strtok (NULL, "= ");
		}

		sscanf(pch2, format, *inbuffer);
		return 1;//ok
	}
	else{
	       return 0; 
	}
}/*}}}*/

void printProgress(double fx,int now,int total){/*{{{*/
	int i;
	double percent = (double)(100*now)/total;
	int i_percent = (int)percent;
	for(i=0;i<i_percent;i+=10){
			printf("=");
	}
	printf("> %.2f%% %.20f\r" ,percent,fx);
}/*}}}*/

void printvector(double* x, const int DIM, char* seq){/*{{{*/
	int i;
	for(i=0;i<DIM;i++){
		printf(" %s[%d]=%f",seq,i,x[i]);
	}
	printf("\n");
}/*}}}*/

void neighborhoodChange(double* x, double* y, double *fx, const int DIM, int *k, const int FUNCTION, int *best_aval, int *aval, double *fy, int lb, int ub){ /*{{{*/

	checkisLbUb(y,DIM,lb,ub);

	*fy = objfunc(y,&FUNCTION,&DIM,aval);
	if( *fy < *fx){
		*best_aval = *aval;
		memcpy(x,y,DIM*sizeof(double));//set as start value
		*fx = *fy;
		*k = 1;

	}else{
		*k += 1;
	}
}/*}}}*/

void getRange(double lb, double ub,double *lower, double *upper, double new_value){/*{{{*/

	if(ub > new_value){
		*upper = new_value;
	}else{
		*upper = ub;
	}

	if(lb < new_value){
		*lower = new_value;
	}else{
		*lower = lb;
	}
}/*}}}*/

void shake(double *x, double *y, double *r, int know, const int DIM,double lb, double ub, int *cont, int RADII_T_FLAG, int p, int max_changes, double delta_percent){/*{{{*/
/**
 * max changes == 1, is shake one
 */ 
	*cont += 1;
	int i,cont_changes=0;;
	double aux;
	int itoChange = randon(0,DIM);
	memcpy(y,x,DIM*sizeof(double));
	
	do
	{
		double amount_ratio = r[know];

		for(i=0;i<DIM;i++){
			if(i == itoChange || randon(0,1) < delta_percent){
				if(cont_changes == max_changes){
					break;		
				}
				cont_changes++;

				aux = y[i];
				if(randon(0,1) > 0.5f){			
					while(aux == y[i] || (aux < lb && aux > ub)) {
						aux = randon(aux,aux+(amount_ratio));
					}
				}else{
					while(aux == y[i] || (aux < lb && aux > ub)){
						aux = randon(aux-(amount_ratio),aux);
					}
				}

				if ( aux < lb){
					aux = lb;
					
				}else if(aux > ub ){
					aux = ub;
				}

				y[i] = aux;
				amount_ratio = amount_ratio - fabs(fabs(x[i]) - fabs(aux));

				if(amount_ratio < 0.00001) {
					i=DIM;
					break;
				}
			}
		}
	}while(!lp(x,y,r,know,DIM,RADII_T_FLAG,p));
}/*}}}*/

bool lp(double* x, double* y, double* r, int know, const int DIM, int RADII_T_FLAG, int p){/*{{{*/

	int i;
	double sum = 0;
	double dist,aux;
	if(p == 0){
		dist = fabs ((double) fabs(x[0])-fabs(y[0]));
		for(i=1;i<DIM;i++){
			aux = fabs ((double) fabs(x[i])- fabs(y[i]));

			if(dist < aux){
				dist = aux;
			}
		}
	}
	else
	{
		for(i=0;i<DIM;i++){
			sum += pow( fabs (x[i]-y[i]) , p );
		}

		if(p == 2){
			dist = sqrt(sum);
		}else{
			dist = sum;
		}
	}

	//check de distance, if is between the correct radii value
	switch (RADII_T_FLAG) {
		case 0: 
			return radiiBetween(&dist,r,know);	

			break;
		case 1: 
			return radiiLess(&dist,r,know);	

			break;
	}

	return false; 
}/*}}}*/

bool radiiBetween(double *dist, double *r, int know){/*{{{*/

	int k_bef = know;
	k_bef--;

	if(*dist > r[k_bef] && *dist < r[know]){
	       return true;
	}else{
 	       return false;
	}		
}/*}}}*/

bool radiiLess(double *dist, double *r, int know){/*{{{*/


	if(*dist < r[know]){
	       return true;
	}else{
 	       return false;
	}		
}/*}}}*/

//hooke use
bool checkisLbUb(double *x, const int DIM,double lb, double ub){/*{{{*/
	
	int i=0;
	bool flag = false;
	for(;i<DIM;i++){
		if(x[i] < lb){
			flag = true;
		       x[i] = lb;
		}else if(x[i] > ub){
			flag = true;
		       x[i] = ub;
		}
		/* if(x[i] < lb || x[i] > ub){ */
		/* 	flag = false; */
		/* 	break; */
		/* } */
	}
	return flag;
}/*}}}*/

void *PRVNS(void *arg){//Populational Reduced VNS/*{{{*/

	pVNS *vns = arg;

	int i,j, t=0,iter=0;
	int r1,r2,p, best_aval, best_index;

	double *y = malloc(sizeof(double) * vns->DIM);

	double fy;
	double bestfo;

	//GRAFICOS
	//INICIO Grafico Convergencia inicializacao
#ifdef GRAFICO
	char *vns_plot_k;
	/* char *vns_plot_best; */
	double sum=0;
	char file_name[100];	

	if(GRAFICO == 1){
		/* vns_plot_best = (char *) malloc(sizeof(char) * vns->TMAX * 1000); */
		sprintf(file_name,"Convergencia/VNS_%d_%d_%d_%d",vns->METHOD,vns->DIM,vns->FUNCTION,vns->RUN);
	}else if(GRAFICO == 2){
		vns_plot_k = (char *) malloc(sizeof(char) * vns->TMAX * 1000);
		sprintf(file_name,"Convergencia/VNS_K_%d_%d_%d_%d",vns->METHOD,vns->DIM,vns->FUNCTION,vns->RUN);
	}
#endif
	//FIM Grafico Convergencia inicializacao
	

	//definir raio
	vns->r = malloc(sizeof(double) * vns->KMAX);
	vns->r[0] = 0.1f;
	vns->r[1] = 0.3f;
	vns->r[2] = 0.5f;
	vns->r[3] = 0.7f;
	vns->r[4] = 0.9f;

	//inicializar populacao
	pVNS **x = malloc (sizeof (pVNS*) * (vns->VNS_POP +1) );

	//faz o primeira antes para pegar o melhor
	x[0] =  malloc (sizeof (pVNS));

	x[0]->know = 0;
	x[0]->best = malloc(sizeof(double) * vns->DIM);
	//the start point
	for (j=0; j<vns->DIM;j++) //each dimension
	{
		x[0]->best[j] = randon(vns->LB,vns->UB);
	}

	//avaliar individuo
	x[0]->bestfo = objfunc(x[0]->best,&vns->FUNCTION,&vns->DIM,&t);


	//inicializa o melhor
	bestfo = x[0]->bestfo;
	best_aval = t;
	best_index = 0;


	for(i=1;i<vns->VNS_POP;i++){
		x[i] =  malloc (sizeof (pVNS));

		x[i]->know = 0;
		x[i]->best = malloc(sizeof(double) * vns->DIM);
		//the start point
		for (j=0; j<vns->DIM;j++) //each dimension
		{
			x[i]->best[j] = randon(vns->LB,vns->UB);
		}

		//evaluate indivi
		x[i]->bestfo = objfunc(x[i]->best,&vns->FUNCTION,&vns->DIM,&t);


		if(x[i]->bestfo < bestfo){
			bestfo = x[i]->bestfo;
			best_aval = t;
			best_index = i;
		}
	}


#ifdef DEBUG
	if(DEBUG == 2){
		printf("\nGeracao=> %d\n",iter);
		for(i=0;i<vns->VNS_POP;i++){
			printf("x[%d]=> [",i);
			for(j=0;j<vns->DIM;j++){
				printf(" %f",x[i]->best[j]);

				if(j < vns->DIM -1) printf(","); 
			}
			printf("] FO=> %f", x[i]->bestfo);

			if(best_index == i){
				printf("*\n");
			}else{
				printf("\n");
			}
		}
	}
#endif

	//GRAFICOS
#ifdef GRAFICO
	if(GRAFICO == 1){
		geracao = 0;
		
		//melhor fo
		bestfo_geracoes[execucao][geracao] = bestfo;

		/* //calcular a media para gerar grafico */
		fo_geracoes[execucao][geracao] = 0.0f;
		for (j=0; j<vns->VNS_POP; j++)
		{
			fo_geracoes[execucao][geracao] += x[j]->bestfo;
		}
		fo_geracoes[execucao][geracao] = fo_geracoes[execucao][geracao]/(double)vns->VNS_POP;

		if(DEBUG){
			printf("\nMedia Geracao[%d] =>%f\n",geracao,fo_geracoes[execucao][geracao]);
			printf("\nMelhor[%d] =>%f\n",geracao,bestfo_geracoes[execucao][geracao]);
		}

	}else if(GRAFICO == 2){
		//calcular a media de K para gerar grafico
		sum = 0;
		for(j=0;j<vns->VNS_POP;j++){
			sum += x[j]->know + 1;
		}
		sprintf(vns_plot_k, "%s %d %lf\n ",vns_plot_k,iter,(double)sum/vns->VNS_POP);
	}
#endif
	
	

    	time (&(vns->solv.stime));
	while((iter + 1 < vns->G_MAX) && ((t + vns->VNS_POP) <= vns->TMAX)){

		
#ifdef GRAFICO
		if(GRAFICO){
			geracao++;
		}
#endif
		

		//etapa populacional
		for(i =0; i < vns->VNS_POP;i++){

			//usar os whiles separados para tentar melhorar o tempo
			//seleciona um por vez pode ser melhor que selecionar todos
			//e testar no mesmo while
			do{
				r1 = randon(0,vns->VNS_POP);
			}while(r1 == i);

			do{
				r2 = randon(0,vns->VNS_POP);
			}while( r1 == r2 || r2 == i );
			

			//pos na dim para sempre perturbar
			p = randon(0,vns->DIM);

			
#ifdef DEBUG
			if(DEBUG == 2){
				printf("r1=> %d r2=> %d p=> %d\n",r1,r2,p);
			}
#endif
			

			//etapa da evolucao diferencial
			for(j=0;j<vns->DIM;j++){

				if( (j == p) || (randon(0,1) < vns->PC)){

					y[j] = x[r2]->best[j] + ( randon(-1.0f * vns->r[x[i]->know], vns->r[x[i]->know]) * x[r1]->best[j] ); 
					/* y[j] = x[i]->best[j] + ( randon(-1.0f * vns->r[x[i]->know], vns->r[x[i]->know]) * x[r1]->best[j] );  */
				}else{
					y[j] = x[i]->best[j];
				}	
				
				//trata se esta no intervalo permitido da FO
				if(vns->UB < y[j]){
					y[j] = vns->UB;
				}else if(vns->LB > y[j]){
					y[j] = vns->LB;
				}
			}

			//etapa de avaliacao e troca de vizinhanca
			fy = objfunc(y,&vns->FUNCTION,&vns->DIM,&t);

			
#ifdef DEBUG
			if(DEBUG == 2){
				printf("Depois:\n");
				printf("x[%d]=> [",i);
				for(j=0;j<vns->DIM;j++){
					printf(" %f",y[j]);

					if(j < vns->DIM -1) printf(","); 
				}
				printf("] FO=> %f", fy);

				if(best_index == i){
					printf("*\n");
				}else{
					printf("\n");
				}
			}
#endif

			if(fy < x[i]->bestfo){
				
				memcpy(x[i]->best,y,sizeof(double) * vns->DIM);
				/* x[i]->best = y; */
				x[i]->bestfo = fy;

				x[i]->know = 0;


				if(x[i]->bestfo < bestfo){
					bestfo = x[i]->bestfo;
					best_aval = t;
					best_index = i;
				}
			}else{

				if(x[i]->know <= vns->KMAX -1){
					x[i]->know += 1;
				}else{
					x[i]->know = 0; //reinicia k
				}
			}
		}
		iter++;

			
#ifdef DEBUG
		if(DEBUG == 2){
			printf("\nGeracao=> %d\n",iter);
			for(i=0;i<vns->VNS_POP;i++){
				printf("x[%d]=> [",i);
				for(j=0;j<vns->DIM;j++){
					printf(" %f",x[i]->best[j]);

					if(j < vns->DIM -1) printf(","); 
				}
				printf("] FO=> %f", x[i]->bestfo);

				if(best_index == i){
					printf("*\n");
				}else{
					printf("\n");
				}
			}
		}
#endif

#ifdef GRAFICO
		//GRAFICOS
		if(GRAFICO == 1){

			//melhor fo
			bestfo_geracoes[execucao][geracao] = bestfo;
			
			/* //calcular a media para gerar grafico */
			fo_geracoes[execucao][geracao] = 0.0f;
			for (j=0; j<vns->VNS_POP; j++)
			{
				fo_geracoes[execucao][geracao] += x[j]->bestfo;
			}
			fo_geracoes[execucao][geracao] = fo_geracoes[execucao][geracao]/(double)vns->VNS_POP;

			if(DEBUG){
				printf("\nMedia Geracao[%d] =>%f\n",geracao,fo_geracoes[execucao][geracao]);
				printf("\nMelhor Geracao[%d] =>%f\n",geracao,bestfo_geracoes[execucao][geracao]);
			}


		}else if(GRAFICO == 2){
			//calcular a media de K para gerar grafico
			sum = 0;
			for(j=0;j<vns->VNS_POP;j++){
				sum += x[j]->know + 1;
			}
			sprintf(vns_plot_k, "%s %d %lf\n ",vns_plot_k,iter,(double)sum/vns->VNS_POP);
		}
#endif

	}
    	time (&(vns->solv.etime));
	vns->solv.t_total = difftime(vns->solv.etime,vns->solv.stime);

	printf("\n==RUN: %d\n",vns->RUN);
	/* printf("VNS total number of avaliations : %d\n",vns_cont); */
	/* printf("Hooke total number of avaliations: %d\n",hooke_geral_cont); */
	printf("Total number of avaliation: %d\n",t);
	printf("Best solution found == %g\n",bestfo);
	printf("Time: == %.0f seconds\n",vns->solv.t_total);

	vns->solv.c_aval=t;
	vns->solv.c_aval_best=best_aval;
	vns->solv.bestfo=bestfo;

	printf("\n");
	
	/* free(f_name); */
	/* free(vns_plot); */

	/* free(x); */
	/* free(y); */
	/* free(r); */


#ifdef GRAFICO
	if(GRAFICO == 2){
		grafico_linhas_x_y(vns_plot_k,"Geracoes","K","PRVNS","K",file_name);
	}
	/* grafico_linhas_x_y(vns_plot_best,"Geracoes","FO","PRVNS","PRVNS",file_name); */
#endif
	
	//when return, return the best of the best
	return (void*)arg;
}/*}}}*/

//PPRVNS START/*{{{*/

void *PPRVNS_Master(void *arg){//Populational Reduced VNS/*{{{*/

/*
typedef struct _VNS_ISLAND
{
	int id;
	int neighborhood_id;
	double **sol_migration; //sol_migration[MIGRATION_PERCENT]
	double **pop; //pop[POP_SIZE + MIGRATION_PERCENT]
	double *pop_fo; //pop[POP_SIZE + MIGRATION_PERCENT]
	double *best;
	int best_index;
	double bestfo;
}pVNS_ISLAND;
*/
	pVNS *vns = arg;
	int i,j,t,rc;

	int MIGRATION_INTERVAL = 1;
	int MIGRATION_PERCENT=5;
	int ISLANDS = 3;
	pthread_t islands_thread[ISLANDS];
	void *status;
    
    vns->MIGRATION_SIZE = (vns->MIGRATION_PERCENT * 100)/vns->VNS_POP +1;//round up
    vns->islands_mutex = (pthread_mutex_t *) malloc (sizeof (pthread_mutex_t)*(vns->ISLANDS);
    vns->islands = (pVNS **) malloc (sizeof (pVNS*)*vns->ISLANDS);
    
    pVNS islands_sol[ISLANDS];// = (pVNS **) malloc (sizeof (pVNS*)*ISLANDS);
			
	//definir raio
	vns->r = malloc(sizeof(double) * vns->KMAX);
	vns->r[0] = 0.1f;
	vns->r[1] = 0.3f;
	vns->r[2] = 0.5f;
	vns->r[3] = 0.7f;
	vns->r[4] = 0.9f;
	vns->know = 0;
	
	//TODO randon start on thread
	for(t=0; t<ISLANDS; t++){
		islands_sol[t] = (pVNS *) malloc (sizeof (pVNS));

		islands_sol[t].vns_params = (pVNS *) malloc (sizeof (pVNS));
		memccpy(vns->islands[t].vns_params,vns,sizeof(pVNS));

		islands_sol[t].id = t;
		islands_sol[t].neighborhood_id = t+1;
		
		islands_sol[t].sol_migration = (double **) malloc (sizeof (double*)*vns->MIGRATION_SIZE);
		for(i=0; i<vns->MIGRATION_SIZE; i++){
			islands_sol[t].sol_migration[i] = (double *) malloc (sizeof (double)*vns->DIM);
		}
		
		islands_sol[t].pop = (pVNS **) malloc (sizeof (pVNS*)*vns->VNS_POP);
		for(i=0; i<vns->VNS_POP; i++){
			islands_sol[t].pop[i] = malloc (sizeof (pVNS));
		}
				
		islands_sol[t].pop[i]->best = malloc(sizeof(double) * vns->DIM);
		//the start point
		for (j=0; j<vns->DIM;j++) //each dimension
		{//TODO continue from her, adjust to pop[i]
			islands_sol[t].pop[i]->best[j] = randon(vns->LB,vns->UB);
		}

		//avaliar individuo
		islands_sol[t].bestfo = objfunc(islands_sol[t].best,&vns->FUNCTION,&vns->DIM,&t);
		islands_sol[t].best_index = 0;

		for (j=1; j<vns->VNS_POP;j++) //each dimension
		{
				for(i=0;i<vns->DIM;i++){
					islands_sol[t].pop[j][i] = randon(vns->LB,vns->UB);
				}

				//evaluate indivi
				islands_sol[t].pop_fo[j] = objfunc(islands_sol[t].pop[j],&vns->FUNCTION,&vns->DIM,&t);


				if(islands_sol[t].pop_fo[j] < bestfo){
					islands_sol[t].bestfo = islands_sol[t].pop_fo[j];
					islands_sol[t].best_index = j;
				}
		}
	}
	islands_sol[ISLANDS-1].neighborhood_id = 0;//the last points to first

	/* for(i=0;i<MIGRATION_INTERVAL;i++){ */
		for(t=0; t<ISLANDS; t++){

			printf("In PPRVNS_Master: creating thread %d\n", t);
			rc = pthread_create(&islands_thread[t], NULL, PPRVNS_Island, (void *)&islands_sol[t]);
			if (rc){
				printf("ERROR; return code from pthread_create() is %d\n", rc);
				exit(-1);
			}
		}

		for(t=0; t<ISLANDS; t++) {
			printf("joining tread %d\n",t);
			rc = pthread_join(islands_thread[t], NULL);
			if (rc) {
				printf("ERROR; return code from pthread_join() is %d\n", rc);
				exit(-1);
			}
			
			printf("Main: completed join with thread %d having a best of %g\n",t, islands_sol[t].solv.bestfo);
		}

	/* } */

	//TODO return only the best
	return vns;
}/*}}}*/

void randomIndexPRVNS(int *r1, int *r2, int *p, pVNS *vns, int i){/*{{{*/
	//usar os whiles separados para tentar melhorar o tempo
			//seleciona um por vez pode ser melhor que selecionar todos
			//e testar no mesmo while
			do{
				*r1 = randon(0,vns->VNS_POP);
			}while(*r1 == i);

			do{
				*r2 = randon(0,vns->VNS_POP);
			}while( *r1 == *r2 || *r2 == i );
			

			//pos na dim para sempre perturbar
			*p = randon(0,vns->DIM);
}/*}}}*/

void evolucaoVizinhancaPRNS(int *r1, int *r2, int *p, pVNS *vns, int i, double *y, pVNS **x){/*{{{*/
	int j;
	//etapa da evolucao diferencial
	for(j=0;j<vns->DIM;j++){

		if( (j == *p) || (randon(0,1) < vns->PC)){

			y[j] = x[*r2]->best[j] + ( randon(-1.0f * vns->r[x[i]->know], vns->r[x[i]->know]) * x[*r1]->best[j] ); 
			/* y[j] = x[i]->best[j] + ( randon(-1.0f * vns->r[x[i]->know], vns->r[x[i]->know]) * x[r1]->best[j] );  */
		}else{
			y[j] = x[i]->best[j];
		}	

		//trata se esta no intervalo permitido da FO
		if(vns->UB < y[j]){
			y[j] = vns->UB;
		}else if(vns->LB > y[j]){
			y[j] = vns->LB;
		}
	}
}/*}}}*/

void trocaVizinhancaPRNS(pVNS *vns, int i, int *t, int *best_aval, int *best_index, double *bestfo, double *y, double *fy, pVNS **x){/*{{{*/

	//etapa de avaliacao e troca de vizinhanca
	*fy = objfunc(y,&vns->FUNCTION,&vns->DIM,t);

	if(*fy < x[i]->bestfo){

		memcpy(x[i]->best,y,sizeof(double) * vns->DIM);
		/* x[i]->best = y; */
		x[i]->bestfo = *fy;

		x[i]->know = 0;


		if(x[i]->bestfo < *bestfo){
			*bestfo = x[i]->bestfo;
			*best_aval = *t;
			*best_index = i;
		}
	}else{

		if(x[i]->know <= vns->KMAX -1){
			x[i]->know += 1;
		}else{
			x[i]->know = 0; //reinicia k
		}
	}
}/*}}}*/

void etapaPopulacionalPRNS(pVNS *vns,int *t, int *best_aval, int *best_index, double *bestfo, double *y, double *fy, pVNS **x){/*{{{*/
	int i;
	int r1,r2,p;
	r1 = 0;
	r2=0;
	p=0;
	//etapa populacional
	for(i =0; i < vns->VNS_POP;i++){
		randomIndexPRVNS(&r1,&r2,&p,vns,i);
		evolucaoVizinhancaPRNS(&r1,&r2,&p,vns,i,y,x);
		trocaVizinhancaPRNS(vns,i,t, best_aval,best_index, bestfo,y, fy, x);

	}
}/*}}}*/

void *PPRVNS_Island(void *arg){//Populational Reduced VNS/*{{{*/

	pVNS_ISLAND *vns_island = arg;
	pVNS vns = vns_island.vns;
	
	int i,j, t=0,iter=0;
	int  best_aval, best_index;

	double *y = malloc(sizeof(double) * vns->DIM);

	double fy;
	double bestfo;
	
    time (&(vns->solv.stime));
	while((iter + 1 < vns->G_MAX) && ((t + vns->VNS_POP) <= vns->TMAX)){
		
		etapaPopulacionalPRNS(vns,&t, &best_aval, &best_index, &bestfo, y, &fy, vns_island->pop);

		iter++;
		
		//TODO MIGRATION
		/*
		   Lock a mutex prior to updating the value in the shared
		   structure, and unlock it upon updating.
		   
		   pthread_mutex_lock (&mutex[neigh_id]);
		   pthread_mutex_unlock (&mutex[neigh_id]);
		 */
		
	}
    time (&(vns->solv.etime));
	vns->solv.t_total = difftime(vns->solv.etime,vns->solv.stime);

	printf("\n==RUN: %d\n",vns->RUN);
	/* printf("VNS total number of avaliations : %d\n",vns_cont); */
	/* printf("Hooke total number of avaliations: %d\n",hooke_geral_cont); */
	printf("Total number of avaliation: %d\n",t);
	printf("Best solution found == %g\n",bestfo);
	printf("Time: == %.0f seconds\n",vns->solv.t_total);

	vns->solv.c_aval=t;
	vns->solv.c_aval_best=best_aval;
	vns->solv.bestfo=bestfo;

	//TODO salvar em island struct the best

	printf("\n");
	
	/* free(f_name); */
	/* free(vns_plot); */

	/* free(x); */
	/* free(y); */
	/* free(r); */

	//when return, return the best of the best

	pthread_exit((void*) arg);
	/* return (void*)arg; */
}/*}}}*/

//PPRVNS END/*}}}*/

void grafico_linhas_x_y(char *data, char *xtitle, char *ytitle, char *title, char *legend, char *filename){/*{{{*/
// Funcao para gerar um grafico de uma linha (apenas uma mesmo) relacionando x/y
// Junto ao grafico eh salvo um arquivo filename.data com os dados usados para montar o grafico
//

//data: 	string com os dados na forma: "x y \n"
//		Deve ser separado por 'espaco/tab' seguindo a ordem: "valor_para_x <espaco> valor_para_y \n"
//xtitle: 	titulo do eixo X
//ytitle:	titulo do eixo Y
//title:	titulo do grafico
//legend:	legenda do grafico
//filename:	nome do arquivo de saida (o grafico) sem extensao. Caso precisar, pode passar caminhos relativos ou completos
//		aceita caminhos para subdiretorios. Ex graficos/meugrafico

	FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
	if(!gnuplotPipe){
	       printf("ERROR in gnuplotPipe!\n");
	       return;
	}

	char filename_data[255];
	sprintf(filename_data,"%s.data",filename);
	FILE *f_data = fopen(filename_data, "w");
	fputs(data,f_data);
	fflush(f_data);
	fclose(f_data);

	fputs(" set bmargin 7 \n",gnuplotPipe); 
	fputs(" unset colorbox \n",gnuplotPipe);
	fputs(" set terminal png enhanced font 'Verdana,10' \n",gnuplotPipe);
	fputs(" set output '",gnuplotPipe);
	fputs(filename,gnuplotPipe);
	fputs(".png' \n",gnuplotPipe );
	fputs(" set style line 2  lc rgb '#0404E9' lt 1 pt 7 \n",gnuplotPipe); //blue

	char paramaxis[1024];
	sprintf(paramaxis," set xlabel '%s' \n set ylabel '%s' \n",xtitle,ytitle);
	fputs(paramaxis,gnuplotPipe);

	char paramtitle[1024];
	sprintf(paramtitle," set title '%s'\n ",title);
	fputs(paramtitle,gnuplotPipe);

	char paramlegend[1024];
	sprintf(paramlegend," plot '-' title '%s' ls 2 with lines\n",legend);
	fputs(paramlegend,gnuplotPipe);
	fputs( data,gnuplotPipe);
	fputs("e",gnuplotPipe);
	fputs("\nquit",gnuplotPipe);
	fflush(gnuplotPipe);
	fclose(gnuplotPipe);
}/*}}}*/

void grafico_duas_linhas(char *data1,char *data2, char *xtitle, char *ytitle, char *title,  char *legend1, char *legend2, char *filename){/*{{{*/
// Funcao para gerar um grafico de duas linhas relacionando x/y
// Junto ao grafico eh salvo um arquivo filename_legend1.data e filename_legend2.data com os dados usados para montar o grafico
//

//data1 e data2: 	string com os dados na forma: "x y \n"
//			Deve ser separado por 'espaco/tab' seguindo a ordem: "valor_para_x <espaco> valor_para_y \n"
//xtitle: 		titulo do eixo X
//ytitle:		titulo do eixo Y
//title:		titulo do grafico
//legend1:		legenda do grafico para o data1
//legend2:		legenda do grafico para o data2
//filename:		nome do arquivo de saida (o grafico) sem extensao. Caso precisar, pode passar caminhos relativos ou completos
//			aceita caminhos para subdiretorios. Ex graficos/meugrafico
//
	FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
	if(!gnuplotPipe){
		printf("ERROR in gnuplotPipe!\n");
		printf("Something went wrong with popen()! %s\n", strerror(errno));
	       return;
	}

	char filename_data1[255];
	sprintf(filename_data1,"%s_%s.data",filename,legend1);
	FILE *f_data = fopen(filename_data1, "w");
	fputs(data1,f_data);
	fflush(f_data);
	fclose(f_data);

	free(data1);
	
	char filename_data2[255];
	sprintf(filename_data2,"%s_%s.data",filename,legend2);
	FILE *f_data2 = fopen(filename_data2, "w");
	fputs(data2,f_data2);
	fflush(f_data2);
	fclose(f_data2);

	free(data2);

	fputs(" set bmargin 7 \n",gnuplotPipe); 
	fputs(" unset colorbox \n",gnuplotPipe);
	fputs(" set terminal png enhanced font 'Verdana,10' \n",gnuplotPipe);
	fputs(" set output '",gnuplotPipe);
	fputs(filename,gnuplotPipe);
	fputs(".png' \n",gnuplotPipe );
	fputs(" set style line 2  lc rgb '#0404E9' lt 1 pt 7 \n",gnuplotPipe); //blue
	fputs(" set style line 3  lc rgb '#B22C2C' lt 1 pt 7 \n",gnuplotPipe); //blue


	char paramaxis[1024];
	sprintf(paramaxis," set xlabel '%s' \n set ylabel '%s' \n",xtitle,ytitle);
	fputs(paramaxis,gnuplotPipe);

	char paramtitle[1024];
	sprintf(paramtitle," set title '%s'\n ",title);
	fputs(paramtitle,gnuplotPipe);

	char paramplot1[1024];
	sprintf(paramplot1," plot '%s' using 1:2 title '%s' ls 2 with lines, ",filename_data1,legend1);
	char paramplot2[1024];
	sprintf(paramplot2," '%s' using 1:2 title '%s' ls 3 with lines\n ",filename_data2,legend2);
	fputs(paramplot1,gnuplotPipe);
	fputs(paramplot2,gnuplotPipe);

	fputs("\nquit",gnuplotPipe);
	fflush(gnuplotPipe);
	fclose(gnuplotPipe);

}/*}}}*/

