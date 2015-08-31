/* 
        Description: The ANSI C code of the PRVNS approach
        Programmer:  Wesklei Migliorini
        E-mail:      wesklei.m@gmail.com
        Date:	     04/11/2014
        Lisence:     Free
        Note:        The system was developed using Linux.
        To compile:  Type: make
        To run: ./algorithm input.in
 */
#ifndef functions_h
#define functions_h

#include <assert.h>
#include "data.h"
#define PI 3.14159265

//aux methods
void prepararObjFunc(int* FUNCTION, double* lb, double* ub);
double objfunc(double sol[], const int* FUNCTION, const int* DIM, int *cont);
char *getFunctionName(int FUNCTION);

//functions
double rastrigin( double SOL[],  int DIM);
double schaffer( double SOL[],  int DIM);
double griewank( double SOL[],  int DIM);
double ackley( double SOL[],  int DIM);
double rosenbrock( double SOL[],  int DIM);
double sphere( double SOL[],  int DIM);
double schaffer_f6( double sol[],  int DIM);
double g_schwefels( double sol[],  int DIM);
double mpe( double sol[],  int DIM);
double shubert(double SOL[],  int DIM);
double powell( double sol[], int DIM);
double michalewitz( double sol[], int DIM);
double levy( double sol[], int DIM);
double zakharov( double sol[], int DIM);
double egg_holder( double sol[], int DIM);
double generalized_holzman( double sol[], int DIM);
double rana( double sol[], int DIM);
double holzman( double sol[], int DIM);
double schwefels222( double sol[], int DIM);
double stretchedV( double sol[], int DIM);
double step(double sol[], int DIM);
double penalized1(double sol[], int DIM);
double penalized2(double sol[], int DIM);
double TempValue(double x,int a,int k,int m);
double multimod(double sol[], int DIM);

//shifted functions
double shifted_sphere( double sol[], int DIM);
double shifted_schwefel_221( double sol[], int DIM);
double shifted_rosenbrock( double sol[], int DIM);
double shifted_rastrigin( double sol[], int DIM);
double shifted_griewank( double sol[], int DIM);
double shifted_ackley( double sol[], int DIM);
double shifted_schwefel_222(double sol[], int dim); 
double shifted_schwefel_12(double sol[], int dim);
double f_10(double x, double y);// auxiliar
double shifted_extended_f10(double sol[], int dim);
double shifted_bohachevsky(double sol[], int dim);
double shifted_schaffer(double sol[], int dim);

//hybrid functions
static void divideFunctions(double sol[], int DIM, double *part1, double *part2, double m, int *psize1, int *psize2); 
double Extended_f_10NoDesplazamiento(double sol[], int DIM);
double f_BohachevskyNoDesplazamiento(double sol[], int DIM);
double f_Schwefel2_22NoDesplazamiento(double sol[], int DIM);
double hybrid_1(double sol[], int DIM); 
double hybrid_2(double sol[], int DIM); 
double hybrid_3(double sol[], int DIM); 
double hybrid_4(double sol[], int DIM); 
double hybrid_5(double sol[], int DIM); 
double hybrid_6(double sol[], int DIM); 
double hybrid_7(double sol[], int DIM); 
double hybrid_8(double sol[], int DIM); 

#endif
