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

#include "functions.h"

void prepararObjFunc(int* FUNCTION, double* lb, double* ub)/*{{{*/
{
	switch (*FUNCTION) {
		case 0: //Rastrigin
			*lb = -5.12;
			*ub = 5.12;
			break;
		case 1: //Schaffer
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 2: //Griewank
			*lb = -600.00;
			*ub = 600.00;
			break;
		case 3: //Ackley
			*lb = -32.00;
			*ub = 32.00;
			break;
		case 4: //Rosenbrock
			*lb = -30.00;
			*ub = 30.00;
			break;
		case 5: //Sphere
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 6: //MPE
			*lb = 0;
			*ub = 5;
			break;
		case 7: //SCHAFFER_F6
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 8: //Generalized Schwefel's function 2.26 
			*lb = -500.00;
			*ub = 500.00;
			break;
		case 9: // Step function
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 11: // Generalized Penalized function #1
			*lb = -500.00;
			*ub = 500.00;
			break;
		case 12: // Levy Function
			*lb = -10.00;
			*ub = 10.00;
			break;
		case 13: // Zakharov
			*lb = -5.00;
			*ub = 10.00;
			break;
		case 14: // EggHolder
			*lb = -512.00;
			*ub = 512.00;
			break;
		case 15: // Holzman
			*lb = -10.00;
			*ub = 10.00;
			break;
		case 16: // Michalewitz
			*lb = 0.00;
			*ub = PI;
			break;
		case 17: // Generalized penalized function #2
			*lb = -50.00;
			*ub = 50.00;
			break;
		case 18: // Powell
			*lb = -4.00;
			*ub = 5.00;
			break;
		case 19: // Rana
			*lb = -512.00;
			*ub = 512.00;
			break;
		case 20: // Shubert
			*lb = -10.00;
			*ub = 10.00;
			break;
		case 21: //StretchedV
			*lb = -10.00;
			*ub = 10.00;
			break;
		case 22: // Multimod 
			*lb = -10.00;
			*ub = 10.00;
			break;
		case 23: // Schwefels 222
			*lb = -10.00;
			*ub = 10.00;
			break;

			//Shifted functions

		case 25: //Shifted Sphere
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 26: //Shifted Schwefel's Problem 2.21
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 27: //Shifted Rosenbrock
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 28: //Shifted Rastrigin
			*lb = -5.12;
			*ub = 5.12;
			break;
		case 29: //Shifted Griewank
			*lb = -600.00;
			*ub = 600.00;
			break;
		case 30: //Shifted Ackley
			*lb = -32.00;
			*ub = 32.00;
			break;
		case 31: //Shifted Schwefel 2.22
			*lb = -10.00;
			*ub = 10.00;
			break;
		case 32: //Shifted Schwefel 1.2
			*lb = -65.536;
			*ub = 65.536;
			break;
		case 33: //Shifted Extended_f10
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 34: //Shifted Bohachevsky
			*lb = -15.00;
			*ub = 15.00;
			break;
		case 35: //Shifted schaffer
			*lb = -100.00;
			*ub = 100.00;
			break;

			//Hybrid functions

		case 36: //Hybrid 1
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 37: //Hybrid 2
			*lb = -5.00;
			*ub = 5.00;
			break;
		case 38: //Hybrid 3
			*lb = -10.00;
			*ub = 10.00;
			break;
		case 39: //Hybrid 4
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 40: //Hybrid 5
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 41: //Hybrid 6
			*lb = -100.00;
			*ub = 100.00;
			break;
		case 42: //Hybrid 7
			*lb = -5.00;
			*ub = 5.00;
			break;
		case 43: //Hybrid 8
			*lb = -10.00;
			*ub = 10.00;
			break;

		default:
			printf("Info: Invalid function.\n") ;
			exit(0);
	}
}/*}}}*/

double objfunc(double sol[],const int* FUNCTION, const int* DIM, int *cont)/*{{{*/
{
	*cont+=1;

	switch (*FUNCTION) {

		case 0: //Rastrigin
			return rastrigin(sol,*DIM);
		case 1: //Schaffer
			return schaffer(sol,*DIM);
		case 2: //Griewank
			return griewank(sol,*DIM);
		case 3: //Ackley
			return ackley(sol,*DIM);
		case 4: //Rosenbrock
			return rosenbrock(sol,*DIM);
		case 5: //Sphere
			return sphere(sol,*DIM);
		case 6: //MPE
			return mpe(sol,*DIM);
		case 7: //SCHAFFER_F6
			return schaffer_f6(sol,*DIM);
			break;
		case 8: //Generalized Schwefel's function 2.26
			return g_schwefels(sol,*DIM);
			break;
		case 9: // Step Function
			return step(sol,*DIM);
			break;
		case 11: // Generalized Penalized function #1
			return penalized1(sol,*DIM);
			break;
		case 12: // Levy Function
			return levy(sol,*DIM);
			break;
		case 13: // Zakharov
			return zakharov(sol,*DIM);
			break;
		case 14: // EggHolder
			return egg_holder(sol,*DIM);
			break;
		case 15: // Holzman
			return holzman(sol,*DIM);
			break;
		case 16: // Michalewitz
			return michalewitz(sol,*DIM);
			break;
		case 17: // Generalized Penalized function #2
			return penalized2(sol,*DIM);
			break;
		case 18: // Powell
			return powell(sol,*DIM);
			break;
		case 19: // Rana
			return rana(sol,*DIM);
			break;
		case 20: // Shubert
			return shubert(sol,*DIM);
			break;
		case 21: // StretchedV
			return stretchedV(sol,*DIM);
			break;
		case 22: // Multimod
			return multimod(sol,*DIM);
			break;
		case 23: //  Schwefel's function 2.22
			return schwefels222(sol,*DIM);
			break;

			//Shifted functions

		case 25: //Shifted_Sphere
			return shifted_sphere(sol,*DIM);
			break;
		case 26: //Shifted Schwefel Problem 2.21
			return shifted_schwefel_221(sol,*DIM);
			break;
		case 27: //Shifted_Rosenbrock
			return shifted_rosenbrock(sol,*DIM);
			break;
		case 28: //Shifted Rastrigin
			return shifted_rastrigin(sol,*DIM);
			break;
		case 29: //Shifted Griewank
			return shifted_griewank(sol,*DIM);
			break;
		case 30://Shifted Ackley
			return shifted_ackley(sol,*DIM);
			break;
		case 31: //Shifted Schwefel 2.22
			return shifted_schwefel_222(sol,*DIM);
			break;
		case 32: //Shifted Schwefel 1.2
			return shifted_schwefel_12(sol,*DIM);
			break;
		case 33: //Shifted Extended_f10
			return shifted_extended_f10(sol,*DIM);
			break;
		case 34: //Shifted Bohachevsky
			return shifted_bohachevsky(sol,*DIM);
			break;
		case 35: //Shifted schaffer
			return shifted_schaffer(sol,*DIM);
			break;

			//Hybrid functions
			
		case 36: //Hybrid 1
			return hybrid_1(sol,*DIM);
			break;
		case 37: //Hybrid 2
			return hybrid_2(sol,*DIM);
			break;
		case 38: //Hybrid 3
			return hybrid_3(sol,*DIM);
			break;
		case 39: //Hybrid 4
			return hybrid_4(sol,*DIM);
			break;
		case 40: //Hybrid 5
			return hybrid_5(sol,*DIM);
			break;
		case 41: //Hybrid 6
			return hybrid_6(sol,*DIM);
			break;
		case 42: //Hybrid 7
			return hybrid_7(sol,*DIM);
			break;
		case 43: //Hybrid 8
			return hybrid_8(sol,*DIM);
			break;
		default:
			printf("Info: Invalid function.\n") ;
			exit(0);
	}
}/*}}}*/

char *getFunctionName(int FUNCTION){/*{{{*/

	switch (FUNCTION) {
		case 0: //Rastrigin
			return "Rastrigin";
		case 1: //Schaffer
			return "Schaffer F7";
		case 2: //Griewank
			return "Griewank";
		case 3: //Ackley
			return "Ackley";
		case 4: //Rosenbrock
			return "Rosenbrock";
		case 5: //Sphere
			return "Sphere";
		case 6: //MPE
			return "Molecular Potential Energy";
		case 7: //SCHAFFER_F6
			return "Schaffer F6";
			break;
		case 8: //G_SCHWEFELS
			return "Generalized Schwefels 2.26";
			break;
		case 9: // Step function
			return "Step";
			break;
		case 11: // Generalized Penalized function #1
			return "Generalized Penalized function #1";
			break;
		case 12: // Levy
			return "Levy";
			break;
		case 13: // Zakharov
			return "Zakharov";
			break;
		case 14: // Egg holder
			return "Egg holder";
			break;
		case 15: // Generalized Holzman
			return "Generalized Holzman";
			break;
		case 16: // Michalewitz
			return "Michalewitz";
			break;
		case 17: // Generalized penalized function #2
			return "Generalized penalized function #2";
			break;
		case 18: // Powell
			return "Powell";
			break;
		case 19: // Rana
			return "Rana";
			break;
		case 20: // Shubert
			return "Shubert";
			break;
		case 21: //StretchedV
			return "StretchedV";
			break;
		case 22: //Multimod
			return "Multimod";
			break;
		case 23: //  Schwefel's function 2.22
			return "Schwefel's function 2.22";
			break;

			//Shifted functions

		case 25: //Shifted Sphere
			return "Shifted Sphere";
			break;
		case 26: //Shifted Schwefel's Problem 2.21
			return "Shifted Schwefel's Problem 2.21";
			break;
		case 27: //Shifted Rosenbrock
			return "Shifted Rosenbrock";
			break;
		case 28: //Shifted Rastrigin
			return "Shifted Rastrigin";
			break;
		case 29: //Shifted Griewank
			return "Shifted Griewank";
			break;
		case 30: //Shifted Ackley
			return "Shifted Ackley";
			break;
		case 31: //Shifted Schwefel 2.22
			return "Shifted Schwefel 2.22";
			break;
		case 32: //Shifted Schwefel 1.2
			return "Shifted Schwefel 1.2";
			break;
		case 33: //Shifted Extended_f10
			return "Shifted Extended_f10";
			break;
		case 34: //Shifted Bohachevsky
			return "Shifted Bohachevsky";
			break;
		case 35: //Shifted schaffer
			return "Shifted Schaffer";
			break;

			//Hybrid functions

		case 36: //Hybrid 1
			return "Hybrid 1";
			break;
		case 37: //Hybrid 2
			return "Hybrid 2";
			break;
		case 38: //Hybrid 3
			return "Hybrid 3";
			break;
		case 39: //Hybrid 4
			return "Hybrid 4";
			break;
		case 40: //Hybrid 5
			return "Hybrid 5";
			break;
		case 41: //Hybrid 6
			return "Hybrid 6";
			break;
		case 42: //Hybrid 7
			return "Hybrid 7";
			break;
		case 43: //Hybrid 8
			return "Hybrid 8";
			break;

		default:
			printf("Info: Invalid function.\n") ;
			exit(0);
	}
}/*}}}*/

double rastrigin( double sol[],  int DIM){/*{{{*/
	//sol[] -> the population set
	//DIM the dimension of sol[]
	int j;
	double top = 0.00;

	for(j=0;j<DIM;j++)
	{
		top=top+(pow(sol[j],(double)2)-10*cos(2*M_PI*sol[j])+10);
	}
	return top;
}/*}}}*/

double schaffer( double sol[],  int DIM){/*{{{*/
	//sol[] -> the population set
	//DIM the dimension of sol[]

	int j;
	double top = 0.00 , top1 = 0.00;

	for(j=0;j<DIM;j++)
	{
		top=top+(pow(sol[j],(double)2));
	}

	top = pow(top,(double)0.25);

	for(j=0;j<DIM;j++)
	{
		top1=top1+(pow(sol[j],(double)2));
	}

	top1=pow(top1,(double)0.1);
	top1 = pow(sin(50*top1),(double)2) +1.0;

	return top*top1;
}/*}}}*/

double griewank( double sol[],  int DIM){/*{{{*/
	//sol[] -> the population set
	//DIM the dimension of sol[]

	int j;
	double top = 0.00 , top1 = 0.00, top2 = 0.00;

	for(j=0;j<DIM;j++)
	{
		top1=top1+pow((sol[j]),(double)2);
		top2=top2*cos((((sol[j])/sqrt((double)(j+1)))*M_PI)/180);
	}
	top=(1/(double)4000)*top1-top2+1;

	return top;
}/*}}}*/

double ackley( double sol[],  int DIM){/*{{{*/
	//sol[] -> the population set
	//DIM the dimension of sol[]

	int i;	
	double aux = 0.0;
	double aux1 = 0.0;

	for (i = 0; i < DIM; i++)
	{
		aux += sol[i]*sol[i];
	}
	for (i = 0; i < DIM; i++)
	{
		aux1 += cos(2.0*M_PI*sol[i]);
	}

	return (-20.0*(exp(-0.2*sqrt(1.0/(float)DIM*aux)))-exp(1.0/(float)DIM*aux1)+20.0+exp(1));
}/*}}}*/

double rosenbrock( double sol[],  int DIM){/*{{{*/
	//sol[] -> the population set
	//DIM the dimension of sol[]

	int i;
	double top = 0.00;

	for (i = 0; i < DIM-1; i++)
	{
		top=top+100.*pow((sol[i+1] - pow(sol[i],2.)),2) + pow((1. - sol[i]),2);
	}

	return top;
}/*}}}*/

double sphere( double sol[],  int DIM){/*{{{*/
	//sol[] -> the population set
	//DIM the dimension of sol[]

	int j;
	double top = 0.00;

	for(j=0;j<DIM;j++)
	{
		top=top+sol[j]*sol[j];
	}

	return top;
}/*}}}*/

double shubert(double sol[],  int DIM){/*{{{*/
	//sol[] -> the population set
	//DIM the dimension of sol[]
	//Shubert
	/*
	   -   Domain  |x| <= 10.0
	   -   Number of local minimum = 400
	   -   Global minimum fmin = -24.062499 at the ff. points
	   -    (-6.774576, -6.774576), ..., (5.791794, 5.791794)
	   */
	double sum = 0.0;
	int i;
	for (i = 0; i < DIM; i++) {
		sum += -sin(2.0*sol[i]+1.0)
			-2.0*sin(3.0*sol[i]+2.0)
			-3.0*sin(4.0*sol[i]+3.0)
			-4.0*sin(5.0*sol[i]+4.0)
			-5.0*sin(6.0*sol[i]+5.0);
	}

	return sum/(DIM/2.0);
}/*}}}*/

double schaffer_f6( double sol[],  int DIM){/*{{{*/
	double top=0;
	double top1=0;
	double top2=0;
	int j;
	for(j=0;j<DIM-1;j++)
	{
		top1 = pow( sin( sqrt( pow(sol[j],(double)2) + pow(sol[j+1],(double)2) ) ), (double)2) - 0.5;
		top2 = pow( 1.0 + 0.001* ( pow(sol[j],(double)2) + pow(sol[j+1],(double)2) ), (double)2);
		top = top + (0.5 + top1/top2);
	}
	return top;
}/*}}}*/

double g_schwefels( double sol[],  int DIM){/*{{{*/
	//known_optimal = -418.982887272433 at sol(i)=420.9687
	int i;
	double aux = 0.0;
	for (i=0;i<DIM;i++)
	{
		aux += (-1 * sol[i]) * sin(sqrt(fabs(sol[i])));
	}
	return(aux);
}/*}}}*/

double mpe( double sol[], int DIM){/*{{{*/
	//[0,5] fmin=-0.411183034 * n
	int i;
	double aux1=0.0,aux2=0.0;
	for (i=0;i<DIM;i++)
	{
		aux1 += 1 + cos(3*sol[i]);
		aux2 += (pow(-1,i) / (sqrt(fabs(10.60099896-4.141720682 * cos(sol[i])))));
	}
	return aux1 + aux2;
}/*}}}*/

double michalewitz( double sol[], int DIM){/*{{{*/
	int i;
	double aux;
	//Michalewitz
	/*
Dimension: n
Domain: 0< | x[i] | <= PI
Global minimum: x[] = -0.966*n
*/
	aux=0;
	for (i=0;i<DIM;i++) {
		aux = aux + sin(sol[i]) * pow(sin((i+1)*sol[i]*sol[i]/(float)PI), 2.0 * 10.0);
	}
	return(-1*aux);
}/*}}}*/

double powell( double sol[], int DIM){/*{{{*/
	int j;
	double aux; //Powell
	/*
Dimension: n > 4
Domain: -4<= x[i] <= 5
Global minimum: at (3, -1, 0, 1, ..., 3, -1, 0, 1) with fmin = 0
*/
	aux = 0.0;
	for (j = 1; j <= (int)DIM/4; j++) {
		aux +=  pow(sol[4*j-4] + 10 * sol[4*j-3],2.0)
			+ 5 * pow(sol[4*j-2] - sol[4*j-1],2.0)
			+ pow(sol[4*j-3] - 2 * sol[4*j-2], 4.0)
			+ 10 * pow(sol[4*j - 4] - sol[4*j-1], 4.0);
	}
	return aux;
}/*}}}*/

double levy( double sol[], int DIM)/*{{{*/
{
	//Levy Function
	double aux,*y;
	int i;
	//x[i] = 1 f(x[i])=0
	y = (double*) malloc (DIM * sizeof(double));
	for (i = 0; i< DIM; i++)
		y[i] = 1+(sol[i]-1)/4.0;
	aux = pow(sin(PI*y[0]),2.0);
	for (i = 0; i<DIM-1;i++)
		aux = aux + pow(y[i]-1,2.0)*(1+10*pow(sin(PI*y[i]+1),2.0));

	aux = aux+pow(y[DIM-1]-1,2.0)*(1+pow(sin(2*PI*y[DIM-1]),2.0) );

	free (y);
	return ( aux );
}/*}}}*/

double zakharov( double sol[], int DIM)/*{{{*/
{
	// Zakharov function  //x[i] = 0 f(x[i])=0
	//
	double aux,aux1;
	int j;
	aux = aux1 = 0.0;
	for (j = 0; j< DIM; j++)
	{
		aux = aux + pow(sol[j],2.0);
		aux1 = aux1+0.5*j*sol[j];
	}


	return ( aux+pow(aux1,2.0)+pow(aux1,4.0) );
}/*}}}*/

double egg_holder( double sol[], int DIM)/*{{{*/
{
	//Egg holder
	/*
	   - Dimension: n
	   - Domain:  -512 < | x_i | < 512
	   - Minimum for n=2 fmin(512, 404.2319) = -959.641
	   */
	double aux;
	int i;
	aux = 0.0;
	for (i = 0; i < DIM-1; i++) 
	{
		aux += -(sol[i+1] + 47.0) * sin(sqrt(fabs(sol[i+1] + sol[i] * 0.5 + 47.0))) + sin(sqrt(fabs(sol[i] - (sol[i+1] + 47.0)))) * (-sol[i]);
	}
	return (aux);
}/*}}}*/

double rana( double sol[], int DIM)/*{{{*/
{
	//Rana
	/*
Dimension: n
Domain: -520<= x[i] <= 520
Global minimum: ???
*/
	double sum,t1,t2;
	int i;
	sum = 0.0;
	for (i = 0; i < DIM-1; i++) {
		/* if(sol[i] < -520.00f || sol[i] > 520.00f) */
		/* 	printf("rana %f\n",sol[i]); */

		t1 = sqrt(fabsf(sol[i+1] + sol[i] + 1.0));
		t2 = sqrt(fabsf(sol[i+1] - sol[i] + 1.0));
		sum += (sol[i+1] + 1.0) * cos(t2) * sin(t1) + cos(t1) * sin(t2) * sol[i];
	}
	return sum/(double)(DIM-1);
}/*}}}*/

double holzman( double sol[], int DIM){/*{{{*/
	//Generalized Holzman
	/*
Dimension: n
Domain: | x[i] | <= 10
Global minimum: 0 at x[i] = 0
*/
	int i;
	double aux = 0.0;
	for (i = 0; i < DIM; i++) 
	{
		aux += i * pow(sol[i] , 4);
	}
	return aux;
}/*}}}*/

double schwefels222( double sol[], int DIM){/*{{{*/
	// Schwefel's function 2.22
	/*
	   - Domain:  | x_i | <= 10.0
	   Global minimum is 0.0 at x_i = 0.00
	   */

	int i;
	double aux = 0.0,aux1=0.0;
	for (i=0;i<DIM;i++)
	{
		aux += fabs(sol[i]);
		aux1 *= fabs(sol[i]);
	}

	return (aux+aux1);
}/*}}}*/

double stretchedV( double sol[], int DIM){/*{{{*/
	//StretchedV
	/*
	   - Domain:  | x_i | <= 10.0
	   Global minimum is 0.0 at x_i = 0.00
	   */
	double sum = 0.0;
	double aux;
	int i;
	for (i = 0; i < DIM-1; i++) {
		aux = sol[i+1]*sol[i+1] + sol[i]*sol[i];
		sum += pow(aux, 0.25) * (pow(sin(50.0 * pow(aux, 0.1)), 2.0)+1.0);
	}
	return sum;
}/*}}}*/

double step(double sol[], int DIM){/*{{{*/
	// Step function
	/*
	   -   Domain: | x_i  | < 100.0
	   -   Global minimum is 0.0 at x_i = 0.5
	   */
	double aux,aux1=0.0;
	int i;
	for (i=0;i<DIM;i++)
	{
		aux = (sol[i]+0.5);
		aux1 += aux*aux; 
	}
	return (aux1);
}/*}}}*/

double tempValue(double x,int a,int k,int m)/*{{{*/
{
	double temp = 0.0;
	if( x > a)
	{
		temp = k*pow(x-a,m);
	}
	else if( x <= a && x >= -a)
	{
		temp = 0.0;
	}
	else
	{
		temp = k*pow(-x-a,m);
	}
	return temp;
}/*}}}*/

double penalized1(double sol[], int DIM){/*{{{*/
	//Generalized Penalized function #1
	// -500 <= xi <= 500
	//known_optimal = 1.57044103551786e-032;
	
	double aux=0.0;
	double aux1=0.0;
	int i;
	double *y = (double*) malloc (DIM * sizeof(double));
	for (i=0;i<DIM;i++)
	{
		y[i]=0.0;
	}

	for (i=0;i<DIM;i++)
	{
		y[i]=1+(sol[i]+1)/4.0;
	}

	for (i=0;i<DIM-1;i++)
	{
		aux += pow(y[i]-1,2.0)*(1.0+10.0*pow(sin(PI*y[i+1]),2.0)); 
	}
	for (i=0;i<DIM;i++)
	{
		aux1 += tempValue(sol[i],10,100,4);
	}

	aux = (10.0*pow(sin(PI*y[0]),2.0)+aux+pow(y[DIM-1]-1,2))*PI/DIM+aux1;

	free(y);

	return ( aux );
}/*}}}*/

double penalized2(double sol[], int DIM){/*{{{*/
	//Generalized Penalized function #2
	//known_optimal = 1.34969464963992e-032;
	double aux= 0.0;
	double aux1=0.0;
	int i;
	for (i=0;i<DIM-1;i++)
	{
		aux += pow(sol[i]-1,2.0)*(1.0+10.0*pow(sin(3*PI*sol[i+1]),2.0)); 
	}
	for (i=0;i<DIM;i++)
	{
		aux1 += tempValue(sol[i],5,100,4);
	}

	return ( (pow(sin(3.0*PI*sol[0]),2.0)+aux+pow(sol[DIM-1]-1,2.0)
				*(1.0+pow(sin(2.0*PI*sol[DIM-1]),2.0)))/10.0+aux1 );
}/*}}}*/

double multimod(double sol[], int DIM){/*{{{*/
	//Multimod
	/*
	Dimension: n
	Domain: -10<= x[i] <= 10
	Global minimum: x[0] = 0
	*/
	double s,p,t;
	int i;
	s = p = fabs(sol[0]);
	for (i = 1; i < DIM; i++) {
		t = fabs(sol[i]);
		s += t;
		p *= t;
	}
	return s + p;
}/*}}}*/

//=== Shifted Functions

double shifted_sphere( double sol[], int DIM){/*{{{*/
	// x* = o , F(x*) = f_bias1 = - 450
	int i;
	double Fx = 0.0;
	double z = 0.0;
	for(i=0;i<DIM;i++){
		z = sol[i] - sphere_data[i];
		Fx += z*z;
	}
	return Fx + f_bias[0];
}/*}}}*/

double shifted_schwefel_221( double sol[], int DIM){/*{{{*/
	//Shifted Schwefel Problem 2.21
	// x* = o , F(x*) = f_bias1 = - 450
	int i;
	double Fx = 0.0;
	double z = 0.0;
	Fx = abss(sol[0]);
	for(i=1;i<DIM;i++){
		z = sol[i] - schwefel_data[i];
		Fx = max(Fx , abss(z));
	}
	return Fx + f_bias[1]; 
}/*}}}*/

double shifted_rosenbrock( double sol[], int DIM){/*{{{*/
	int i;
	double Fx = 0.0;
	double zx[DIM];
	//Shifted_Rosenbrock
	// x* = o , F(x* ) = f_bias3 = 390

	for(i=0;i<DIM;i++) zx[i] = sol[i] - rosenbrock_data[i] + 1;   

	for(i=0;i<DIM-1;i++){    
		Fx = Fx + 100*( pow((pow(zx[i],2)-zx[i+1]) , 2) ) + pow((zx[i]-1) , 2);
	}
	return Fx + f_bias[2]; 
}/*}}}*/

double shifted_rastrigin( double sol[], int DIM){/*{{{*/
	int i;
	double Fx = 0.0;
	double z = 0.0;
	double pi = acos(-1.0);
	//Shifted Rastrigin
	//x* = o , F( x * ) = f_bias4 = - 330
	for(i=0;i<DIM;i++){  
		z = sol[i] - rastrigin_data[i];
		Fx = Fx + ( pow(z,2) - 10*cos(2*pi*z) + 10);
	}
	return Fx + f_bias[3];
}/*}}}*/

double shifted_griewank( double sol[], int DIM){/*{{{*/
	int i;
	double z = 0.0;
	double top1 = 0.00, top2 = 0.00;

	//Shifted Griewank
	//x* = o , F(x* ) = f_bias5 = -180
	top1 = 0;
	top2 = 1;
	for(i=0;i<DIM;i++){       
		z = sol[i] - griewank_data[i];
		top1 = top1 + ( pow(z,2) / 4000 );
		top2 = top2 * ( cos(z/sqrt(i+1)));

	}
	return (top1 - top2 + 1 + f_bias[4]);
}/*}}}*/

double shifted_ackley( double sol[], int DIM){/*{{{*/
	int i;
	double Fx = 0.0;
	double z = 0.0;
	double top1 = 0.00, top2 = 0.00;
	double pi = acos(-1.0);
	double e = exp(1.0);

	//Shifted Ackley
	// x* = o , F(x* ) = f_bias6 = - 140
	top1 = 0;
	top2 = 0;
	Fx = 0;
	for(i=0;i<DIM;i++){   
		z = sol[i] - ackley_data[i];
		top1 = top1 + pow(z , 2 );
		top2 = top2 + cos(2*pi*z);
	}
	Fx = -20*exp(-0.2*sqrt(top1/DIM)) -exp(top2/DIM) + 20 + e + f_bias[5];
	return Fx;
}/*}}}*/

double shifted_schwefel_222(double sol[], int DIM) /*{{{*/
{
	double sum, currentGen, prod;
	int i;

	sum = 0.0;
	prod = 1.0;

	for (i = 0; i < DIM; i++) 
	{
		currentGen = fabs(sol[i]-data_shif_schwefels_222[i]);
		sum += currentGen;
		prod *= currentGen;
	}

	return sum + prod;
}/*}}}*/

double shifted_schwefel_12(double sol[], int DIM)/*{{{*/
{
	double Sum=0.0, Val=0.0;
	int i;

	for (i = 0; i < DIM; i++)
	{  
		Val += sol[i]-data_shif_schwefels_12[i];
		Sum += Val * Val;
	}

	return Sum;
}/*}}}*/

double f_10(double x, double y)/*{{{*/
{
	double p, z, t;

	p=(x*x+y*y);

	z=pow(p, 0.25);
	t=sin(50.0*pow(p, 0.1));
	t=t*t+1.0;

	return z*t;
}/*}}}*/

double shifted_extended_f10(double sol[], int DIM)/*{{{*/
{
	double suma=0.0;
	int i;

	for(i=0; i<DIM-1; i++)
		suma+=f_10(sol[i]-data_shif_extendedf10[i], sol[i+1]-data_shif_extendedf10[i+1]);

	suma+=f_10(sol[DIM-1]-data_shif_extendedf10[DIM-1], sol[0]-data_shif_extendedf10[0]);

	return suma;
}/*}}}*/

double shifted_bohachevsky(double sol[], int DIM) /*{{{*/
{   
	double sum = 0.0;
	int i;
	double currentGen;
	double nextGen;

	currentGen = sol[0]-data_shif_bohachevsky[0];

	for (i = 1; i < DIM; i++) 
	{
		nextGen = sol[i]-data_shif_bohachevsky[i];
		sum += currentGen * currentGen + 2.0 * nextGen * nextGen;
		sum += -0.3 * cos(3.0 * PI * currentGen) -0.4 * cos(4.0 * PI * nextGen) + 0.7;
		currentGen = nextGen;
	}

	return sum;
}/*}}}*/

double shifted_schaffer(double sol[], int DIM) /*{{{*/
{
	int i=0;
	double sum;
	double aux, aux2;
	double currentGen, nextGen;

	sum = 0.0;
	currentGen = sol[0]-data_shif_schaffer[i];
	currentGen = currentGen * currentGen;

	for (i = 1; i < DIM; i++) 
	{
		nextGen = sol[i]-data_shif_schaffer[i];
		nextGen = nextGen * nextGen;
		aux = currentGen + nextGen;
		currentGen = nextGen;
		aux2 = sin(50. * pow(aux, 0.1));
		sum += pow(aux, 0.25) * (aux2 * aux2 + 1.0);
	}

	return sum;
}/*}}}*/

//=== Hybrid Functions

static void divideFunctions(double sol[], int DIM, double *part1, double *part2, double m, int *psize1, int *psize2) {/*{{{*/
	int shared;
	int rest, i, total;
	double *partrest;

	if (m <= 0.5) {
		partrest = part2;
	}
	else {
		partrest = part1;
		m = 1-m;
	}

	shared = (int) floor(DIM*m);
	rest = 2*shared;

	for (i = 0; i < shared; i++) {
		part1[i] = sol[i*2];
		part2[i] = sol[i*2+1];
	}
	total = DIM-shared;

	for (i = 0; i < total-shared; i++) {
		partrest[i+shared] = sol[i+rest];
	}

	*psize1 = shared;
	*psize2 = DIM-shared;

	if (partrest == part1) {
		int temp = *psize1;
		*psize1 = *psize2;
		*psize2 = temp;
	}
}/*}}}*/

double Extended_f_10NoDesplazamiento(double sol[], int DIM)/*{{{*/
{
	double suma=0.0;

	int i;
	for(i=0; i<DIM-1; i++)
		suma+=f_10(sol[i], sol[i+1]);

	suma+=f_10(sol[DIM-1], sol[0]);

	return suma;
}/*}}}*/

double f_BohachevskyNoDesplazamiento(double sol[], int DIM) /*{{{*/
{   
    double sum = 0.0;
    int i;
    double currentGen;
    double nextGen;

    currentGen = sol[0];

    for (i = 1; i < DIM; i++) 
    {
        nextGen = sol[i];
        sum += currentGen * currentGen + 2.0 * nextGen * nextGen;
	double c1=cos(3.0 * PI * currentGen);
	double c2=cos(4.0 * PI * nextGen);
        sum += 0.7 -(0.3*c1+0.4*c2);
        currentGen = nextGen;
    }

    return sum;
}/*}}}*/

double f_Schwefel2_22NoDesplazamiento(double sol[], int DIM) /*{{{*/
{
    double sum, currentGen, prod;

    sum = 0.0;
    prod = 1.0;

    int i;
    for (i = 0; i < DIM; i++) 
    {
        currentGen = fabs(sol[i]);
        sum += currentGen;
        prod *= currentGen;
    }

    return sum + prod;
}/*}}}*/

double hybrid_1(double sol[], int DIM) /*{{{*/ //12
{
	double part1[DIM], part2[DIM];
	int size1, size2;
	double f1, f2;
	divideFunctions(sol, DIM, part1, part2, 0.25, &size1, &size2);

	f1=Extended_f_10NoDesplazamiento(part1,size1);
	f2=shifted_sphere(part2,size2)-f_bias[0];
	assert(f1 >= 0);
	assert(f2 >= 0);

	return f1+f2;
}/*}}}*/

double hybrid_2(double sol[], int DIM) /*{{{*/ //13
{
	double part1[DIM], part2[DIM];
	int size1, size2;
	double f1, f2;

	divideFunctions(sol, DIM, part1, part2, 0.25, &size1, &size2);
	f1=Extended_f_10NoDesplazamiento(part1,size1);
	f2=shifted_rosenbrock(part2,size2)-f_bias[2];

	return f1+f2;
}/*}}}*/

double hybrid_3(double sol[], int DIM) /*{{{*/ //14
{
	double part1[DIM], part2[DIM];
	int size1, size2;
	double f1, f2;
	divideFunctions(sol, DIM, part1, part2, 0.25, &size1, &size2);

	f1=Extended_f_10NoDesplazamiento(part1,size1);
	f2=shifted_rastrigin(part2,size2)-f_bias[3];

	return f1+f2; 
}/*}}}*/

double hybrid_4(double sol[], int DIM) /*{{{*/ //15
{
	double part1[DIM], part2[DIM];
	double desp[DIM];
	int size1, size2;
	double f1, f2;
	int i;

	for (i = 0; i < DIM; ++i) {
		desp[i] = sol[i] - f15[i];
	}

	divideFunctions(desp, DIM, part1, part2, 0.25, &size1, &size2);

	f1=f_BohachevskyNoDesplazamiento(part1, size1);
	f2=f_Schwefel2_22NoDesplazamiento(part2, size2);
	return f1+f2; 
}/*}}}*/

double hybrid_5(double sol[], int DIM) /*{{{*/ //16new
{
	double part1[DIM], part2[DIM];
	int size1, size2;
	double f1, f2;
	divideFunctions(sol, DIM, part1, part2, 0.5, &size1, &size2);

	f1=Extended_f_10NoDesplazamiento(part1,size1);
	assert(f1 >= 0);
	f2=shifted_sphere(part2, size2)-f_bias[0];
	assert(f2 >= 0);

	return f1+f2;
}/*}}}*/

double hybrid_6(double sol[], int DIM) /*{{{*/ //17new
{
	double part1[DIM], part2[DIM];
	int size1, size2;
	double f1, f2;

	divideFunctions(sol, DIM, part1, part2, 0.75, &size1, &size2);
	f1=Extended_f_10NoDesplazamiento(part1, size1);
	f2=shifted_rosenbrock(part2, size2)-f_bias[2];

	return f1+f2;
}/*}}}*/

double hybrid_7(double sol[], int DIM) /*{{{*/ //18new
{
	double part1[DIM], part2[DIM];
	int size1, size2;
	double f1=0;
	double f2=0;

	divideFunctions(sol, DIM, part1, part2, 0.75, &size1, &size2);

	f1=Extended_f_10NoDesplazamiento(part1, size1);
	f2=shifted_rastrigin(part2, size2)-f_bias[3];
	assert(isfinite(f1));
	assert(isfinite(f2));

	return f1+f2;
}/*}}}*/

double hybrid_8(double sol[], int DIM) /*{{{*/ //19new
{
	double part1[DIM], part2[DIM];
	int size1, size2;
	double desp[DIM];
	double f1, f2;

	int i;
	for (i = 0; i < DIM; ++i) {
		desp[i] = sol[i] - f19[i];
	}

	divideFunctions(desp, DIM, part1, part2, 0.75, &size1, &size2);

	f1=f_BohachevskyNoDesplazamiento(part1, size1);
	f2=f_Schwefel2_22NoDesplazamiento(part2, size2);

	return f1+f2;
}/*}}}*/

