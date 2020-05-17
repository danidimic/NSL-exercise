#include <armadillo>
#include "random.h"
using namespace std;
using namespace arma;

#define m 1
#define hbar 1

Random rnd;
int nblock = 100, nstep = 5E5, nbins = 150;
double delta, mu, sigma; 

//Minimo tra due valori
double min(double a, double b){
	double Min;
	if(a<b) Min = a;
	else	Min = b;
	return Min;
}

//Funzione d'onda tentativo
double Psi_trial(double x){
	double a = exp( -(pow(x-mu,2))/(2*sigma*sigma) );
	double b = exp( -(pow(x+mu,2))/(2*sigma*sigma) );
	return a + b;
}

//Densità di probabilità
double pdf(double x){
	double f = Psi_trial(x);
	return f*f;
}

//Hamiltoniana del sistema
double Hamiltonian(double x){
	double K, V, a, b;

	//Termine cinetico
	a = exp( -(pow(x-mu,2))/(2*sigma*sigma) );
	b = exp( -(pow(x+mu,2))/(2*sigma*sigma) );
	K = pow( (x-mu)/(sigma*sigma), 2 )*a + pow( (x+mu)/(sigma*sigma), 2 )*b - (a+b)/(sigma*sigma);
	K *= -pow(hbar,2)/(2*m);

	//Potenziale
	V = pow(x,4) - 5.0/2.0*pow(x,2);
	V *= Psi_trial(x);

	return K+V;
}


void Input();
vec gethisto(vec, int);
double MCintegral(vec x);


