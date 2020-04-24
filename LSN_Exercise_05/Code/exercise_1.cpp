#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>
#include "header.h"
#include "random.h"
#include "statistics.h"
using namespace std;
using namespace arma;


int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	vec  x = ones<vec>(nstep), xnext(nstep);
	mat T = symmatu( randu<mat>(nstep, nstep) );
	//mat U = symmatu( randn<mat>(nstep, nstep) );
	cout<<x*T<<endl;
	
	double sum;
	for(int i=0; i<T.n_rows; i++){
		sum = 0;
		for(int j=0; j<T.n_cols; j++)	sum +=T(i,j);

		for(int j=0; j<T.n_cols; j++)	T(i,j) /= sqrt(sum);
	}
	cout<<T<<endl;

	
	
		
	rnd.SaveSeed();
	return 0;
}

//densità di probabilità target	
double pdf(double x){
	double r;
	return r=gauss(5, 1, x);
}
//minimo tra due valori
double min(double a, double b){
	double r;
	if(a<b)	r=a;
	else	r=b;
	return r;
}
//gaussiana normalizzata
double gauss(double mu, double sigma, double x){
	double r = exp( -pow(mu-x,2) / (2*pow(sigma,2)) );
	return r/(sqrt(2*M_PI)*sigma);
}
