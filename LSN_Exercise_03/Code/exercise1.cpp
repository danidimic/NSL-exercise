#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "statistics.h"
using namespace std;

#define S0 100	//prezzo iniziale
#define T 1		//tempo finale
#define K 100	//prezzo di mercato
#define r 0.1		//interesse privo di rischio
#define sigma 0.25	//volatilità

//massimo tra due valori
double max(double a, double b){
	double Max;
	if (a>b) Max=a;
	else Max=b;
	return Max;
}

int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	int nstep = 1E4, nblock = 100, nint = 100;

	double t, t1, x, z, S;
	double sumC_direct, sumP_direct;
	double sumC_discr, sumP_discr;
	vector<double> C_direct(nblock), P_direct(nblock);
	vector<double> C_discr(nblock),  P_discr(nblock);

	for(int i=0; i<nblock; i++){
		sumC_direct = 0;
		sumP_direct = 0;
		sumC_discr = 0;
		sumP_discr = 0;

		for(int j=0; j<nstep; j++){
			//simulazione diretta del prezzo
			z = rnd.Gauss(0, T);
			S = S0*exp( (r-pow(sigma,2)/2)*T + sigma*z );

			x = exp(-r*T)*max(0, S-K);	//profitto Call-option
			sumC_direct += x;
			x = exp(-r*T)*max(0, K-S);	//profitto Put-option
			sumP_direct += x;

			//simulazione discretizzata del prezzo
			S = S0;
			t = 0;
			for(int k=0; k<nint; k++){
				t1 = t + 1./nint;
				z = rnd.Gauss(0,1);
				S = S*exp( (r-pow(sigma,2)/2)*(t1-t) + sigma*z*sqrt(t1-t) );
				t = t1;
			}
			x = exp(-r*T)*max(0, S-K);	//profitto Call-option
			sumC_discr += x;

			x = exp(-r*T)*max(0, K-S);	//profitto Put-option
			sumP_discr += x;
		}
		//calcolo diretto
		C_direct[i] = sumC_direct/nstep;	//media call-option
		P_direct[i] = sumP_direct/nstep;	//media put-option
		//calcolo discretizzato
		C_discr[i] = sumC_discr/nstep;		//media call-option
		P_discr[i] = sumP_discr/nstep;		//media put-option
	}

	data_blocking(C_direct, "../Files/Cdirect.out");
	data_blocking(P_direct, "../Files/Pdirect.out");
	data_blocking(C_discr, "../Files/Cdiscr.out" );
	data_blocking(P_discr, "../Files/Pdiscr.out" );

	rnd.SaveSeed();
	return 0;
}
