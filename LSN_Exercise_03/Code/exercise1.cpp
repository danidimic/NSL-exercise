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

double max(double, double);


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

	//medie progressive simulazione diretta
	vector<double> aveC_direct = cumulative_average( C_direct );	//media progressiva call-option
	vector<double> errC_direct = cumulative_error(   C_direct );	//errore progressivo call-option
	vector<double> aveP_direct = cumulative_average( P_direct );	//media progressiva call-option
	vector<double> errP_direct = cumulative_error(   P_direct );	//errore progressivo call-option
	//medie progressive simulazione discretizzata
	vector<double> aveC_discr = cumulative_average( C_discr );		//media progressiva call-option
	vector<double> errC_discr = cumulative_error(   C_discr );		//errore progressivo call-option
	vector<double> aveP_discr = cumulative_average( P_discr );		//media progressiva call-option
	vector<double> errP_discr = cumulative_error(   P_discr );		//errore progressivo call-option

	ofstream Direct("../Files/direct.out");		//output calcolo diretto
	ofstream Discr ("../Files/discr.out");		//output calcolo discretizzato
	for(int i=0; i<nblock; i++){
		Direct<<aveC_direct[i]<<"  "<<errC_direct[i]<<"  "<<aveP_direct[i]<<"  "<<errP_direct[i]<<endl;
		Discr<<aveC_discr[i]<<"  "<<errC_discr[i]<<"  "<<aveP_discr[i]<<"  "<<errP_discr[i]<<endl;
	}
	Direct.close();
	Discr.close();

	rnd.SaveSeed();
	return 0;
}
//massimo tra due valori
double max(double a, double b){
	double Max;
	if (a>b) Max=a;
	else Max=b;
	return Max;
}
