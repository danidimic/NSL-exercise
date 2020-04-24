#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "statistics.h"
using namespace std;

#define M 10000  //numero di prezzi degli asset
#define N 100   //numero di blocchi

#define S0 100	//prezzo iniziale
#define T 1		//tempo finale
#define K 100	//prezzo di mercato
#define r 0.1		//interesse privo di rischio
#define sigma 0.25	//volatilità


double max(double a, double b){
	double Max;
	if (a>b) Max=a;
	else Max=b;
	return Max;
}

int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	double x, z, S;
	double sumC, sumP;
	vector<double> C(N);  //call-option e relative incertezze
	vector<double> P(N);  //put-option e relative incertezze


	//--------------------------------------------------------------------------//
	// Esercizio 03.1 parte 1 : stima diretta di S(T)							//
	//--------------------------------------------------------------------------//

	for(int i=0; i<N; i++){
		sumC = 0;
		sumP = 0;
		for(int j=0; j<M; j++){
			z = rnd.Gauss(0, T);
			S = S0*exp( (r-pow(sigma,2)/2)*T + sigma*z );	//simulazione diretta del prezzo S(T) 

			x = exp(-r*T)*max(0, S-K);	//profitto call-option
			sumC += x;

			x = exp(-r*T)*max(0, K-S);	//profitto put-option
			sumP += x;
		}
		C[i] = sumC/M;	//media call-option
		P[i] = sumP/M;	//media put-option
	}

	vector<double> ave_progC = cumulative_average( C );	//media progressiva call-option
	vector<double> err_progC = cumulative_error( C );		//errore progressivo call-option

	vector<double> ave_progP = cumulative_average( P );	//media progressiva call-option
	vector<double> err_progP = cumulative_error( P );		//errore progressivo call-option

	ofstream data("../Files/exe1_direct.out");		//genero l'output per il calcolo diretto
	if (data.is_open()){
		for(int i=0; i<N; i++){
			data<<i+1<<"  "<<ave_progC[i]<<"  "<<err_progC[i]<<"  "<<ave_progP[i]<<"  "<<err_progP[i]<<endl;
		}
	} else cerr << "PROBLEM: Unable to open exe1_direct.out" << endl;
	data.close();


	//--------------------------------------------------------------------------//
	// Esercizio 03.1 parte 2 : stima discretizzata di S(T)						//
	//--------------------------------------------------------------------------//

	double t, t1;

	for(int i=0; i<N; i++){
		sumC = 0;
		sumP = 0;
		for(int j=0; j<M; j++){
			S = S0;
			t = 0;
			for(int k=0; k<N; k++){		//simulazione discretizzata del prezzo S(T) 
				t1 = t + 1./N;
				z = rnd.Gauss(0,1);
				S = S*exp( (r-pow(sigma,2)/2)*(t1-t) + sigma*z*sqrt(t1-t) );
				t = t1;
			}
			x = exp(-r*T)*max(0, S-K);	//profitto call-option
			sumC += x;

			x = exp(-r*T)*max(0, K-S);	//profitto put-option
			sumP += x;
		}
		C[i] = sumC/M;	//media call-option
		P[i] = sumP/M;	//media put-option
	}

	ave_progC = cumulative_average( C );	//media progressiva call-option
	err_progC = cumulative_error( C );		//errore progressivo call-option

	ave_progP = cumulative_average( P );	//media progressiva call-option
	err_progP = cumulative_error( P );		//errore progressivo call-option

	data.open("../Files/exe1_discretiz.out");		//genero l'output per il calcolo discretizzato
	if (data.is_open()){
		for(int i=0; i<N; i++){
			data<<i+1<<"  "<<ave_progC[i]<<"  "<<err_progC[i]<<"  "<<ave_progP[i]<<"  "<<err_progP[i]<<endl;
		}
	} else cerr << "PROBLEM: Unable to open exe1_discretiz.out" << endl;
	data.close();

	rnd.SaveSeed();
	return 0;
}
