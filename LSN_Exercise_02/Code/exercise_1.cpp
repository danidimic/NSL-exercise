#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "statistics.h"
using namespace std;

#define N 100		//numero di blocchi
#define M 10E4		//valori per il calcolo dell'integrale

#define a 0		//estremo inferiore intervallo integrazione
#define b 1		//estremo superiore intervallo integrazione


double integranda(double x){	//funzione integranda
	return M_PI/2*cos( M_PI*x/2);
}

double distr_probab(double y){	//distribuzione di probabilità per importance sampling
	return 1-sqrt(1-y);
}

int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	double r;
	double sum_uniform, sum_importance;
	vector<double> I_uniform(N);		//risultati dell'integrale uniform sampling
	vector<double> I_importance(N);		//risultati dell'integrale importance sampling
	
	for(int i=0; i<N; i++){
		sum_uniform = 0;
		sum_importance = 0;
		for(int j=0; j<M; j++){
			sum_uniform += integranda( rnd.Rannyu() );	//uniform sampling

			r = distr_probab( rnd.Rannyu()) ;	//importance sampling
			sum_importance += integranda(r)/(2*(1-r));
		}
		I_uniform[i] = (b-a)*sum_uniform/M;
		I_importance[i] = (b-a)*sum_importance/M;
	}

	vector<double> ave1(N), err1(N);  	//vettori per immagazzinare i risultati uniform sampling
	vector<double> ave2(N), err2(N);  	//vettori per immagazzinare i risultati importance sampling

	ave1 = cumulative_average(I_uniform);	//media cumulativa integrale uniform sampling
	err1 = cumulative_error(I_uniform);				//errore statistico integrale uniform sampling

	ave2 = cumulative_average(I_importance);	//media cumulativa integrale importance sampling
	err2 = cumulative_error(I_importance);			//errore statistico integrale importance sampling

	ofstream data("../Files/exe1.out");		//genero l'output su un unico file di testo
	if (data.is_open()){
		for(int i=0; i<N; i++){
			data<<i+1<<"  "<<ave1[i]<<"  "<<err1[i]<<"  "<<ave2[i]<<"  "<<err2[i]<<endl;
		}
	} else cerr << "PROBLEM: Unable to open exe1.out" << endl;
	data.close();

	rnd.SaveSeed();
	return 0;
}
