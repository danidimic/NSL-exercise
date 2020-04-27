#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "statistics.h"
using namespace std;


int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	int nrand = 1E6;	//valori casuali generati
	int nblock = 100;	//numero di blocchi
	int nvalues = nrand/nblock;	//valori in ciascun blocco

	int nint = 100;		//intervalli di [0,1] per calcolo chi2

	double sum1, sum2, r, min, max;
	vector<double> ave(nblock);		//vettore per calcolo dei val medi
	vector<double> var(nblock);		//vettore per calcolo varianza
	vector<double> oss(nint);		//vettore per calcolo chi quadro
	vector<double> chi2(nblock);				//vettore per il chi quadro

	for(int i=0; i<nblock; i++){
		sum1=0;
		sum2=0;

		for(int j=0; j<nint; j++) oss[j]=0;	//azzera il contatore di valori osservati

		for(int j=0; j<nvalues; j++) {
			r = rnd.Rannyu();
			sum1 += r;					//calcolo del valore medio
			sum2 += pow( (r-0.5), 2);	//calcolo della varianza

			for(int l=0; l<nint; l++){		//determino a che intervallo appartiene il valore casuale r
				min = l*1./(double)nint;			//minimo dell'intervallo l-esimo
				max = (l+1)*1./(double)nint;		//massimo dell'intervallo l-esimo

				if( r >= min && r < max){
					oss[l]++;
					break;
				}
			}
			//medie sui blocchi
			ave[i] = sum1/nvalues;
			var[i] = sum2/nvalues;
			//calcolo del chi quadro
			chi2[i] = chi_quadro(oss, (double)nvalues/nint );
		}
	}

	data_blocking(ave, "../Files/ave.out");
	data_blocking(var, "../Files/var.out");

	ofstream Chi;
	Chi.open("../Files/chi.out");
	for(int i=0; i<nblock; i++) 	Chi<<chi2[i]<<endl;

	rnd.SaveSeed();
	return 0;
}
