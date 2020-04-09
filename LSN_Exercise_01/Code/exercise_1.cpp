#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "statistics.h"
using namespace std;

#define N 100		//numero di blocchi
#define n 10E4		//valori per il calcolo del chi quadrato
#define M 10E5		//valori per il calcolo di media e varianza


int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	double sum;
	int L = M/N;		//numero di valori in ciascun blocco
	vector<double> ave(N);	//vettore per le medie ottenute sui singoli blocchi

	vector<double> r1(N), r1_err(N);	//vettori per i risultati della parte 1
	vector<double> r2(N), r2_err(N);	//vettori per i risultati della parte 2

	//--------------------------------------------------------------------------//
	// Esercizio 01.1 parte 1 : calcolo del valore medio e della sua incertezza //
	//--------------------------------------------------------------------------//

	for(int i=0; i<N; i++){
		sum=0;
		for(int j=0; j<L; j++) {
			sum += rnd.Rannyu();
		}
		ave[i] = sum/L;		//valore medio degli L elementi di ciascun blocco
	}

	r1 = media_progressiva(ave);	//medie progressive
	r1_err = errore(ave);		//incertezze statistica sulle medie progressive

	//----------------------------------------------------------------------------//
	// Esercizio 01.1 parte 2 : calcolo della varianza			  				  //
	//----------------------------------------------------------------------------//

	for(int i=0; i<N; i++){
		sum=0;
		for(int j=0; j<L; j++){
			sum += pow( (rnd.Rannyu()-0.5), 2);
		}
		ave[i] = sum/L;		//valore medio degli L elementi di ciascun blocco
	}

	r2 = media_progressiva(ave);	//medie progressive
	r2_err = errore(ave);		//incertezze statistica sulle medie progressive

	//----------------------------------------------------------------------------//
	// Esercizio 01.1 parte 3 : test chi quadrato								  //
	//----------------------------------------------------------------------------//

	int R = 100;		//intervalli di [0,1] con probabilità 1/R
	double r, max, min;
	vector<double> oss(R);
	vector<double> chi2(N);

	for(int i=0; i<N; i++){		//ciclo per calcolo di N chi2
		
		for(int j=0; j<R; j++) oss[j]=0;	//azzera il contatore di valori osservati

		for(int k=0; k<n; k++){		//genero n valori casuali per ciascun blocco
			r = rnd.Rannyu();
			for(int l=0; l<R; l++){		//determino a che intervallo appartiene il valore
				min = l*1./R;		//minimo dell'intervallo
				max = (l+1)*1./R;	//massimo dell'intervallo

				if( r >= min && r < max){ 
					oss[l]++;
					break;
				}
			}
		}
		chi2[i] = chi_quadro(oss, (double)n/R );	//calcolo il chi2 tramite i valori osservati	
	}


	ofstream data("../Files/exe1.out");		//genero l'output su un unico file di testo
	if (data.is_open()){
		for(int i=0; i<N; i++){
			data<<i*L<<"  "<<r1[i]<<"  "<<r1_err[i]<<"  ";		//dati relativi alla parte 1
			data<<r2[i]<<"  "<<r2_err[i]<<"  ";					//dati relativi alla parte 2
			data<<i<<"  "<<chi2[i]<<endl;						//dati relativi alla parte 3
		}
	} else cerr << "PROBLEM: Unable to open exe1.out" << endl;
	data.close();

	rnd.SaveSeed();
	return 0;
}
