#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "statistics.h"
using namespace std;

#define N 100		//numero di blocchi
#define R 10000		//numero di semi lanciati in ogni blocco

#define d 0.75		//spaziatura tra le righe del campo
#define L 0.5		//lunghezza dell'ago


int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali
	
	int M;
	double x,y;		//(x,y) coppia di punti per generare un angolo uniforme in (0,pi/2)
	double z;		//z = distanza dall riga più vicina
	double theta;	//theta = angolo acuto formato con le righe
	vector<double> pi(N);	//valori di pi_greco ottenuti sui singoli blocchi

	for(int i=0; i<N; i++){
		M=0;
		for(int j=0; j<R; j++){
			z = rnd.Rannyu(0, d/2);	//distanza dalla riga più vicina

			while ( true ){		
				x = rnd.Rannyu(0,1);
				y = rnd.Rannyu(0,1);
				if ( pow(x,2)+pow(y,2) <= 1 ){	//rigetto la coppia (x,y) se non è interna alla circonferenza
					theta = atan( y/x );
					break;
				}
			}

			if( z <= 0.5*L*sin(theta) ) M++;	//condizione per avere intersezione tra semi e linee
		}
		pi[i] = (2*L*R)/(M*d);
	}

	vector<double> ave(N), err(N);  //vettori necessari per immagazzinare i risultati finali
	ave = media_progressiva(pi);	//medie progressive sui blocchi
	err = errore(pi);	//incertezze statistiche sui singoli blocchi

	ofstream data("../Files/exe3.out");		//genero l'output su un unico file di testo
	if (data.is_open()){
		for(int i=0; i<N; i++){
			data<<i+1<<"  "<<ave[i]<<"  "<<err[i]<<endl;
		}
	} else cerr << "PROBLEM: Unable to open exe3.out" << endl;
	data.close();


	rnd.SaveSeed();
	return 0;
}
