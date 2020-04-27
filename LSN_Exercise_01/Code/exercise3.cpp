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
	
	int nblock = 100, nthrow = 1E4, M;
	double d = 0.75;	//spaziatura tra le righe del campo
	double L = 0.5;		//lunghezza dell'ago
	double x, y, z, theta;
	vector<double> pi(nblock);	//valori di pi_greco ottenuti sui singoli blocchi

	for(int i=0; i<nblock; i++){
		M=0;
		for(int j=0; j<nthrow; j++){
			z = rnd.Rannyu(0, d/2);	//distanza dalla riga più vicina

			while ( true ){	//genero l'angolo in modo uniforme		
				x = rnd.Rannyu(0,1);
				y = rnd.Rannyu(0,1);
				if ( pow(x,2)+pow(y,2) <= 1 ){	//rigetto la coppia (x,y) se non è interna alla circonferenzaunitaria
					theta = atan( y/x );
					break;
				}
			}

			if( z <= 0.5*L*sin(theta) ) M++;	//condizione per avere intersezione tra semi e linee
		}
		pi[i] = (2*L*nthrow)/(M*d);
	}

	data_blocking(pi, "../Files/pi.out");

	rnd.SaveSeed();
	return 0;
}
