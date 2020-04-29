#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "header.h"
#include "random.h"
#include "statistics.h"


int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali
	
	int nthrows = 1E7, nblock = 500, nvalues = nthrows/nblock;
	int nacceptU100 = 0, nacceptU210 = 0, nacceptN100 = 0, nacceptN210 = 0;

	double sU100, sU210, sN100, sN210;
	double alfaU100, alfaU210, alfaN100, alfaN210;
	double Du100 = 1.2, Du210 = 3.15;	//delta transizione uniforme
	double Dg100 = 0.75, Dg210 = 1.95;	//delta transizione gaussiana

	//Transizione uniforme
	vec current_unif100(3), proposed_unif100(3);	//orbitale Psi100
	vec current_unif210(3), proposed_unif210(3);	//orbitale Psi210
	//punto di partenza sul massimo della pdf
	current_unif100.fill(1.5/sqrt(3));
	current_unif210.fill(5.0/sqrt(3));

	//Transizione normale
	vec current_norm100(3), proposed_norm100(3);	//orbitale Psi100
	vec current_norm210(3), proposed_norm210(3);	//orbitale Psi210
	//punto di partenza sul massimo della pdf
	current_norm100.fill(1.5/sqrt(3));
	current_norm210.fill(5.0/sqrt(3));

	vec Tuniform100(3), Tuniform210(3);		//Probabilità transizione uniforme
	vec Tnormal100(3),  Tnormal210(3);		//Probabilità transizione normale

	vector<double> rU100(nblock), rU210(nblock), rN100(nblock), rN210(nblock);

	for(int j=0; j<nblock; j++){
		sU100 = 0;
		sU210 = 0;
		sN100 = 0;
		sN210 = 0;
		for(int i=0; i<nvalues; i++){
			//PROBABILITA UNIFORME
			Tuniform100.imbue( [&]() {return rnd.Rannyu(-Du100, Du100);} );	
			Tuniform210.imbue( [&]() {return rnd.Rannyu(-Du210, Du210);} );
			//posizione proposta
			proposed_unif100 = current_unif100 + Tuniform100;	//posizione proposta per Psi n=1, l=0, m=0
			proposed_unif210 = current_unif210 + Tuniform210;	//posizione proposta per Psi n=2, l=1, m=0
			//calcolo alfa per accept-reject
			alfaU100 = min(1, pdf(Psi100, proposed_unif100)/pdf(Psi100, current_unif100) );
			alfaU210 = min(1, pdf(Psi210, proposed_unif210)/pdf(Psi210, current_unif210) );

			if( rnd.Rannyu()<alfaU100 ){	//accept
				current_unif100 = proposed_unif100;
				nacceptU100++;
			}
			if( rnd.Rannyu()<alfaU210 ){	//accept
				current_unif210 = proposed_unif210;
				nacceptU210++;
			}
			sU100 += norm(current_unif100);		//r dei vettori per Psi n=1, l=0, m=0 
			sU210 += norm(current_unif210);		//r dei vettori per Psi n=2, l=1, m=0

			//PROBABILITA NORMALE
			Tnormal100.imbue( [&]() {return rnd.Gauss(0, Dg100);} );	
			Tnormal210.imbue( [&]() {return rnd.Gauss(0, Dg210);} );	
			//posizione proposta
			proposed_norm100 = current_norm100 + Tnormal100;	//posizione proposta per Psi n=1, l=0, m=0
			proposed_norm210 = current_norm210 + Tnormal210;	//posizione proposta per Psi n=2, l=1, m=0
			//calcolo alfa per accept-reject
			alfaN100 = min(1, pdf(Psi100, proposed_norm100)/pdf(Psi100, current_norm100) );
			alfaN210 = min(1, pdf(Psi210, proposed_norm210)/pdf(Psi210, current_norm210) );

			if( rnd.Rannyu()<alfaN100 ){	//accept
				current_norm100 = proposed_norm100;
				nacceptN100++;
			}
			if( rnd.Rannyu()<alfaN210 ){	//accept
				current_norm210 = proposed_norm210;
				nacceptN210++;
			}
			sN100 += norm(current_norm100);		//r dei vettori per Psi n=1, l=0, m=0 
			sN210 += norm(current_norm210);		//r dei vettori per Psi n=2, l=1, m=0

		}
		rU100[j] = sU100/nvalues;
		rU210[j] = sU210/nvalues;
		rN100[j] = sN100/nvalues;
		rN210[j] = sN210/nvalues;
	}

	//risultati medie progressive su file (PROBABILITA UNIFORME)
	data_blocking(rU100, "../Files/unif100.out");
	data_blocking(rU210, "../Files/unif210.out");

	//risultati medie progressive su file (PROBABILITA NORMALE)
	data_blocking(rN100, "../Files/gauss100.out");
	data_blocking(rN210, "../Files/gauss210.out");
	
	cout<<"Efficenza transizione uniforme n=1, l=0, m=0 = "<<nacceptU100/(double)nthrows<<endl;
	cout<<"Efficenza transizione uniforme n=2, l=1, m=0 = "<<nacceptU210/(double)nthrows<<endl;
	cout<<"Efficenza transizione normale n=1, l=0, m=0 = "<<nacceptN100/(double)nthrows<<endl;
	cout<<"Efficenza transizione normale n=2, l=1, m=0 = "<<nacceptN210/(double)nthrows<<endl;
	
	rnd.SaveSeed();
	return 0;
}

