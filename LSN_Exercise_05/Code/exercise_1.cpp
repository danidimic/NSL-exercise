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
	int naccept100 = 0, naccept210 = 0;
	double alfa100, alfa210, s100, s210;

	double Du100 = 1.2, Du210 = 3.15;	//delta transizione uniforme
	double Dg100 = 0.75, Dg210 = 1.95;	//delta transizione gaussiana

	vec current_pos100(3), proposed_pos100(3);	//vettori per orbitale Psi100
	vec current_pos210(3), proposed_pos210(3);	//vettori per orbitale Psi210
	vec transition100(3), transition210(3);		//vettori per probabilità transizione
	//punto di partenza sul massimo della pdf
	current_pos100.fill(1.5/sqrt(3));
	current_pos210.fill(5.0/sqrt(3));

	vector<double> r100(nblock), r210(nblock);	//vettori per le distanze r dei punti dall'origine

	for(int j=0; j<nblock; j++){
		s100 = 0;
		s210 = 0;
		for(int i=0; i<nvalues; i++){
			//probabilità uniforme
			//transition100.imbue( [&]() {return rnd.Rannyu(-Du100, Du100);} );	
			//transition210.imbue( [&]() {return rnd.Rannyu(-Du210, Du210);} );
			//probabilità gaussiana	
			transition100.imbue( [&]() {return rnd.Gauss(0, Dg100);} );	
			transition210.imbue( [&]() {return rnd.Gauss(0, Dg210);} );	
			
			proposed_pos100 = current_pos100 + transition100;	//posizione proposta per Psi n=1, l=0, m=0
			proposed_pos210 = current_pos210 + transition210;	//posizione proposta per Psi n=2, l=1, m=0

			alfa100 = min(1, pdf(Psi100, proposed_pos100)/pdf(Psi100, current_pos100) );
			alfa210 = min(1, pdf(Psi210, proposed_pos210)/pdf(Psi210, current_pos210) );

			if( rnd.Rannyu()<alfa100 ){	//accept
				current_pos100 = proposed_pos100;
				naccept100++;
			}
			if( rnd.Rannyu()<alfa210 ){	//accept
				current_pos210 = proposed_pos210;
				naccept210++;
			}

			s100 += norm(current_pos100);		//r dei vettori per Psi n=1, l=0, m=0 
			s210 += norm(current_pos210);		//r dei vettori per Psi n=2, l=1, m=0
		}
		r100[j] = s100/nvalues;
		r210[j] = s210/nvalues;
	}

	//risultati medie progressive su file (probabilità uniforme)
	//data_blocking(r100, "../Files/unif100.out");
	//data_blocking(r210, "../Files/unif210.out");

	//risultati medie progressive su file (probabilità uniforme)
	data_blocking(r100, "../Files/gauss100.out");
	data_blocking(r210, "../Files/gauss210.out");
	
	cout<<"Efficenza n=1, l=0, m=0 = "<<naccept100/(double)nthrows<<endl;
	cout<<"Efficenza n=2, l=1, m=0 = "<<naccept210/(double)nthrows<<endl;
	
	rnd.SaveSeed();
	return 0;
}
