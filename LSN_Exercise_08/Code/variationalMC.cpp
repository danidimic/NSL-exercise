#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "header.h"
#include "random.h"
#include "statistics.h"


int main (int argc, char *argv[]){

	Input();

	int naccept = 0;
	double current_pos = 0, proposed_pos = 0, alfa, c, xmax=3, xmin=-3;
	vec x(nstep), PsiT(nbins, fill::zeros);
	vector<double> H(nblock);
		
	for(int j=0; j<nblock; j++){
		for(int i=0; i<nstep; i++){
			proposed_pos = current_pos + rnd.Rannyu(-delta, delta);		//posizione proposta

			alfa = min( 1, pdf(proposed_pos) / pdf(current_pos) );
			if( rnd.Rannyu()<alfa ){	//accept
				current_pos = proposed_pos;
				naccept++;
			}
			x[i] = current_pos;
		}
		H[j] = MCintegral(x);					//calcolo l'integrale MC importance sampling

		//Definisco gli estremi dell'istogramma
		x.resize(nstep+2);
		x[nstep-2] = xmax;
		x[nstep-1] = xmin;
		PsiT += gethisto(x, nbins);		//istogramma delle posizioni
	}
	data_blocking(H, "results/Ene.out");
	cout<<"Accettazione complessiva Metropolis: "<<naccept/((double)nblock*nstep)<<endl;

	//Ottengo l'istogramma di PsiT calcolato tramite data blocking
	PsiT /= (double)nblock;
	ofstream Out("results/Psi.out");
	for(int i=0; i<(int)PsiT.size(); i++){
		c =  xmin + i*(xmax-xmin)/nbins;
		Out<<c<<"   "<<PsiT[i]<<endl;
	}
	Out.close();

/*
	//Valori da stampare per simulated annealing
	double ene = cumulative_average(H).back();
	ofstream Out("results/minimize.out");
	Out<<ene<<endl;
	Out.close();
*/

	rnd.SaveSeed();
	return 0;
}

//Input dei valori 
void Input(){

	rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	ifstream ReadInput("input.dat");
	ReadInput >> delta;
	cout << "Transizione uniforme di delta = " << delta<<endl;

	ReadInput >> mu;
	ReadInput >> sigma;
	cout<<"Gaussiane di parametri mu = " << mu << ", sigma = " << sigma << endl << endl;
}

//Integrale MC importance sampling
double MCintegral(vec x){
	int n = x.size();
	double s = 0, Eloc;
	for (int i=0; i<n; i++){
		Eloc = Hamiltonian(x[i]) / Psi_trial(x[i]);
		s += Eloc;
	}
	return s/(double)n;
}

//Ottieni funzione campionata come istogramma
vec gethisto(vec x, int nbins){
	uvec hist_psi = hist(x, nbins);
	
	vec histo(nbins), histo_norm(nbins);
	for (int i=0; i<nbins; i++)
		histo[i] = (double)hist_psi[i];

	histo_norm = normalise(histo);
	return histo_norm;
}

