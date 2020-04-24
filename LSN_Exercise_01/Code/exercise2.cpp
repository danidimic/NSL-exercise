#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
using namespace std;


int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	int n, nvalues = 1E5;
	double s,e,l;
	vector<int> N = {1,2,10,100};

	vector<double> stand(4*nvalues);			//dati per S_N uniforme
	vector<double> expo(4*nvalues);			//dati per S_N esponenziale
	vector<double> lorentz(4*nvalues);		//dati per S_N lorentziani

	for(int i=0; i<nvalues; i++){ 	//ciclo per nvalues valori casuali di ciascun istogramma

		for(int j=0; j<4; j++){		//ciclo sui possibili valori N=1,2,10,100
			n = 4*i+j;
			s = 0;			//somme sui valori uniformi
			e = 0;			//somme sui valori esponenziali
			l = 0;			//somme sui valori lorentziani

			for(int k=0; k<N[j]; k++){		//ciclo per il calcolo delle medie degli N valori casuali
				s += rnd.Rannyu();	//valori distribuiti in modo uniforme
				e += rnd.Expo(1);	//valori distribuiti in modo esponenziale	
				l += rnd.Lorentz(0,1);  //valori distribuiti in modo lorentziano
			}
			stand[n] = s/N[j];		//media di N[j] valori casuali uniformi
			expo[n] = e/N[j];		//media di N[j] valori casuali distribuiti esponenzialmente
			lorentz[n] = l/N[j];	//media di N[j] valori casuali distribuiti lorentzianamente
		}
	}

	ofstream Stand("../Files/SNstand.out");		//output per gli istogrammi dei valori uniformi
	ofstream Expo("../Files/SNexpo.out");		//output per gli istogrammi dei valori distribuiti esponenzialmente
	ofstream Lore("../Files/SNlorentz.out");	//output per gli istogrammi dei valori distribuiti lorentzianamente
	for(int i=0; i<nvalues; i++){
		n=4*i;
		Stand<<stand[n]<<"  "<<stand[n+1]<<"  "<<stand[n+2]<<"  "<<stand[n+3]<<endl;
		Expo<<expo[n]<<"  "<<expo[n+1]<<"  "<<expo[n+2]<<"  "<<expo[n+3]<<endl;
		Lore<<lorentz[n]<<"  "<<lorentz[n+1]<<"  "<<lorentz[n+2]<<"  "<<lorentz[n+3]<<endl;
	}
	Stand.close();
	Expo.close();
	Lore.close();

	rnd.SaveSeed();

	return 0;
}
