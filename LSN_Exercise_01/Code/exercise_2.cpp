#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
using namespace std;

#define M 10000		//numero di valori per ciascun istogramma S_N


int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	int n;
	double s,e,l;
	vector<int> N = {1,2,10,100};
	vector<double> stand(4*M);
	vector<double> expo(4*M);
	vector<double> lorentz(4*M);		//vettori per immagazzinare i dati

	for(int i=0; i<M; i++){ 	//ciclo per M valori casuali di ciascun istogramma

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

	ofstream Sdata("../Files/exe2_stand.out");		//output per gli istogrammi dei valori uniformi
	ofstream Edata("../Files/exe2_expo.out");		//output per gli istogrammi dei valori distribuiti esponenzialmente
	ofstream Ldata("../Files/exe2_lorentz.out");	//output per gli istogrammi dei valori distribuiti lorentzianamente
	for(int i=0; i<M; i++){
		n=4*i;
		Sdata<<stand[n]<<"  "<<stand[n+1]<<"  "<<stand[n+2]<<"  "<<stand[n+3]<<endl;
		Edata<<expo[n]<<"  "<<expo[n+1]<<"  "<<expo[n+2]<<"  "<<expo[n+3]<<endl;
		Ldata<<lorentz[n]<<"  "<<lorentz[n+1]<<"  "<<lorentz[n+2]<<"  "<<lorentz[n+3]<<endl;
	}
	Sdata.close();
	Edata.close();
	Ldata.close();

	rnd.SaveSeed();

	return 0;
}
