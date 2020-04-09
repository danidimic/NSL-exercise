#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "statistics.h"
using namespace std;

#define N 10000		//numero di random walk simulati
#define M 100		//numero di passi in ogni random walk

#define d 3			//numero di dimensioni
#define a 1			//costante reticolare


//coordinate sferiche
double X(double theta, double phi){		//coordinata x
	return sin(theta)*sin(phi);
}

double Y(double theta, double phi){		//coordinata y
	return sin(theta)*cos(phi);
}

double Z(double theta, double phi){		//coordinata z
	return cos(theta);
}

//incertezze statistiche
double error(double sum, double sum2){
	double e = sqrt( (sum2-pow(sum,2))/(N-1) );
	return e/(2*sqrt(sum));	
}

int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	int r;
	double p;
	double* angles = new double[2];		//theta e phi in coordinate sferiche
	
	vector<double> pos(d);				//vettore d-dimensionale per la posizione dell'oggetto
	vector<double> mean(M), err(M);		//medie e incertezze da plottare
	vector<double> sum(M), sum2(M);


	//--------------------------------------------------------------------------//
	// Esercizio 02.2 parte 1 : random walk su reticolo cubico					//
	//--------------------------------------------------------------------------//


	for(int i=0; i<N; i++){
		for(int k=0; k<d; k++)	//random walk parte dall'origine
			pos[k] = 0;

		for(int j=0; j<M; j++){

			r = rnd.Rannyu(0,d);	//direzione di spostamento casuale
			p = rnd.Rannyu();		//spostamento "avanti" o "indietro" nella direzione r
			if( p<0.5 ) 
				pos[r] += a;
			else
				pos[r] -= a;

			sum[j] += pow( modulo(pos), 2 );
			sum2[j] += pow( modulo(pos), 4 );
		}
	}

	for(int i=0; i<M; i++){
		sum[i] = sum[i]/N;		//media dei moduli quadri
		sum2[i] = sum2[i]/N;	//media dei quadrati dei moduli quadri
	}

	ofstream data("../Files/exe2_lattice.out");		//genero l'output relativo al random walk nel reticolo
	if (data.is_open()){
		for(int i=0; i<M; i++){
			data<<i+1<<"  "<<sqrt(sum[i])<<"  "<<error(sum[i], sum2[i])<<endl;
		}
	} else cerr << "PROBLEM: Unable to open exe2_lattice.out" << endl;
	data.close();


	//--------------------------------------------------------------------------//
	// Esercizio 02.2 parte 2 : random walk su reticolo continuo				//
	//--------------------------------------------------------------------------//

	sum.clear();
	sum.resize(M, 0);
	sum2.clear();
	sum2.resize(M, 0);

	for(int i=0; i<N; i++){
		for(int k=0; k<d; k++)	//random walk parte dall'origine
			pos[k] = 0;

		for(int j=0; j<M; j++){

			angles = rnd.Sphere();
			pos[0] += X(angles[0], angles[1]);
			pos[1] += Y(angles[0], angles[1]);
			pos[2] += Z(angles[0], angles[1]);

			sum[j] += pow( modulo(pos), 2 );
			sum2[j] += pow( modulo(pos), 4 );
		}
	}

	for(int i=0; i<M; i++){
		sum[i] = sum[i]/N;		//media dei moduli quadri
		sum2[i] = sum2[i]/N;	//media dei quadrati dei moduli quadri
	}

	data.open("../Files/exe2_continuum.out");		//genero l'output relativo al random walk nel continuo
	if (data.is_open()){
		for(int i=0; i<M; i++){
			data<<i+1<<"  "<<sqrt(sum[i])<<"  "<<error(sum[i], sum2[i])<<endl;
		}
	} else cerr << "PROBLEM: Unable to open exe2_continuum.out" << endl;
	data.close();


	rnd.SaveSeed();

	return 0;
}
