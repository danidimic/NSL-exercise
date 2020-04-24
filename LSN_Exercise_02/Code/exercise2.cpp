#include <iostream>
#include <fstream>
#include <string>
#include "header.h"
#include "random.h"
#include "statistics.h"


int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	int nstep = 100, nrandwalk = 1E4, r;
	double* angles = new double[3];

	//Random walk cubic lattice
	vector<double> RWcube(3);	//posizione RW
	vector<double> mean_cube(nstep), err_cube(nstep);
	vector<double> sum_cube(nstep), sum2_cube(nstep);
	//Random walk continuum
	vector<double> RWcont(3);	//posizione RW
	vector<double> mean_cont(nstep), err_cont(nstep);
	vector<double> sum_cont(nstep), sum2_cont(nstep);

	for(int i=0; i<nrandwalk; i++){
		for(int k=0; k<3; k++){	//entrambi i random walk partono dall'origine
			RWcube[k] = 0;
			RWcont[k] = 0;
		}

		for(int j=0; j<nstep; j++){
			//Random walk cubic lattice
			r = rnd.Rannyu(0,3);	//direzione di spostamento casuale
			if( rnd.Rannyu()<0.5 ) 
				RWcube[r] += a;
			else
				RWcube[r] -= a;

			sum_cube[j] += pow( modulo(RWcube), 2 );
			sum2_cube[j] += pow( modulo(RWcube), 4 );

			//Random walk continuum
			angles = rnd.Sphere();
			RWcont[0] += X(a, angles[0], angles[1]);
			RWcont[1] += Y(a, angles[0], angles[1]);
			RWcont[2] += Z(a, angles[0], angles[1]);

			sum_cont[j] += pow( modulo(RWcont), 2 );
			sum2_cont[j] += pow( modulo(RWcont), 4 );
		}
	}

	for(int i=0; i<nstep; i++){
		//Random walk cubic lattice
		sum_cube[i] = sum_cube[i]/nrandwalk;		//media dei moduli quadri
		sum2_cube[i] = sum2_cube[i]/nrandwalk;		//media dei quadrati dei moduli quadri
		//Random walk continuum
		sum_cont[i] = sum_cont[i]/nrandwalk;		//media dei moduli quadri
		sum2_cont[i] = sum2_cont[i]/nrandwalk;		//media dei quadrati dei moduli quadri

		//cout<<sum_cube[i]<<"  "<<sum_cont[i]<<endl;
	}

	ofstream Cube("../Files/RWlattice.out");
	ofstream Cont("../Files/RWcontinuum.out");
	for(int i=0; i<nstep; i++){
		Cube<<sqrt(sum_cube[i])<<"  "<<error(sum_cube[i], sum2_cube[i], nrandwalk)<<endl;
		Cont<<sqrt(sum_cont[i])<<"  "<<error(sum_cont[i], sum2_cont[i], nrandwalk)<<endl;
	}
	Cube.close();
	Cont.close();

	rnd.SaveSeed();

	return 0;
}
