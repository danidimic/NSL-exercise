#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "statistics.h"
using namespace std;

//funzione integranda
double integranda(double x){
	return M_PI/2*cos( M_PI*x/2);
}
//distribuzione di probabilità per importance sampling
double pdf(double y){
	return 1-sqrt(1-y);
}


int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();

	int nblock = 100, nrand = 1E6;
	double a = 0, b = 1;	//estremi di integrazione

	double r, sum_uniform, sum_importance;
	vector<double> I_uniform(nblock);		//integrale uniform sampling
	vector<double> I_importance(nblock);		//integrale importance sampling
	
	for(int i=0; i<nblock; i++){
		sum_uniform = 0;
		sum_importance = 0;
		for(int j=0; j<nrand; j++){
			//uniform sampling
			sum_uniform += integranda( rnd.Rannyu() );
			//importance sampling
			r = pdf(rnd.Rannyu());
			sum_importance += integranda(r)/(2*(1-r));
		}
		I_uniform[i] = (b-a)*sum_uniform/nrand;
		I_importance[i] = (b-a)*sum_importance/nrand;
	}

	vector<double> ave1(nblock), err1(nblock);
	vector<double> ave2(nblock), err2(nblock);

	ave1 = cumulative_average(I_uniform);		//media cumulativa integrale uniform sampling
	err1 = cumulative_error(I_uniform);			//errore statistico integrale uniform sampling
	ave2 = cumulative_average(I_importance);	//media cumulativa integrale importance sampling
	err2 = cumulative_error(I_importance);		//errore statistico integrale importance sampling

	ofstream MCuniform("../Files/uniform.out");		//genero l'output su un unico file di testo
	ofstream MCimportance("../Files/importance.out");	
	for(int i=0; i<nblock; i++){
		MCuniform<<setprecision(9)<<ave1[i]<<"  "<<err1[i]<<endl;
		MCimportance<<setprecision(9)<<ave2[i]<<"  "<<err2[i]<<endl;
	}
	MCuniform.close();
	MCimportance.close();

	rnd.SaveSeed();
	return 0;
}
