#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;


//calcolo delle medie progressive
vector<double> cumulative_average( vector<double> v ){

	int N = v.size();
	vector<double> mp(N);
	for(int i=0; i<N; i++){
		mp[i] = 0;
		for(int j=0; j<=i; j++){
			mp[i] += v[j];
		}
		mp[i] /= (i+1);
	}
	return mp;
}
//calcolo delle medie progressive dei quadrati
vector<double> cumulative_average_quad( vector<double> v ){

	int N = v.size();
	vector<double> mpq(N);
	for(int i=0; i<N; i++){
		mpq[i] = 0;
		for(int j=0; j<=i; j++){
			mpq[i] += pow(v[j],2);
		}
		mpq[i] /= (i+1);
	}
	return mpq;
}
//calcolo delle incertezze statistiche sulle medie progressive
vector<double> cumulative_error( vector<double> v ){

	int N = v.size();
	vector<double> mp(N);
	mp = cumulative_average(v);	//vettore con le medie progressive
	vector<double> mpq(N);
	mpq = cumulative_average_quad(v);  //vettore con le medie progressive dei quadrati
	vector<double> err(N);

	for(int i=0; i<N; i++){
		err[i] = sqrt( mpq[i]-pow(mp[i],2) );
		err[i] /= sqrt(i+1);
	}
	return err;
}
//risultati data blocking su file
void data_blocking(vector<double> results, string path){

	int nblock = results.size();
	vector<double> ave(nblock), err(nblock);
	ave = cumulative_average(results);
	err = cumulative_error(results);
	
	ofstream Output;
	Output.open(path);
	for(int i=0; i<nblock; i++){
		Output<<setprecision(9)<<ave[i]<<"  "<<err[i]<<endl;
	}
	Output.close();

	return;
}
//calcolo del chi quadrato
double chi_quadro( vector<double> oss, double E ){

	int N = oss.size();
	double chi = 0;
	for(int i=0; i<N; i++){
		chi += pow( (oss[i]-E),2 )/E;
	}
	return chi;

}
