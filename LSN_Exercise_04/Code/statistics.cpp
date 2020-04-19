#include <vector>
#include <cmath>
using namespace std;

//calcolo del modulo di un vettore
double modulo(vector<double> vet){

	int N = vet.size();
	double s = 0;
	for(int i=0; i<N; i++){
		s += pow( vet[i],2 );
	}
	return pow(s,0.5);
}
//calcolo della media
double media(vector<double> v){

	int N = v.size();
	double sum=0;
	for(int i=0; i<N; i++){
		sum += v[i];
	}
	return sum/N;
}
//calcolo della deviazione standard
double devstd(vector<double> v){

	int N = v.size();
	double sum = 0;
	double mean = media(v);
	for(int i=0; i<N; i++){
		sum += pow( v[i]-mean,2 );
	}
	return sqrt( sum/N );
}
//calcolo della deviazione standard della media
double devstd_media(vector<double> v){
	int N = v.size();
	return devstd(v)/sqrt(N);
}
//calcolo delle medie progressive
vector<double> media_progressiva( vector<double> v ){

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
vector<double> media_progressiva_quad( vector<double> v ){

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
vector<double> errore( vector<double> v ){

	int N = v.size();
	vector<double> mp(N);
	mp = media_progressiva(v);	//vettore con le medie progressive
	vector<double> mpq(N);
	mpq = media_progressiva_quad(v);  //vettore con le medie progressive dei quadrati
	vector<double> err(N);

	for(int i=0; i<N; i++){
		err[i] = sqrt( mpq[i]-pow(mp[i],2) );
		err[i] /= sqrt(i+1);
	}
	return err;
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
