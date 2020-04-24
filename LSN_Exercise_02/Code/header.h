#include <vector>
#include <cmath>
using namespace std;

#define a 1		//lunghezza di ogn step

//coordinate sferiche
double X(double r, double theta, double phi){		//coordinata x
	return r*sin(theta)*sin(phi);
}

double Y(double r, double theta, double phi){		//coordinata y
	return r*sin(theta)*cos(phi);
}

double Z(double r, double theta, double phi){		//coordinata z
	return r*cos(theta);
}
//incertezze statistiche
double error(double sum, double sum2, int N){
	double e = sqrt( (sum2-pow(sum,2))/(N-1) );
	return e/(2*sqrt(sum));	
}
//modulo di un vettore
double modulo(vector<double> v){
	int n=v.size();
	double s=0;
	for(int i=0; i<n; i++)
		s += v[i]*v[i];
	return sqrt(s);
}
