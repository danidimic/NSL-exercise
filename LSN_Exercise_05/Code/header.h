#include <armadillo>
using namespace std;
using namespace arma;


#define a0 1	//raggio di Bohr

//minimo tra due valori
double min(double a, double b){
	double Min;
	if(a<b) Min = a;
	else	Min = b;
	return Min;
}

//coordinata theta a partire da x,y,z
double Theta(double x, double y){
	double t=0;
	if( x==0 and y>0 ) t = M_PI/2;
	else if ( x==0 and y<0 ) t = 3*M_PI/2;
	else if ( x>0 and y>=0 ) t = atan(y/x);
	else if ( x*y<0 ) t = atan(y/x) + 2*M_PI;
	else if ( x<0 and y<=0 ) t = atan(y/x) + M_PI;

	return t;
}
//orbitale atomico n=1, l=0, m=0
double Psi100(double x, double y, double z){
	double r = sqrt( x*x + y*y + z*z );
	double C = pow(a0, -3.0/2.0)/sqrt(M_PI);
	return C*exp(-r/a0);
}
//orbitale atomico n=2, l=1, m=0
double Psi210(double x, double y, double z){
	double r = sqrt( x*x + y*y + z*z );
	double theta = Theta(x,y);

	double C = pow(a0, -5.0/2.0)/8*sqrt(2/M_PI);
	return C*r*exp( -r/(2*a0) )*cos(theta);
}
//densità di probabilità
double pdf( double (*func) (double, double, double), vec v ){
	double p, x, y, z;
	x = v[0];
	y = v[1];
	z = v[2];

	p = abs( func(x, y, z) );
	return pow(p, 2);
}

void write_points(vector<double>, string);
