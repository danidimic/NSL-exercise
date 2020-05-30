#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("../Files/seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  WriteSeed.close();
  return;
}

//Generazione di numeri con distribuzione Gaussiana
double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}
//Generazione di numeri con distribuzione Esponenziale
double Random:: Expo(double lambda){
	double x=0;
	double y;
	if(lambda != 0){
		y=Rannyu();
		x=-1/lambda*log(1-y);
	}
	return x;
}
//Generazione di numeri con distribuzione Lorentziana
double Random:: Lorentz(double mu, double gamma){
	double y=Rannyu();
	double x=mu+gamma*tan( M_PI*(y-0.5) );
	return x;
}
//Generazione di numeri uniformi nell'intervallo [min, max]
double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}
//Generazione di numeri uniformi nell'intervallo [0, 1]
double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}
//Generazioni di valori uniformi sulla sfera
double* Random:: Sphere(){
	double* s = new double[2];
	s[0] = acos( 2*Rannyu()-1 );		//angolo theta
	s[1] = Rannyu(0, 2*M_PI);			//angolo phi
	return s;	
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

Random RandomGenerator(){

	Random rnd;

	int seed[4];
 	int p1, p2;
  	ifstream Primes("../Files/Primes");
  	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
  	} else cerr << "PROBLEM: Unable to open Primes" << endl;
  	Primes.close();

  	ifstream input("../Files/seed.in");
  	string property;
  	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	return rnd;
}	

Random RandomGenerator(int rank){

	Random rnd;

	int seed[4];
 	int p1, p2;
  	ifstream Primes("../Files/Primes");
  	if (Primes.is_open()){
			for(int i=0; i<=rank; i++)		Primes >> p1 >> p2 ;
  	} 
		else cerr << "PROBLEM: Unable to open Primes" << endl;
  	Primes.close();

  	ifstream input("../Files/seed.in");
  	string property;
  	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	return rnd;
}	
