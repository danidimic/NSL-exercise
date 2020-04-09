#ifndef __Statistics__
#define __Statistics__

#include "random.h"
#include <vector>
using namespace std;

//funzioni generiche
double modulo(vector<double> );
//funzioni statistiche
double media(vector<double> );
double devstd(vector<double>);
double devstd_media(vector<double>);
double devstd_media(int, double, double);
//funzioni medie a blocchi
vector<double> media_progressiva( vector<double> );
vector<double> media_progressiva_quad( vector<double> );
vector<double> errore( vector<double> );
double chi_quadro( vector<double>, double );
//integrali Monte Carlo
double uniform_Monte_Carlo(double, double, int, double (*func)(double x), Random);

#endif // __Statistics__
