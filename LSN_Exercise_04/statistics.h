#ifndef __Statistics__
#define __Statistics__

#include <vector>
using namespace std;

double modulo(vector<double> );
double media(vector<double> );
double devstd(vector<double>);
double devstd_media(vector<double>);
vector<double> media_progressiva( vector<double> );
vector<double> media_progressiva_quad( vector<double> );
vector<double> errore( vector<double> );
double chi_quadro( vector<double>, double );

#endif // __Statistics__
