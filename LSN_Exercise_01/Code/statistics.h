#ifndef __Statistics__
#define __Statistics__

#include <vector>
using namespace std;

vector<double> cumulative_average( vector<double> );
vector<double> cumulative_average_quad( vector<double> );
vector<double> cumulative_error( vector<double> );
double chi_quadro( vector<double>, double );

#endif // __Statistics__
