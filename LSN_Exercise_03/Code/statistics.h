#ifndef __Statistics__
#define __Statistics__

#include <vector>
#include <string>
using namespace std;

vector<double> cumulative_average( vector<double> );
vector<double> cumulative_average_quad( vector<double> );
vector<double> cumulative_error( vector<double> );
void data_blocking(vector<double>, string);
double chi_quadro( vector<double>, double );

#endif // __Statistics__
