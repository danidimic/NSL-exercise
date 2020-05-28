#include <armadillo>
#include <cstdlib>
#include "random.h"
using namespace std;
using namespace arma;

Random rnd;

//Genetic algorithm
mat cities;
vec INDECES;
rowvec currentpath, proposedpath;
double pmpp, pmsh, pmrev, side;
int power = 1, ncities, ndim, nstep;

//Simulated Annealing
vec temp;
double beta, Tmin = 0.01, Tmax = 5.; 

vector<double> lenght;

//Funzioni
void Input(void);
void CreateCities(void);
void UpdateProbabilities(void);
bool CheckPath(rowvec);
int Pbc(int);
double CostFunction(rowvec);
double min(double, double);
rowvec GeneratePath();
rowvec Mutation(rowvec);
rowvec Shift(rowvec);
rowvec Reverse(rowvec);
rowvec PairPermutation(rowvec);

