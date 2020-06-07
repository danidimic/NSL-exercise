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
double pmpp, pmsh, pmrev, pmul, side;
int power = 1, ncities, ndim, nstep;

//Simulated Annealing
vec temp;
int ntemp = 500;
double beta, Tmin = 0.01, Tmax = 5.;

vector<double> lenght, PMpp, PMsh, PMrev, PMmul;

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
rowvec MultiPermutation(rowvec, int);

