#include <armadillo>
#include <cstdlib>
#include "random.h"
using namespace std;
using namespace arma;

Random rnd;

//Genetic algorithm
rowvec bestpath;
vec INDECES, fitness, bestpop;
mat cities, population, newgeneration;
double pcross, pmpp, pmsh, pmrev, elite, side;
int power = 1, ncities, ndim, npop, nstep, elsize;

vector<double> lenght, avelenght;

//Funzioni
void Input(void);
void FitnessFunc(void);
void Generation(void);
void CreateCities(void);
bool CheckPath(rowvec);
int Pbc(int);
int Select(void);
double CostFunction(rowvec);
double BestHalf();
mat Crossover();
rowvec GeneratePath();
rowvec Mutation(rowvec);
rowvec Shift(rowvec);
rowvec Reverse(rowvec);
rowvec PairPermutation(rowvec);




