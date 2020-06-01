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

//Parallel Program
int nmigr;
vector<int> change;


//Funzioni
void Input(int);
void FitnessFunc(void);
void Generation(void);
void CreateCities(void);
void ExchangeBestPath(int, int);
int Pbc(int);
int Select(void);
bool CheckPath(rowvec);
double CostFunction(rowvec);
double BestHalf();
mat Crossover();
rowvec GeneratePath();
rowvec Mutation(rowvec);
rowvec Shift(rowvec);
rowvec Reverse(rowvec);
rowvec PairPermutation(rowvec);




