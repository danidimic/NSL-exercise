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
int swapindex[4], nmigr;
vector<int> change;

//Funzioni
void Input(int);
void FitnessFunc(void);
void Generation(void);
void RandomExchange(void);
void CreateCities(int);
void ExchangeBestPath(MPI_Status, int);
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




