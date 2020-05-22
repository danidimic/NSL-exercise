#include <armadillo>
#include <cstdlib>
#include "random.h"
using namespace std;
using namespace arma;

Random rnd;

//Genetic algorithm
int power;
rowvec bestpath;
vec INDECES, fitness;
mat cities, population;

double side;
int ncities, ndim, npop, nstep;

vector<double> lenght, avelenght;

//Funzioni
void Input(void);
void FitnessFunc(int p);
void Crossover(void);
void CreateCities(void);
void Mutation(rowvec);
bool CheckPath(rowvec);
double CostFunction(rowvec, int n);
int Pbc(int);
int Select(void);
rowvec GeneratePath();
rowvec Shift(rowvec, int m);
rowvec PairPermutation(rowvec);
rowvec PairPermutation(rowvec, int, int);
rowvec Reverse(rowvec, int, int);




