#include <armadillo>
#include <cstdlib>
#include "random.h"
using namespace std;
using namespace arma;

Random rnd;

//Genetic algorithm
int ncities, ndim, npop;
double side = 15;
vec INDECES;
mat cities, population;

//Funzioni
void Input(void);
bool CheckPath(vec);
double CostFunction(vec, int n);
int Pbc(int);
vec GeneratePath();
vec Shift(vec, int m);
vec PairPermutation(vec);
vec PairPermutation(vec, int, int);
vec Reverse(vec, int, int);




