#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "header.h"
#include "random.h"
#include "statistics.h"
using namespace std;

int main (int argc, char *argv[]){
	
	unsigned int ibest;
	double pc, L;
	ofstream Cities, Path, Lenght;

	Input();

	ibest = fitness.index_min();
	L = CostFunction(population.row(ibest), power);
	lenght.push_back(L);

	for(int i=0; i<nstep; i++){
		pc = rnd.Rannyu();		
		if ( pc>0.5 ){
			Crossover();
			ibest = fitness.index_min();
			L = CostFunction(population.row(ibest), power);
			lenght.push_back(L);
		}
	}


	Cities.open("results/cities.out");
	for(int i=0; i<ncities; i++)
		Cities<<cities.row(i)<<endl;
	Cities.close();

	Path.open("results/best_path.out"); 
	for(int i=0; i<ncities; i++)
		Path<<population.row(ibest)<<endl;
	Path.close();

	Lenght.open("results/lenght.out");
	for(int i=0; i<(int)lenght.size(); i++)
		Lenght<<i<<"  "<<lenght[i]<<endl;

	rnd.SaveSeed();
	return 0;
}

//Inizializzo le variabili
void Input(void){
	
	rnd = RandomGenerator();
	ifstream ReadInput("input.dat");

	cout << "Algoritmo genetico per la soluzione del 'Traveling Salesman Problem'"<<endl<<endl;
	ReadInput >> ncities;
	cout<<"Numero di città = "<<ncities<<endl;

	ReadInput >> ndim;
	ReadInput >> side;
	CreateCities();		//inizializzo l'insieme di città
	INDECES = linspace<vec>(1, ncities-1, ncities-1);

	ReadInput >> npop;
	cout<<"Popolazione di "<<npop<<" individui"<<endl;
	ReadInput >> nstep;
	cout<<"Iterazioni dell'algoritmo genetico = "<<nstep<<endl;
	ReadInput >> power;

	//Creo una popolazione di npop individui
	rowvec path(ncities);
	for(int i=0; i<npop; i++){

		for(;;){		//Genero la popolazione
			path = GeneratePath();
			if( CheckPath(path)==1 ) break;
		}
		population.insert_rows(i, path);
	}
	//Calcolo la funzione fitness
	fitness.resize(npop);
	FitnessFunc(power);
}

//Inizializzo casualmente un insieme di città
void CreateCities(){

	rowvec city(2);

	//Città lungo la circonferenza
	if(ndim==1){
		double theta, r = side;
		cout<<"Città disposte lungo una circonferenza di raggio = "<<side<<endl<<endl;

		for(int i=0; i<ncities; i++){
			theta = rnd.Rannyu(0, 2*M_PI);
			city[0] = r*cos(theta);
			city[1] = r*sin(theta);
			cities.insert_rows(i, city);
		}
	}
	//Città all'interno del quadrato
	if(ndim==2){
		cout<<"Città disposte all'interno di un quadrato di lato = "<<side<<endl<<endl;
		cities.resize(ncities, 2);
		cities.imbue( [&]() {return rnd.Rannyu(0, side);} );
	}
}

//Definisco la funzione costo
double CostFunction(rowvec path, int p){
	int istart, istop;
	double s = 0, distance;
	rowvec cstart, cstop;
	
	//Calcolo la lunghezza del percorso
	for(int i=0; i<ncities; i++){
		istart = path[i];
		istop  = Pbc(path[i+1]);

		cstart = cities.row(istart);		//città iniziale
		cstop  = cities.row(istop);			//città finale
		distance = norm(cstop-cstart);
		s += pow(distance, p);
	}
	return s;
}

//Calcolo la funzione Fitness
void FitnessFunc(int p){
	rowvec individual(ncities);
	for(int i=0; i<npop; i++){
		individual = population.row(i);
		fitness[i] = CostFunction(individual, p);
	}
	fitness /= (double)fitness.max();
	fitness = 1-fitness;
	fitness /= (double)fitness.max();
}

//Seleziona un individuo
int Select(void){
	int r;
	double p;
	for(;;){
		r = rnd.Rannyu(0, npop);
		p = rnd.Rannyu(0.5, 1.);
		if( p<fitness[r] )	break;
	}
	return r;
}

//Crossover per nuovi individui
void Crossover(){

	int ifather, imother;
	rowvec father, mother;

	ifather = Select();
	for(;;){
		imother = Select();
		if( imother!=ifather )
			break;
	}
	father = population.row(ifather);
	mother = population.row(imother);
	
	int crosspoint;
	for(;;){
		crosspoint = rnd.Gauss( (int) ncities/2, (double) ncities/4 );
		if(crosspoint>1 && crosspoint<ncities-1) break;
	}
	//Stacco i pezzi fino al crosspoint
	rowvec son1(ncities), son2(ncities);
	for(int i=0; i<crosspoint; i++){
		son1[i] = father[i];
		son2[i] = mother[i];
	}
	//Ottengo le parti restanti nell'ordine
	int i1 = crosspoint, i2 = crosspoint, f, m;
	for(int i=0; i<ncities; i++){
		f = father[i];
		m = mother[i];

		if( count(son1.begin(), son1.end(), m)==0 ){
			son1[i1] = m;
			i1++;
		}
		if( count(son2.begin(), son2.end(), f)==0 ){
			son2[i2] = f;
			i2++;
		}
	}

	//Modifico la popolazione
	int imin;
	if( CheckPath(son1)==1 && CheckPath(son2)==1 ){
		//Eventuali mutazioni
		Mutation(son1);
		Mutation(son2);

		//Aggiungo il primo figlio
		imin = fitness.index_min();
		population.row(imin) = son1;
		FitnessFunc(power);
		//Aggiungo il secondo figlio
		imin = fitness.index_min();
		population.row(imin) = son2;
		FitnessFunc(power);
	}
	else	Crossover();
}

//Inserisco eventuali mutazioni 
void Mutation(rowvec individual){

	int pm = rnd.Rannyu(0, 10);
	int shift, start, stop;
	switch(pm){
		case 0:
			individual = PairPermutation(individual);
			break;
		case 1:
			shift = rnd.Rannyu(1, ncities);
			individual = Shift( individual, shift );
			break;
		case 2:
			start = rnd.Rannyu(1, ncities);
			stop  = rnd.Rannyu(1, ncities);
			individual = Reverse( individual, start, stop );
			break;
		default:
			break;
	}
}

//Genero un possibile percorso
rowvec GeneratePath(){
	rowvec newpath(ncities);
	rowvec index = shuffle(INDECES).as_row();

	newpath[0] = 0; //fisso la prima città del percorso
	for(int i=1; i<ncities; i++)
		newpath[i] = index[i-1];
	return newpath;
}

//Controllo che il percorso sia corretto
bool CheckPath(rowvec path){
	path = sort(path);
	bool control = true;
	for(int i=0; i<(int)path.size(); i++)
		if( path[i]!=i ){
			control = false;
			break;
		}
	return control;
}

//Pair permutation di due città
rowvec PairPermutation(rowvec path, int i1, int i2){
	path.swap_cols(i1, i2);		//scambio le due città
	return path;
}

//Pair permutation casuale di due città
rowvec PairPermutation(rowvec path){
	int i1, i2;

	i1 = (int) rnd.Rannyu(1, ncities-1);
	for(;;){
		i2 = (int) rnd.Rannyu(1, ncities-1);
		if(i1 != i2)	break;
	}
	return PairPermutation(path, i1, i2);
}

//Shift delle città di m posti nel percorso
rowvec Shift(rowvec path, int m){
	if(m >= ncities-1)	
		return path;

	else{
		path.shed_col(0);
		path = shift(path, m);
		path.insert_cols(0, 1);
		return path;
	}
}

//Inversione dell'ordine delle città in un dato intervallo
rowvec Reverse(rowvec path, int start, int end){

	if(start>end) swap(start, end);

	rowvec newpath(end-start);
	for(int i=0; i<end-start; i++)
		newpath[i] = path[i+start];

	newpath = reverse(newpath);
	for(int i=0; i<end-start; i++)
		path[i+start] = newpath[i];

	return path;
} 

//Periodic boundary condition
int Pbc(int index){
	
	if(index >= ncities)
		return index-ncities;
	else
		return index;
}




