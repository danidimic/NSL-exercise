#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "header.h"
#include "random.h"
#include "statistics.h"
using namespace std;


int main(int argc, char *argv[]){

	int ibest;
	double L, aveL;
	vec copyfit;
	ofstream Cities, Path, Lenght;

	Input();

	for(int i=0; i<nstep; i++){
		
		if(i%50==0){
			cout<<"Iterazione "<<i<<endl;
			ibest = fitness.index_max();
			cout<<"L min = "<<CostFunction( population.row(ibest) )<<endl;
			cout<<sum(fitness)<<endl<<endl;
		}

		Generation();

		copyfit = fitness;
		for(int j=0; j<(int)npop/2; j++){
			ibest = copyfit.index_max();
			bestpop[j] = CostFunction(population.row(ibest));
			copyfit[ibest] = 0;
		}

		//Lunghezza del percordo migliore
		ibest = fitness.index_max();
		L = CostFunction(population.row(ibest));
		lenght.push_back(L);
		//Lunghezza mediata sulla popolazione migliore
		aveL = mean(bestpop);
		avelenght.push_back(aveL);
	}

	Path.open("results/best_path.out");
	ibest = fitness.index_max();
	Path<<population.row(ibest)<<endl;
	Path.close();

	Cities.open("results/cities.out");
	for(int i=0; i<ncities; i++)
		Cities<<cities.row(i)<<endl;
	Cities.close();

	Lenght.open("results/lenght.out");
	for(int i=0; i<(int)lenght.size(); i++)
		Lenght<<i<<"  "<<lenght[i]<<endl;
	Lenght.close();

	Lenght.open("results/avelenght.out");
	for(int i=0; i<(int)avelenght.size(); i++)
		Lenght<<i<<"  "<<avelenght[i]<<endl;
	Lenght.close();

	rnd.SaveSeed();
	return 0;

}


//Inizializzo le variabili
void Input(void){
	rnd = RandomGenerator();
	ifstream ReadInput("input.dat");

	cout << "Algoritmo genetico per la soluzione del 'Travelling Salesman Problem'"<<endl<<endl;
	ReadInput >> ncities;
	cout<<"Numero di città = "<<ncities<<endl;

	ReadInput >> ndim;
	ReadInput >> side;
	CreateCities();		//inizializzo l'insieme di città
	INDECES = linspace<vec>(1, ncities-1, ncities-1);

	ReadInput >> npop;
	cout<<"Individui nella popolazione = "<<npop<<endl;
	elsize = 0.1*npop;	//dimensioni elite
	cout<<"Dimensioni dell'elite = "<<elsize<<endl;
	ReadInput >> nstep;
	cout<<"Iterazioni dell'algoritmo genetico = "<<nstep<<endl<<endl;
	ReadInput >> power;

	//Probabilità crossover e mutazioni
	ReadInput >> pcross;
	ReadInput >> pmpp;
	ReadInput >> pmsh;
	ReadInput >> pmrev;

	ReadInput.close();

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
	FitnessFunc();
	bestpop.resize((int)npop/2);
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

//Definisco la funzione costo come lunghezza del percorso
double CostFunction(rowvec path){
	int istart, istop;
	double s = 0, distance;
	rowvec cstart, cstop;
	
	//Calcolo la lunghezza del percorso
	for(int i=0; i<ncities; i++){
		istart = path[i];
		istop  = path[Pbc(i+1)];

		cstart = cities.row(istart);		//città iniziale
		cstop  = cities.row(istop);			//città finale
		distance = norm(cstop-cstart);	//distanza tra città
		s += pow(distance, power);
	}
	return s;
}

//Calcolo la funzione Fitness
void FitnessFunc(){
	rowvec individual(ncities);
	for(int i=0; i<npop; i++){
		individual = population.row(i);
		fitness[i] = 1. / CostFunction(individual);
	}
	//fitness /= fitness.max();
	//fitness  = 1 - fitness;
	//fitness = normalise(fitness);
}

//Seleziona un individuo
int Select(void){
	int j;
	double F, alpha, isum = 0;
	vec fit = fitness;

	F = sum(fitness);
	alpha = rnd.Rannyu(0, F);
	do{
		j = fit.index_max();
		isum += fitness[j];
		fit[j] = 0;
	}while(isum<alpha && j<npop-1);
	return j;
}

//Ricavo l'elite della popolazione
void Elite(){
	int ibest;
	vec copyfit = fitness;
	newgeneration.resize(0,0);

	for(int i=0; i<elsize; i++){
		ibest = copyfit.index_max();
		newgeneration.insert_rows(i, population.row(ibest));
		copyfit[ibest] = 0;
	}
}

//Genero una nuova popolazione
void Generation(){
	double pc;
	mat sons;
	rowvec son1, son2;

	Elite();

	for(int j=elsize; j<npop; j+=2){	//Creo una nuova generazione

		pc = rnd.Rannyu();
		if( pc<pcross ){
			sons = Crossover();
			son1 = sons.row(0);
			son2 = sons.row(1);
			//Inserisco eventuali mutazioni
			son1 = Mutation(son1);
			son2 = Mutation(son2);
		}else{
			son1 = population.row( Select() );
			son2 = population.row( Select() );
		}

		newgeneration.insert_rows(j, son1);
		newgeneration.insert_rows(j+1, son2);
	}
	population.resize(0,0);
	population = newgeneration;
	FitnessFunc();
	//cout<<"qui ci sono"<<endl;
}

//Genero due figli dalla popolazione
mat Crossover(){
	rowvec father, mother;
	int ifather = Select(), imother = Select();
	father = population.row(ifather);
	mother = population.row(imother);

	//Stacco i pezzi fino al crosspoint
	int crosspoint = rnd.Rannyu(1, ncities-1), l = ncities-crosspoint;
	rowvec son1(ncities, fill::zeros), son2(ncities, fill::zeros);
	for(int i=0; i<crosspoint; i++){
		son1[i] = father[i];
		son2[i] = mother[i];
	}
	rowvec restF(l), restM(l);
	for(int i=0; i<l; i++){
		restF[i] = father[i+crosspoint];
		restM[i] = mother[i+crosspoint];
	}

	int f, m;
	rowvec restSon1(l, fill::zeros), restSon2(l, fill::zeros);
	//Completo figlio 1
	int i=0;
	for(int j=0; j<ncities; j++){
		m = mother[j];
		if( count(restF.begin(), restF.end(), m)==1 ){
			son1[crosspoint+i] = m;
			i++;
		}
	}
	//Completo figlio 2
	i = 0;
	for(int j=0; j<ncities; j++){
		f = father[j];
		if( count(restM.begin(), restM.end(), f)==1 ){
			son2[crosspoint+i] = f;
			i++;
		}
	}

	mat sons;
	sons.insert_rows(0, son1);
	sons.insert_rows(1, son2);
	return sons;
}

//Inserisco eventuali mutazioni 
rowvec Mutation(rowvec individual){
	double pm = rnd.Rannyu();

	if(pm<pmpp) individual = PairPermutation(individual);

	else if( pm>pmpp and pm<pmpp+pmsh) individual = Shift(individual);
	
	else if( pm>pmpp+pmsh and pm<pmpp+pmsh+pmrev) individual = Reverse(individual);

	return individual;
}

//Shift delle città del percorso
rowvec Shift(rowvec path){
	//shift di 2/3/4 posti nel percorso
	int m = rnd.Rannyu(2, 5);

	path.shed_col(0);
	path = shift(path, m);
	path.insert_cols(0, 1);
	return path;
}

//Pair permutation casuale di due città
rowvec PairPermutation(rowvec path){
	int i1, i2;
	i1 = (int) rnd.Rannyu(1, ncities-1);
	i2 = (int) rnd.Rannyu(1, ncities-1);

	path.swap_cols(i1, i2); 	//scambio le due città nel percorso
	return path;
}

//Inversione dell'ordine delle città in un dato intervallo
rowvec Reverse(rowvec path){
	int start, end;
	start = rnd.Rannyu(1, ncities-10);
	end = start + (int)rnd.Rannyu(2, 5);

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
