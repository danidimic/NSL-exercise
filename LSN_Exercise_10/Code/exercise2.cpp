#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "mpi.h"
#include "random.h"
#include "statistics.h"
#include "genetic_algorithm.h"
using namespace std;


int main(int argc, char *argv[]){
	int size, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat;

	int ibest, n=1;
	double L[4], aveL, l;
	ofstream Lenght, Avelenght, Cities, Path;

	Input(rank);

	Lenght.open("results/lenght_rank" + to_string(rank) + ".out");
	Avelenght.open("results/avelenght_" + to_string(rank) + ".out");

	for(int i=0; i<nstep; i++){

		Generation();		//Nuova generazione di percorsi

		if(i%10==0){
			//lunghezza del miglior percorso di ciascun nodo
			ibest = fitness.index_max();
			Lenght<<i<<"  "<<CostFunction(population.row(ibest))<<endl;
			//lunghezza media della migliore metà popolazione di ciascun nodo
			aveL = BestHalf();
			Avelenght<<i<<"  "<<aveL<<endl;
		}

		if(i%nmigr==0){	//Scambio dei percorsi migliori ogni nmigr generazioni
			if(rank==0) cout<<endl<<"Migrazione numero = "<<n<<" / "<<nstep/nmigr<<endl;
			n++;

			for(int j=0; j<size; j++){
				//condivido il percorso migliore
				ShareBestPath(size, j);
				//inserisco nella popolazione il percorso migliore
				if(rank == change.back()) GetPath();
				change.resize(0);
			}
		}
	}
	Lenght.close();
	Avelenght.close();

	ibest = fitness.index_max();
	l = CostFunction( population.row(ibest) );
	MPI_Gather(&l, 1, MPI_DOUBLE, L, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(rank==0){
		ibest = 0;
		l = L[0];
		for(int i=0; i<size; i++) 
			if(L[i]<l){
				l = L[i];
				ibest = i;
			}
	}
	MPI_Bcast(&ibest, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

	if(rank==ibest){
		Path.open("results/best_path.out");
		ibest = fitness.index_max();
		Path<<population.row(ibest)<<endl;
		Path.close();
	}

	rnd.SaveSeed();
	MPI_Finalize();
	return 0;

}


//Inizializzo le variabili
void Input(int rank){
	rnd = RandomGenerator(rank);
	ifstream ReadInput("input.dat");

	ReadInput >> ncities;
	ReadInput >> ndim;
	ReadInput >> side;
	INDECES = linspace<vec>(1, ncities-1, ncities-1);

	ReadInput >> npop;
	ReadInput >> elite;
	elsize = elite*npop;	//dimensioni dell'elite
	ReadInput >> nstep;
	ReadInput >> nmigr;

	if(rank == 0){
		cout << "Algoritmo genetico per la soluzione del 'Travelling Salesman Problem'"<<endl<<endl;

		cout<<"Individui nella popolazione = "<<npop<<endl;
		cout<<"Dimensioni dell'elite = "<<elsize<<endl;
		cout<<"Iterazioni dell'algoritmo genetico = "<<nstep<<endl<<endl;

		cout<<"Migrazioni tra individui migliori ogni "<<nmigr<<" generazioni"<<endl;
	}
	CreateCities();			//inizializzo l'insieme di città

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

//Inizializzo un insieme di città da file
void CreateCities(){
	rowvec city(2);
	ifstream Cities("cities.in");
	for(int i=0; i<ncities; i++){
		Cities >> city[0];
		Cities >> city[1];
		cities.insert_rows(i, city);
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

//Lunghezza media della metà migliore della popolazione
double BestHalf(){
	int ibest;
	vec copyfit = fitness;
	for(int j=0; j<(int)npop/2; j++){
		ibest = copyfit.index_max();
		bestpop[j] = CostFunction(population.row(ibest));
		copyfit[ibest] = 0;
	}
	return mean(bestpop);
}

//Scambio dei percorsi migliori
void ShareBestPath(int size, int root){

	int index, r;
	index = fitness.index_max();
	rowvec path = population.row(index);
	//metto il percorso in uno std vector
	for(int i=0; i<ncities; i++)
		change.push_back( path[i] );
	//seleziono il processo ricevente
	for(;;){
		r = rnd.Rannyu(0, size);
		if(r!=root) break;
	}
	change.push_back(r);

	//Passo il vettore change a tutti i processi 
	MPI_Bcast(&change.front(), change.size(), MPI_INTEGER, root, MPI_COMM_WORLD);
}

//Inserisco change nella popolazione
void GetPath(){
	rowvec path(ncities);
	for(int i=0; i<ncities; i++)
		path[i] = change[i];

	int index = Select();
	population.row(index) = path;
}

//Periodic boundary condition
int Pbc(int index){
	
	if(index >= ncities)
		return index-ncities;
	else
		return index;
}
