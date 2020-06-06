#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "statistics.h"
#include "simulated_annealing.h"
using namespace std;


int main(int argc, char *argv[]){

	rowvec bestpath;
	int naccept = 0, ntot = 0;
	double delta, alpha, r, L, bestL;
	ofstream Cities, Path, Lenght, Prob;

	Input();
	
	bestpath = currentpath;
	bestL = CostFunction(currentpath);

	for(int i=0; i<(int)temp.size(); i++){
		beta = 1. / temp[i]; 		//temperatura fittizia
		cout<<endl<<"Temperatura fittizia T = "<<temp[i]<<endl<<endl;

		for(int j=0; j<nstep; j++){

			proposedpath = Mutation(currentpath);

			r = rnd.Rannyu();
			delta = CostFunction(proposedpath) - CostFunction(currentpath);
			alpha = min( 1, exp(-beta*delta) );

			if(r<alpha){
				currentpath = proposedpath;
				naccept++;

				L = CostFunction(currentpath);
				if(L < bestL){
					bestL = L;
					bestpath = currentpath;
				}
			}

			if(j%100000==0){
				cout<<"Generazione = "<<ntot<<endl;
				cout<<"Lunghezza minima = "<<bestL<<endl;

				lenght.push_back(bestL);
			}

			ntot++;
		}
		UpdateProbabilities();
	}

	cout<<endl<<"Generazioni totali = "<<ntot<<endl;
	cout<<"Accettazione del metropolis = "<<(double) naccept/ntot<<endl;
	
	Cities.open("results/cities.out");
	for(int i=0; i<ncities; i++)
		Cities<<cities.row(i)<<endl;
	Cities.close();

	Lenght.open("results/lenght.out");
	for(int i=0; i<(int)lenght.size(); i++)
		Lenght<<i<<"  "<<lenght[i]<<endl;
	Lenght.close();

	Path.open("results/best_path.out");
	Path<<bestpath<<endl;
	Path.close();

	Prob.open("results/probabilities.out");
	for(int i=0; i<(int)temp.size(); i++)
		Prob<<temp[i]<<"  "<<PMpp[i]<<"  "<<PMsh[i]<<"  "<<PMrev[i]<<"  "<<PMmul[i]<<endl;
	Prob.close();

	rnd.SaveSeed();
	return 0;

}


//Inizializzo le variabili
void Input(void){
	rnd = RandomGenerator();
	ifstream ReadInput("input.dat");

	cout << "Simulated Annealing per la soluzione del 'Travelling Salesman Problem'"<<endl<<endl;
	ReadInput >> ncities;
	cout<<"Numero di città = "<<ncities<<endl;

	ReadInput >> ndim;
	ReadInput >> side;
	CreateCities();		//inizializzo l'insieme di città
	INDECES = linspace<vec>(1, ncities-1, ncities-1);

	ReadInput >> nstep;
	cout<<"Iterazioni dell'algoritmo genetico = "<<nstep<<endl;
	temp = linspace<vec>(Tmax, Tmin, 25);
	cout<<"Temperatura fittizia da "<<Tmax<<" a "<<Tmin<<endl<<endl;

	//Probabilità mutazioni
	ReadInput >> pmpp;
	ReadInput >> pmsh;
	ReadInput >> pmrev;
	ReadInput >> pmul;

	ReadInput.close();

	currentpath = GeneratePath();
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

//Inserisco eventuali mutazioni 
rowvec Mutation(rowvec individual){
	double pm = rnd.Rannyu();

	if(pm<pmpp) individual = PairPermutation(individual);

	else if( pm>pmpp and pm<pmpp+pmsh) individual = Shift(individual);
	
	else if( pm>pmpp+pmsh and pm<pmpp+pmsh+pmrev) individual = Reverse(individual);

	else if( pm>pmpp+pmsh+pmrev and pm<pmpp+pmsh+pmrev+pmul) individual = MultiPermutation(individual, (int)rnd.Rannyu(2,8) );

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

rowvec MultiPermutation(rowvec path, int nperm){
	int i1, i2;
	for(int i=0; i<nperm; i++){
		i1 = (int) rnd.Rannyu(1, ncities-1);
		i2 = (int) rnd.Rannyu(1, ncities-1);
		path.swap_cols(i1, i2);
	}
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

//Aggiorna probabilità di mutazione
void UpdateProbabilities(void){

	PMpp.push_back(pmpp);
	PMsh.push_back(pmsh);
	PMrev.push_back(pmrev);
	PMmul.push_back(pmul);

	pmpp *= 0.98;
 	pmsh *= 0.97;
	pmrev *= 0.92;
	pmul *= 0.95;
}


//Minimo tra due valori
double min(double a, double b){
	double Min;
	if(a<b) Min = a;
	else	Min = b;
	return Min;
}
