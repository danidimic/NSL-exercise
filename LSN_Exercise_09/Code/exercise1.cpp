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
	
	Input();

	vec path = GeneratePath();
	cout<<path<<endl;
	cout<<CheckPath(path)<<endl<<endl;

	vec pperm = PairPermutation(path);
	cout<<"Pair permutation"<<endl;
	cout<<CheckPath(pperm)<<endl;
	cout<<pperm<<endl<<endl;

	vec shift = Shift(path, 3);
	cout<<"Shift di 3"<<endl;
	cout<<CheckPath(shift)<<endl;
	cout<<shift<<endl<<endl;

	vec rev = Reverse(path, 2, 5);
	cout<<"Reverse"<<endl;
	cout<<CheckPath(rev)<<endl;
	cout<<rev<<endl<<endl;
	
	rnd.SaveSeed();
	return 0;
}


//Inizializzo le variabili
void Input(void){
	
	rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	ifstream ReadInput("input.dat");
	ReadInput >> ncities;
	ReadInput >> ndim;
	
	cities.resize(ncities, ndim);	//matrice con le coordinate spaziali delle città
	cities.imbue( [&]() {return rnd.Rannyu(-side, side);} );	
	INDECES = linspace<vec>(1, ncities-1, ncities-1);		

	ReadInput >> npop;
	population.resize(npop, ncities);
	for(int i=0; i<npop; i++){
	
	}
}

//Definisco la funzione costo
double CostFunction(vec path, int p){
	int ipath;
	double s = 0, distance;
	rowvec cstart, cstop;
	
	//Calcolo la lunghezza del percorso
	for(int i=0; i<ncities; i++){
		ipath  = path[i];
		cstart = cities.row(ipath);			//città iniziale
		cstop  = cities.row(Pbc(ipath+1));	//città finale
	
		distance = norm(cstop-cstart);
		s += pow(distance, p);
	}
	return s;
}

//Controllo che il percorso sia corretto
bool CheckPath(vec path){
	path = sort(path);
	bool control = true;
	for(int i=0; i<path.size(); i++)
		if( path[i]!=i ){
			control = false;
			break;
		}
	return control;
}

//Genero un possibile percorso
vec GeneratePath(){
	vec newpath(ncities);
	vec index = shuffle(INDECES);

	newpath[0] = 0; //fisso la prima città del percorso
	for(int i=1; i<ncities; i++)
		newpath[i] = index[i-1];
	return newpath;
}

//Pair permutation di due città
vec PairPermutation(vec path, int i1, int i2){
	path.swap_rows(i1, i2);		//scambio le due città
	return path;
}

//Pair permutation casuale di due città
vec PairPermutation(vec path){
	int i1, i2;

	i1 = (int) rnd.Rannyu(1, ncities-1);
	for(;;){
		i2 = (int) rnd.Rannyu(1, ncities-1);
		if(i1 != i2)	break;
	}
	return PairPermutation(path, i1, i2);
}

//Shift delle città di m posti nel percorso
vec Shift(vec path, int m){
	if(m >= ncities-1)	
		return path;

	else{
		path.shed_row(0);
		path = shift(path, m);
		path.insert_rows(0, 1);
		return path;
	}
}

//Inversione dell'ordine delle città in un dato intervallo
vec Reverse(vec path, int start, int end){
	vec newpath(end-start);
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




