/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){

	Input(); //Inizialization
	for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation

		Reset(iblk);   //Reset block averages
		for(int istep=1; istep <= nstep; ++istep){

			Move(metro);
			Measure();
			Accumulate(); //Update block averages
		}
		Averages(iblk);   //Print results for current block
	}
	ConfFinal(); //Write final configuration
  return 0;
}


void Input(void){
  ifstream ReadInput;

	cout << "Classic 1D Ising model             " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
	cout << "Nearest neighbour interaction      " << endl << endl;
	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
	cout << "The program uses k_B=1 and mu_B=1 units " << endl;

	//Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
  
	//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;
  ReadInput >> nstep;
	ReadInput >> restart;
	ReadInput >> vsTemp;		//misure in funzione della temperatura

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

	//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

	//initial configuration from file conf.final
	if( restart==1 ){
		ReadInput.open("config.final");
		for(int i=0; i<nspin; i++)
			ReadInput >> s[i];
	}
	//random initial configuration
	else{
		for (int i=0; i<nspin; ++i){
			if(rnd.Rannyu() >= 0.5) s[i] = 1;
			else s[i] = -1;
		}
	}

	//Evaluate energy etc. of the initial configuration
	Measure();
	//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro){
	int o;
	double A, delta_energy, p;

	for(int i=0; i<nspin; ++i){
		//Select randomly a particle
		o = (int)(rnd.Rannyu()*nspin);

		if(metro==1){ //Metropolis
			delta_energy = Boltzmann(s[o], o);
			A = exp( -beta * delta_energy);

			if( rnd.Rannyu()<A ){
				s[o] *= -1;
	   			accepted++;
			}
			attempted++;			
		}

		else{ //Gibbs sampling
			delta_energy = Gibbs(o);
			p = 1 / ( 1+exp(delta_energy) );

			if( rnd.Rannyu()<p ) 
				s[o] = 1;
			else
				s[o] = -1;
		}
	}
	return;
}

double Gibbs(int ip){
	double ene = -2 * beta * J *( s[Pbc(ip-1)] + s[Pbc(ip+1)] );
	return ene;
}

double Boltzmann(int sm, int ip){
	double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
	return ene;
}

void Measure(){
	double u = 0.0, u2 = 0.0, m = 0.0;

	//cycle over spins
	for (int i=0; i<nspin; ++i){
		u  += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		u2 += pow(u, 2);
		m  += s[i];
  }
	//immagazzino i valori istantanei
	walker[iu] = u;														//energia
	walker[ic] = pow(beta, 2) * (u2 - u*u);		//capacità termica
	walker[im] = m;														//magnetizzazione
	walker[ix] = beta * pow(m,2);							//susciettività
	return;
}


void Reset(int iblk){ //Reset block averages

	if(iblk == 1){
		for(int i=0; i<n_props; ++i){
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i){
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void){ //Update block averages
	for(int i=0; i<n_props; ++i){
		blk_av[i] = blk_av[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){ //Print results for current block
    
	ofstream Ene, Heat, Mag, Chi;
	const int wd=12;
   
	if (vsTemp==0){
		cout << "Block number " << iblk << endl;
		cout << "Acceptance rate " << accepted/attempted << endl;
		cout << "----------------------------" << endl << endl;
	}
    
	Ene.open("results/ave_ene.out", ios::app);
	Heat.open("results/ave_heat.out", ios::app);
	Mag.open("results/ave_mag.out", ios::app);
	Chi.open("results/ave_chi.out", ios::app);

	stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
	glob_av[iu]  += stima_u;
	glob_av2[iu] += stima_u*stima_u;
	err_u=Error(glob_av[iu],glob_av2[iu],iblk);
	Ene << setw(5) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(15) << err_u << endl;

	stima_c = blk_av[ic]/blk_norm/(double)nspin; //Heat
	glob_av[ic]  += stima_c;
	glob_av2[ic] += stima_c*stima_c;
	err_c=Error(glob_av[ic],glob_av2[ic],iblk);
	Heat << setw(5) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(15) << err_c << endl;

	stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m=Error(glob_av[im],glob_av2[im],iblk);
	Mag << setw(5) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(15) << err_m << endl;

	stima_x = blk_av[ix]/blk_norm/(double)nspin; //Susceptivity
	glob_av[ix]  += stima_x;
	glob_av2[ix] += stima_x*stima_x;
	err_x=Error(glob_av[ix],glob_av2[ix],iblk);
	Chi << setw(5) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(15) << err_x << endl;

	Ene.close();
	Heat.close();
	Mag.close();
	Chi.close();

	if( iblk==nblk && vsTemp==1 ){			//scrivi risultati per grafici vs temperatura 
		Ene.open("results/Ene.out", ios::app);
		Heat.open("results/Heat.out", ios::app);
		Mag.open("results/Mag.out", ios::app);
		Chi.open("results/Chi.out", ios::app);

		Ene  << temp << "  " << glob_av[iu]/(double)iblk << "  " << err_u << endl;
		Heat << temp << "  " << glob_av[ic]/(double)iblk << "  " << err_c << endl;
		Mag  << temp << "  " << glob_av[im]/(double)iblk << "  " << err_m << endl;
		Chi  << temp << "  " << glob_av[ix]/(double)iblk << "  " << err_x << endl;

		Ene.close();
		Heat.close();
		Mag.close();
		Chi.close();
	}
	return;
}


void ConfFinal(void){
	ofstream WriteConf;

	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");
	for (int i=0; i<nspin; ++i){
		WriteConf << s[i] << endl;
	}
	WriteConf.close();

	rnd.SaveSeed();
}


int Pbc(int i){  //Algorithm for periodic boundary conditions
	if(i >= nspin) i = i - nspin;
	else if(i < 0) i = i + nspin;
	return i;
}


double Error(double sum, double sum2, int iblk){
	return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}
 

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
