/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iomanip>
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "statistics.h"
#include "MolDyn_NVE.h"
using namespace std;

int main(){ 
	Input();             //Inizialization
	//int nconf = 1;
	iblock = 0;
	int nmeas = 0;

	for(int istep=1; istep <= nstep; ++istep){
		Move();           //Move particles with Verlet algorithm
		if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
		
		if(istep%10 == 0){
			Measure();     			//Properties measurement
			nmeas += 1;

			//ConfXYZ(nconf);		//Write actual configuration in XYZ format
			//nconf += 1;

			if(blocking==true) Accumulate(nmeas);	//somme su ciascun blocco
		}
	}
	
	if(blocking==true) Average();			//Write average values

  ConfFinal();	//Write final configuration to restart
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput, ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl << endl;
	dl = 0.5*box/nbins;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
	cout << "Tail correction for potential energy = " << vtail << endl;
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
	cout << "Tail correction for pressure = " << ptail << endl << endl;

	ReadInput >> restart;
	//ReadConf.open("old.0");
	//restart = ReadConf.is_open();

	ReadInput >> blocking;
	ReadInput >> nblock;
	nvalues = nstep/(10*nblock);
	ReadInput >> instant;

	if(blocking==true){
		ave_epot.resize(nblock);
		ave_ekin.resize(nblock);
		ave_etot.resize(nblock);
		ave_temp.resize(nblock);
		ave_pres.resize(nblock);
		ave_gdir.resize(nblock*nbins);
	}

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

	//Prepare system from only one configuration
	if( restart==false ){
		//Read initial configuration
		cout << "Read initial configuration from file config.0 " << endl << endl;
		ReadConf.open("config.0");
		for (int i=0; i<npart; ++i){
		  ReadConf >> x[i] >> y[i] >> z[i];
		  x[i] = x[i] * box;
		  y[i] = y[i] * box;
		  z[i] = z[i] * box;
		}
		ReadConf.close();

		//Prepare initial velocities
		 cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
		 double sumv[3] = {0.0, 0.0, 0.0};
		 for (int i=0; i<npart; ++i){
		   vx[i] = rand()/double(RAND_MAX) - 0.5;
		   vy[i] = rand()/double(RAND_MAX) - 0.5;
		   vz[i] = rand()/double(RAND_MAX) - 0.5;

		   sumv[0] += vx[i];
		   sumv[1] += vy[i];
		   sumv[2] += vz[i];
		 }
		 for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
		 double sumv2 = 0.0, fs;
		 for (int i=0; i<npart; ++i){
		   vx[i] = vx[i] - sumv[0];
		   vy[i] = vy[i] - sumv[1];
		   vz[i] = vz[i] - sumv[2];

		   sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		 }
		 sumv2 /= (double)npart;

		 fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
		 for (int i=0; i<npart; ++i){
		   vx[i] *= fs;
		   vy[i] *= fs;
		   vz[i] *= fs;

		   xold[i] = Pbc(x[i] - vx[i] * delta);
		   yold[i] = Pbc(y[i] - vy[i] * delta);
		   zold[i] = Pbc(z[i] - vz[i] * delta);
		 }
	}

	//Prepare the system from old configuratons
	if( restart==true ){
		cout << "Read initial configuration from file old.0 " << endl << endl;
		ReadConf.open("old.0");
		for (int i=0; i<npart; ++i){
		  ReadConf >> x[i] >> y[i] >> z[i];
		  x[i] = x[i] * box;
		  y[i] = y[i] * box;
		  z[i] = z[i] * box;
		}
		ReadConf.close();

		cout << "Read previous configuration from file old.final " << endl << endl;
		ReadConf.open("old.final");
		for (int i=0; i<npart; ++i){
		  ReadConf >> xold[i] >> yold[i] >> zold[i];
		  xold[i] = xold[i] * box;
		  yold[i] = yold[i] * box;
		  zold[i] = zold[i] * box;
		}
		ReadConf.close();
		
		Move();

		double fs, sumv2=0.0;
		for (int i=0; i<npart; ++i){
			vx[i] = (x[i]-xold[i])/(delta/2);
			vy[i] = (y[i]-yold[i])/(delta/2);
			vz[i] = (z[i]-zold[i])/(delta/2);

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;

		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = Pbc(x[i] - vx[i] * delta);
			yold[i] = Pbc(y[i] - vy[i] * delta);
			zold[i] = Pbc(z[i] - vz[i] * delta);
		}
	}

   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  double v, p, t, vij, pij;
  double dx, dy, dz, dr;
	double min, max, deltaV;

	//reset the hystogram of g(r)
	for (int k=0; k<nbins; ++k) gdir[k]=0.0;

  v = 0.0; //reset observables
  t = 0.0;
	p = 0.0;

//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){

			dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
			dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
			dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			//update of the histogram of g(r)
			for (int k=0; k<nbins; ++k){
					min = k*dl;
					max =(k+1)*dl;
					
					if( dr>min && dr<max )	gdir[k] += 2;
			}

			if(dr < rcut){
				//Potential energy
				vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
				v += vij;
				//Pressure
				pij = 48.0*( 1.0/pow(dr,12) - 0.5/pow(dr,6) );
				p += pij;
			}
		}          
	}

	//Normalizzo l'istogramma di g(r)
	for (int k=0; k<nbins; ++k){
		min = k*dl;
		max =(k+1)*dl;

		deltaV = 4/3*M_PI * ( pow(max,3)-pow(min,3) );
		gdir[k] /= (rho*npart*deltaV);
	}

	//Kinetic energy
	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
	stima_pot = v/(double)npart; //Potential energy per particle
	stima_kin = t/(double)npart; //Kinetic energy per particle
	stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
	stima_etot = (t+v)/(double)npart; //Total energy per particle
	stima_pres = rho*temp + p/(3*vol);	//Pressure

  ofstream Epot, Ekin, Etot, Temp, Pres;
	
	if(instant==true){
		Epot.open("results/epot.out",ios::app);
		Ekin.open("results/ekin.out",ios::app);
		Temp.open("results/temp.out",ios::app);
		Etot.open("results/etot.out",ios::app);
		Pres.open("results/pres.out",ios::app);

		Epot << stima_pot  << endl;
		Ekin << stima_kin  << endl;
		Temp << stima_temp << endl;
		Etot << stima_etot << endl;
		Pres << stima_pres << endl;

		Epot.close();
		Ekin.close();
		Temp.close();
		Etot.close();
		Pres.close();
	}

    return;
}


void Accumulate(int nmeas){		//accumulo le misure relative ai blocchi
	ave_epot.at(iblock) += stima_pot;
	ave_ekin.at(iblock) += stima_kin;
	ave_temp.at(iblock) += stima_temp;
	ave_etot.at(iblock) += stima_etot;
	ave_pres.at(iblock) += stima_pres;
	for(int k=0; k<nbins; k++)	ave_gdir[iblock*nbins+k] += gdir[k];

	if(nmeas%nvalues==0){		//valori medi su ciascun blocco
		ave_epot.at(iblock) /= nvalues;
		ave_ekin.at(iblock) /= nvalues;
		ave_temp.at(iblock) /= nvalues;
		ave_etot.at(iblock) /= nvalues;
		ave_pres.at(iblock) /= nvalues;

		//Tail corrections
		ave_epot.at(iblock) += vtail;					 		//Potential energy per particle
		ave_pres.at(iblock) += (npart*ptail)/vol;	//Pressure

		for(int k=0; k<nbins; k++)	ave_gdir[iblock*nbins+k] /= nvalues;
		iblock++;
	}
}

void Average(void){
	double r, ag;
	const int wd = 12;
	vector<double> g(nbins), gerr(nbins), s(nbins), s2(nbins);
	ofstream Epot, Ekin, Etot, Temp, Pres, Gdir;

	Epot.open("results/ave_epot.out",ios::app);
	Ekin.open("results/ave_ekin.out",ios::app);
	Temp.open("results/ave_temp.out",ios::app);
	Etot.open("results/ave_etot.out",ios::app);
	Pres.open("results/ave_pres.out",ios::app);
	Gdir.open("results/MD_Gave.out",ios::app);

	vector<double> cum_epot(nblock), cum_ekin(nblock), cum_etot(nblock), cum_temp(nblock), cum_pres(nblock);
	vector<double> err_epot(nblock), err_ekin(nblock), err_etot(nblock), err_temp(nblock), err_pres(nblock);

	//medie progressive
	cum_epot = cumulative_average( ave_epot );
	cum_ekin = cumulative_average( ave_ekin );
	cum_etot = cumulative_average( ave_etot );
	cum_temp = cumulative_average( ave_temp );
	cum_pres = cumulative_average( ave_pres );
	//incertezze statistiche
	err_epot = cumulative_error( ave_epot );
	err_ekin = cumulative_error( ave_ekin );
	err_etot = cumulative_error( ave_etot );
	err_temp = cumulative_error( ave_temp );
	err_pres = cumulative_error( ave_pres );

	for(int i=0; i<nblock; i++){
		Epot << i << setprecision(9) << "  " << cum_epot[i] << "  " << err_epot[i] << endl;
		Ekin << i << setprecision(9) << "  " << cum_ekin[i] << "  " << err_ekin[i] << endl;
		Temp << i << setprecision(9) << "  " << cum_temp[i] << "  " << err_temp[i] << endl;
		Etot << i << setprecision(9) << "  " << cum_etot[i] << "  " << err_etot[i] << endl;
		Pres << i << setprecision(9) << "  " << cum_pres[i] << "  " << err_pres[i] << endl;
	}

	//calcolo g(r) e l'incertezza utilizzando tutti i blocchi
	for(int iblock=0; iblock<nblock; iblock++){
		for(int k=0; k<nbins; k++){
			ag = ave_gdir[iblock*nbins+k];
			s[k]  += ag;
			s2[k] += ag*ag;
		}
	}
	for(int k=0; k<nbins; k++){
		r = (k+0.5)*dl;
		g[k] = s[k] / (double)nblock;
		gerr[k] = Error(s[k], s2[k], nblock);

		Gdir << setw(wd) << r << setw(wd) << g[k] << setw(wd) << gerr[k] << endl;
	}

	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();
	Pres.close();
	Gdir.close();
	return;
}

double Error(double sum, double sum2, int iblk){
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file old.0 " << endl << endl;
  WriteConf.open("old.0");
  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  cout << "Print final configuration to file old.final " << endl << endl;
  WriteConf.open("old.final");
  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();

  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
