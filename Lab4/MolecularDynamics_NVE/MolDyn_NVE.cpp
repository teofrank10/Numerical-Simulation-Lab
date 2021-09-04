/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     
#include <iostream>     
#include <fstream>      
#include <cmath>        
#include "MolDyn_NVE.h"


using namespace std;


MolDyn_NVE :: MolDyn_NVE(){ }

MolDyn_NVE :: ~MolDyn_NVE(){ }

void MolDyn_NVE :: Input(void){ //Prepare the parameters for the simulation
  	ifstream ReadInput;//from the input files we obtain the parameters of our system. from config.0 we read the initial configuration, i.e the positions of the particles at time t. But we know that in order to use the Verlet's algorithm we need also the positions at time t - delta t.

  	cout << "Classic Lennard-Jones fluid        " << endl;
  	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  	cout << "The program uses Lennard-Jones units " << endl;
  
  	ReadInput.open("MolecularDynamics_NVE/input/input.solid"); //Read input(in this particular case we read the file relative to the solid phase)
  	if (!ReadInput.is_open()) {
		cerr << "Error: unable to open input.solid" << endl;
		exit(1);
	}

 	ReadInput >> temp; //Read the value of temperature. All the quantities of the system are in LJ units, thus Temp = 1.1 means that k_b * T / \epsilon equals 1.1
 	cout << "Temperature = " << temp << endl;

  	ReadInput >> npart; //Read the number of particles
  	cout << "Number of particles = " << npart << endl;

  	ReadInput >> rho; //Read the numerical density: number of particles on volume
  	cout << "Density of particles = " << rho << endl;
  	vol = (double)npart/rho;
  	cout << "Volume of the simulation box = " << vol << endl;
  	box = pow(vol,1.0/3.0); //the side value is volume with power 1/3. The cell is cubic.
  	cout << "Edge of the simulation box = " << box << endl;

  	ReadInput >> rcut; // cut radius of the potential
  	ReadInput >> delta; // time integration step 
  	ReadInput >> nstep; //number of step that the algorithm has to do 
  	ReadInput >> iprint; // time in which I write on the screen
  	ReadInput >> restart; // restart from old configuration 
  	ReadInput >> print; // time in which I write to file

  	cout << "The program integrates Newton equations with the Verlet method " << endl;
  	cout << "Time step = " << delta << endl;
  	cout << "Number of steps = " << nstep << endl;
  	cout << "Number of bins = " << nbins << endl;
  	cout << "Restart from old configuration = " << restart << endl;
  	ReadInput.close();
  	
  	//vector inizializations:
  	//only useful in the ex7.3: computation of the radial distribution function
  	histogram.reserve(nbins);
  	g_block.reserve(nbins);
	ave_g.reserve(nbins);
	ave2_g.reserve(nbins);
	g.reserve(nbins);
	g_err.reserve(nbins);
	
	
	bin_size = (box/2.0)/(double)nbins;

  	
  	//at this point, the program, following the instructions in the input file decides if read an old configuration and restart the simulation or read config.0 and initialize the system.
  	if(restart){
  		Restart();
  	} else { Inizialization();}
  	
  	//Tail corrections for potential energy and pressure
  	vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  	ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  	cout << "Tail correction for the potential energy = " << vtail << endl;
  	cout << "Tail correction for the virial           = " << ptail << endl; //only useful in the ex7.3: computation of the radial distribution function
  	
}

void MolDyn_NVE :: Restart(void){

	//restart the simulation reading the old configuration x(t)
	ifstream ReadConf;
	
	//we follow the suggestions for the rescalation of the velocities given in the LSN_Exercises_04
	cout << "Reading initial configuraton from old.0" << endl;
	ReadConf.open("MolecularDynamics_NVE/config/old.0");
  	if (!ReadConf.is_open()) {
		cerr << "Error: unable to open old.0" << endl;
		exit(1);
	}
	
	for (int i=0; i<npart; ++i){
    		ReadConf >> x[i] >> y[i] >> z[i]; //reads the coordinates for each line
    		x[i] = x[i] * box;
    		y[i] = y[i] * box;
    		z[i] = z[i] * box; // multiply the coordinates by the length of the edge in LJ units. Since the coordinates in the file are stored with units of the lenght of the box edge, I obtain the coordinates in LJ units.
  	}
  	ReadConf.close();
  	
  	//read the configuration x(t-delta*t)
  	cout << "Reading initial configuration old.final" << endl;
  	ReadConf.open("MolecularDynamics_NVE/config/old.final");
  	if (!ReadConf.is_open()) {
		cerr << "Error: unable to open old.final" << endl;
		exit(1);
	}
	
  	for (int i=0; i<npart; ++i){
    		ReadConf >> x_old[i] >> y_old[i] >> z_old[i]; //reads the coordinates for each line
    		x_old[i] = x_old[i] * box;
    		y_old[i] = y_old[i] * box;
    		z_old[i] = z_old[i] * box; // multiply the coordinates by the legnth of the edge in LJ units. Since the coordinates in the file are stored with units of the lenght of the box edge, I obtain the coordinates in LJ units.
  	}
  	ReadConf.close();
  	
  	Move(); //calculate the configuration in t+delta*t 
  	
  	//calculate the velocities at t+delta*t/2
  	for(int i = 0; i < npart; i++){
  		vx[i] = Pbc(x[i] - x_old[i])/(delta);
  		vy[i] = Pbc(y[i] - y_old[i])/(delta);
  		vz[i] = Pbc(z[i] - z_old[i])/(delta);
  	}
  	
	//calculate the mean of square velocities
	double avev2 = 0.0;
	for(int i = 0; i < npart; i++){
		//square velocity of single particles
		double sumv2 = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		//mean square velocity for each run of the cycle
		avev2 = i/double(i+1)*avev2 + 1./double(i+1)*sumv2;
	}
	
	//Rescalation of the velocities in order to obtain the desired temperature
	double fs = sqrt(3 * temp / avev2);// fs = velocity scale factor, which fixes the temperature of the system. Thus the kinetic energy of the system corresponds to the right temperature.
   		
   	for (int i=0; i<npart; ++i){
     		vx[i] *= fs;
     		vy[i] *= fs;
     		vz[i] *= fs;

     		//correct the configuration x(t)
     		x_old[i] = Pbc(x[i] - vx[i] * delta);
     		y_old[i] = Pbc(y[i] - vy[i] * delta);
     		z_old[i] = Pbc(z[i] - vz[i] * delta);
     			
     	}
}
		
void MolDyn_NVE :: Inizialization(void){
	//initialize the system by reading the initial configuration of the particles 
	ifstream ReadConf;
	cout << "Reading the intial configuration from config.0 " << endl << endl; 
	ReadConf.open("MolecularDynamics_NVE/config/config.0");	
  	if (!ReadConf.is_open()) {
		cerr << "Error: unable to open config.0" << endl;
		exit(1);
	}
  	
  	for (int i=0; i<npart; ++i){
    		ReadConf >> x[i] >> y[i] >> z[i]; //reads the coordinats for each line
    		x[i] = x[i] * box;
    		y[i] = y[i] * box;
    		z[i] = z[i] * box; // multiply the coordinates by the length of the edge in LJ units. Since the coordinates in the file are stored with units of the lenght of the box edge, I obtain the coordinates in LJ units.
  	}
  	ReadConf.close();
  	
  	//Create a random number generator 
	Random Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in");
	
  	//prepare initial velocities
   	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl; //Here we could sample a Maxwell Boltzmann distribution given a certain temperature.
   	double sumv[3] = {0.0, 0.0, 0.0};
   	for (int i=0; i<npart; ++i){
     		vx[i] = Rand.Rannyu() - 0.5; 
     		vy[i] = Rand.Rannyu() - 0.5;
     		vz[i] = Rand.Rannyu() - 0.5;//velocities along the axis are between -0.5 and 0.5. These velocities must be rescaled in order they are coherent with the given temperature. In order to do this i have to calculate the kinetic energy.

     		sumv[0] += vx[i];
     		sumv[1] += vy[i];
     		sumv[2] += vz[i];
   		}
   		
   	for (int idim=0; idim<3; ++idim) 
   		sumv[idim] /= (double)npart; //divides the sum of the velocities along the three directions by the number of particles.
   		
   	double sumv2 = 0.0, fs;
  	for (int i=0; i<npart; ++i){
     		vx[i] = vx[i] - sumv[0];
     		vy[i] = vy[i] - sumv[1];
     		vz[i] = vz[i] - sumv[2]; // to each component of the velocity of each particle, I subtract the value of sumv for each particle. In this way the centre of mass is going to have velocity zero.

     		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   	}
   	
   	sumv2 /= (double)npart;

   	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor, which fixes the temperature of the system. Thus the kinetic energy of the system corresponds to the right temperature.
   		
   	for (int i=0; i<npart; ++i){
     		vx[i] *= fs;
     		vy[i] *= fs;
     		vz[i] *= fs;

     
     		x_old[i] = Pbc(x[i] - vx[i] * delta);
     		y_old[i] = Pbc(y[i] - vy[i] * delta);
     		z_old[i] = Pbc(z[i] - vz[i] * delta);
     			
     			
     	}
     	
  }
  
void MolDyn_NVE :: Move(void){ //Move particles with Verlet algorithm
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
	
	
	for(int i=0; i<npart; ++i){ //Force acting on particle i
    		fx[i] = Force(i,0);
    		fy[i] = Force(i,1);
    		fz[i] = Force(i,2);
  	}

  	for(int i=0; i<npart; ++i){ //Verlet integration scheme

    		xnew = Pbc( 2.0 * x[i] - x_old[i] + fx[i] * pow(delta,2) );
    		ynew = Pbc( 2.0 * y[i] - y_old[i] + fy[i] * pow(delta,2) );
    		znew = Pbc( 2.0 * z[i] - z_old[i] + fz[i] * pow(delta,2) );

    		vx[i] = Pbc(xnew - x_old[i])/(2.0 * delta);
    		vy[i] = Pbc(ynew - y_old[i])/(2.0 * delta);
		vz[i] = Pbc(znew - z_old[i])/(2.0 * delta);

    		x_old[i] = x[i];
    		y_old[i] = y[i];
    		z_old[i] = z[i];

    		x[i] = xnew;
    		y[i] = ynew;
    		z[i] = znew;
  	}

}


double MolDyn_NVE :: Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
	double f=0.0; //force in the direction idir
  	double dvec[3], dr; // dvec vector joining two particles, dr distance between two particles

  	for (int i=0; i<npart; ++i){
    		if(i != ip){ //ip is the index of the particle under analisys. In this way we perform a cycle on all the particles that are not the particle ip, calculating the force which acts on the ip particle.
      			dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc, with the minimum image convention. It means that we are calculating the distance between the ip particle and another particle, which is the closest one to ip particle.
      			dvec[1] = Pbc( y[ip] - y[i] );
      			dvec[2] = Pbc( z[ip] - z[i] );

      			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      			dr = sqrt(dr);

      			if(dr < rcut){ // condition in which the couple of particles is interacting
        			f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r). Calculating the contribution of the particle i to the force acting on the particle ip
      			}
    		}
  	}
  
  return f;
}

void MolDyn_NVE :: ConfFinal(void){ //Write final configuration
  	ofstream WriteConf;

  	cout << "Print final configuration to file config.final " << endl << endl;
  	WriteConf.open("MolecularDynamics_NVE/config/config.final");
  	if (!WriteConf.is_open()) {
		cerr << "Error: unable to open config.final" << endl;
		exit(1);
	}
	
  	for (int i=0; i<npart; ++i){
    		WriteConf << x[i] / box << "\t" <<  y[i] / box << "\t" << z[i] / box << endl;
  }
  	WriteConf.close();
}

void MolDyn_NVE :: SaveOldConfig(void){
	//Save old configuration obtained during the initialization 
	ofstream WriteConf_old;
	cout << "Print final configuration r(t-dt) to file old.final " << endl << endl;
  	WriteConf_old.open("MolecularDynamics_NVE/config/old.final");
  	if (!WriteConf_old.is_open()) {
		cerr << "Error: unable to open WriteConf_old.final" << endl;
		exit(1);
	}
	
  	for (int i=0; i < npart; i++) {
  		WriteConf_old << x_old[i] / box << "\t" << y_old[i] / box << "\t" << z_old[i] / box << endl;
  		}
  		
  	WriteConf_old.close();
  	
  	cout << "Print final configuration r(t) to file old.0 " << endl << endl;
  	WriteConf_old.open("MolecularDynamics_NVE/config/old.0");
  	if (!WriteConf_old.is_open()) {
		cerr << "Error: unable to open WriteConf_old.0" << endl;
		exit(1);
	}
	for (int i=0; i < npart; i++) {
  		WriteConf_old << x[i] / box << "\t" << y[i] / box << "\t" << z[i] / box << endl;
  		}
  		
  	WriteConf_old.close();
}


void MolDyn_NVE :: ConfXYZ(int nconf){ //Write configuration in .xyz format
  	ofstream WriteXYZ;

  	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  	if (!WriteXYZ.is_open()) {
		cerr << "Error: unable to open frames/config_" + to_string(nconf) + ".xyz" << endl;
		exit(1);
	}
  	WriteXYZ << npart << endl;
  	WriteXYZ << "This is only a comment!" << endl;
  	for (int i=0; i<npart; ++i){
    		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  	}
  	WriteXYZ.close();
}
  

void MolDyn_NVE :: Measure(void){ //Properties measurement
  	//int bin;
  	double v = 0., t = 0., vij = 0., vir = 0., virij; // define total potential energy, total kinetic energy and the interaction potential between particles i and j. The virial vir is useful only for the last point
  	double dx, dy, dz, dr;
 
 	//reset the histrogram to zero.
 	for (int i = 0; i<nbins; i++){
 		histogram[i] = 0;
 	}
 	
  	
  	//cycle over pairs of particles
  	for (int i=0; i<npart-1; ++i){
    		for (int j=i+1; j<npart; ++j){

     			dx = Pbc( x_old[i] - x_old[j] ); 
     			dy = Pbc( y_old[i] - y_old[j] ); 
     			dz = Pbc( z_old[i] - z_old[j] ); 

     			dr = dx*dx + dy*dy + dz*dz;
    			dr = sqrt(dr);
    			
    			//counting the particles and updating the histogram
    			int index = int(dr/box*2*nbins);
						
			if (index < nbins and index >= 0){
				histogram[index] ++;
				histogram[index] ++;
			}

     			if(dr < rcut){
       				vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       				virij = 48./3. * (1./pow(dr,12) - 0.5/pow(dr,6));
       				
				//Potential energy
       				v += vij;
       				vir += virij;
     			}
    		}          
  	}
  	
  	//Kinetic energy
  	for (int i=0; i<npart; ++i) 
  		t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  		
  	double est_pot = v/(double)npart; //Potential energy per particle
    	double est_kin = t/(double)npart; //Kinetic energy per particle
    	double est_temp = (2.0 / 3.0) * est_kin; //Temperature
    	double est_etot = est_pot + est_kin; //Total energy per particle
    	
    	
    	blockEpot = iblock/double(iblock + 1) * blockEpot + 1./double(iblock + 1) * est_pot;
    	blockEkin = iblock/double(iblock + 1) * blockEkin + 1./double(iblock + 1) * est_kin;
    	blockEtot = iblock/double(iblock + 1) * blockEtot + 1./double(iblock + 1) * est_etot;
    	blockTemp = iblock/double(iblock + 1) * blockTemp + 1./double(iblock + 1) * est_temp;
    	blockVir = iblock/double(iblock + 1) * blockVir + 1./double(iblock + 1) * vir;
    	
    	//computation of g(r)
    	for(int i = 0; i < nbins; i++){
  		double normalization = 4./3. * pi * rho * npart *  (pow((i+1)*bin_size, 3) - pow(i*bin_size, 3));
  	 	
  		g[i] = histogram[i] / normalization;
  		//cout << " g : " << g[i] << endl;
  		
  	}
  	
  	for(int i = 0; i < nbins; i++){
  		g_block[i]= iblock/double(iblock+1)*g_block[i]+1./double(iblock+1)*g[i];
  		//cout << " g_block : " << g_block[i] << endl;
  		 
  	}
    	
    	iblock++; //definitions useful in the second and third exercises
    	
    	
    	
    	if(print){
    		ofstream Epot, Ekin, Etot, Temp, Pres;

  		Epot.open("data/output_epot.dat", ios::app);
  		if (!Epot.is_open()) {
			cerr << "Error: unable to open output_epot.dat" << endl;
			exit(1);
		}
  		Ekin.open("data/output_ekin.dat",ios::app);
  		if (!Ekin.is_open()) {
			cerr << "Error: unable to open output_ekin.dat" << endl;
			exit(1);
		}
  		Temp.open("data/output_temp.dat",ios::app);
  		if (!Temp.is_open()) {
			cerr << "Error: unable to open output_temp.dat" << endl;
			exit(1);
		}
  		Etot.open("data/output_etot.dat",ios::app);
  		if (!Etot.is_open()) {
			cerr << "Error: unable to open output_etot.dat" << endl;
			exit(1);
		}
		
		Pres.open("data/output_pres.dat",ios::app);
  		if (!Etot.is_open()) {
			cerr << "Error: unable to open output_pres.dat" << endl;
			exit(1);
		}
  	
  		Epot << est_pot  << endl;
    		Ekin << est_kin  << endl;
    		Temp << est_temp << endl;
    		Etot << est_etot << endl;
    		Pres << rho*temp+vir/vol << endl;
    		
    		Epot.close();
    		Ekin.close();
   		Temp.close();
    		Etot.close();
    		Pres.close();
	}
	
}

    		
void MolDyn_NVE :: Blocking(void){   //useful in the last two exercises 	
	double errorEpot, errorEkin, errorEtot, errorTemp, errorVir;
	
	aveEpot = n_blocks/double(n_blocks + 1) * aveEpot + 1./double(n_blocks + 1) * blockEpot;
	aveEkin = n_blocks/double(n_blocks + 1) * aveEkin + 1./double(n_blocks + 1) * blockEkin;
	aveEtot = n_blocks/double(n_blocks + 1) * aveEtot + 1./double(n_blocks + 1) * blockEtot;
	aveTemp = n_blocks/double(n_blocks + 1) * aveTemp + 1./double(n_blocks + 1) * blockTemp;
	aveVir = n_blocks/double(n_blocks + 1) * aveVir + 1./double(n_blocks + 1) * blockVir;
	
	
	ave2Epot = n_blocks/double(n_blocks + 1) * ave2Epot + 1./double(n_blocks + 1) * blockEpot * blockEpot;
	ave2Ekin = n_blocks/double(n_blocks + 1) * ave2Ekin + 1./double(n_blocks + 1) * blockEkin * blockEkin;
	ave2Etot = n_blocks/double(n_blocks + 1) * ave2Etot + 1./double(n_blocks + 1) * blockEtot * blockEtot;
	ave2Temp = n_blocks/double(n_blocks + 1) * ave2Temp + 1./double(n_blocks + 1) * blockTemp * blockTemp;
	ave2Vir = n_blocks/double(n_blocks + 1) * ave2Vir + 1./double(n_blocks + 1) * blockVir * blockVir;
	
	for(int i = 0; i < nbins; i++){
   		ave_g[i] = n_blocks/(double)(n_blocks+1)*ave_g[i] + 1./(double)(n_blocks+1)*g_block[i];
   		//cout << " media g : " << ave_g[i] << endl;
   		ave2_g[i] = n_blocks/(double)(n_blocks+1)*ave_g[i] + 1./(double)(n_blocks+1)*g_block[i]*g_block[i];
   		//cout << " media2 g : " << ave2_g[i] << endl;
   	}
	
	if(n_blocks == 0) {
		errorEpot = 0;
		errorEkin = 0;
		errorEtot = 0;
		errorTemp = 0;
		errorVir = 0;
		for(int i = 0; i < nbins; i++){
   			g_err[i] = 0;
   		}
	} else { 
		errorEpot = sqrt((ave2Epot - aveEpot * aveEpot)/n_blocks);
		errorEkin = sqrt((ave2Ekin - aveEkin * aveEkin)/n_blocks);
		errorEtot = sqrt((ave2Etot - aveEtot * aveEtot)/n_blocks);
		errorTemp = sqrt((ave2Temp - aveTemp * aveTemp)/n_blocks);
		errorVir = sqrt((ave2Vir - aveVir * aveVir)/n_blocks);
		for(int i = 0; i < nbins; i++){
   			g_err[i] = sqrt(abs(ave2_g[i] - ave_g[i] * ave_g[i])/n_blocks);
	}
	
   		
   	}
	
	n_blocks++;
	
	double P = rho*temp+aveVir/vol;
	double errorP = errorVir/vol;
	
	ofstream ave_Epot, ave_Ekin, ave_Temp, ave_Etot, ave_Pres, Gofr, Gave;
	
  	ave_Epot.open("data/output_ave_epot.dat",ios::app);
  	if (!ave_Epot.is_open()) {
			cerr << "Error: unable to open output_ave_epot.dat" << endl;
			exit(1);
		}
  	ave_Ekin.open("data/output_ave_ekin.dat",ios::app);
  	if (!ave_Ekin.is_open()) {
			cerr << "Error: unable to open output_ave_ekin.dat" << endl;
			exit(1);
		}
  	ave_Temp.open("data/output_ave_temp.dat",ios::app);
  	if (!ave_Temp.is_open()) {
			cerr << "Error: unable to open output_ave_temp.dat" << endl;
			exit(1);
		}
  	ave_Etot.open("data/output_ave_etot.dat",ios::app);
  	if (!ave_Etot.is_open()) {
			cerr << "Error: unable to open output_ave_etot.dat" << endl;
			exit(1);
		}
	
	ave_Pres.open("data/output_ave_pres.dat",ios::app);
  	if (!ave_Etot.is_open()) {
			cerr << "Error: unable to open output_ave_pres.dat" << endl;
			exit(1);
		}
		
	
	ave_Epot << n_blocks << "\t" << aveEpot << "\t" << errorEpot  << endl;
    	ave_Ekin << n_blocks << "\t" << aveEkin << "\t" << errorEkin  << endl;
    	ave_Temp << n_blocks << "\t" << aveTemp << "\t" << errorTemp  << endl;
   	ave_Etot << n_blocks << "\t" << aveEtot << "\t" << errorEtot  << endl;
   	ave_Pres << n_blocks << "\t" << P << "\t" << errorP  << endl;
   	
   	
   	Gofr.open("data/gofr_0.dat",ios::app);
    	if (!Gofr.is_open()) {
		cerr << "Error: unable to open gofr_0.dat" << endl;
		exit(1);
	}
	Gofr << "___NEW___BLOCK___NUMBER___" << iblock << endl;
	for(int i = 0; i < nbins; i++){
		Gofr << iblock << "\t" << i*bin_size << "\t" << (i+1)*bin_size << "\t" << g_block[i] << endl;
		//cout << " GBLOCK : " << g_block[i] << endl;
	}
   	Gave.open("data/gave_0.dat",ios::app);
   	if (!Gave.is_open()) {
		cerr << "Error: unable to open gave.dat" << endl;
		exit(1);
	}
	
	if(n_blocks == nstep/1000){//into gave I put g(r) g_err, rmin and rmax
		for (int i = 0; i < nbins; i++){
			
			Gave << ave_g[i] << "\t" << g_err[i] << "\t" << i*bin_size << "\t" << (i+1)*bin_size << endl;
		}
	}
		
   	
    		
   	ave_Epot.close();
   	ave_Ekin.close();
   	ave_Temp.close();
   	ave_Etot.close();
   	ave_Pres.close();
   	Gofr.close();
   	Gave.close();
   	
   	blockEpot = 0.;
   	blockEkin = 0.;
   	blockEtot = 0.;
   	blockTemp = 0.;
   	blockVir = 0.;
   	iblock = 0;
   	for (int i = 0; i < nbins; i++){
   		g_block[i] = 0;
   	}
   		
}


double MolDyn_NVE :: Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    	return r - box * rint(r/box);
}

int MolDyn_NVE :: Steps(void){ // take the steps of the simulation
	return nstep;
}

int MolDyn_NVE :: Prints(void){ // take the steps of the printing
	return iprint;
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
