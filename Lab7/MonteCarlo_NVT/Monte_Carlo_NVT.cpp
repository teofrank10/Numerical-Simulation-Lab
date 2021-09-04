/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "Monte_Carlo_NVT.h"

using namespace std;

Monte_Carlo_NVT :: Monte_Carlo_NVT(): Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in"){ }

Monte_Carlo_NVT :: ~Monte_Carlo_NVT(){ }

void Monte_Carlo_NVT :: Input(void){

	ifstream ReadInput;

  	cout << "Classic Lennard-Jones fluid        " << endl;
  	cout << "Monte Carlo simulation             " << endl << endl;
  	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  	cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  	cout << "The program uses Lennard-Jones units " << endl;
  
  	string dir = "MonteCarlo_NVT/input/";
  	
  	string file = "input.solid";
  	
  	
  
  	cout << "The program simulates " << file.substr(6,file.length()) << " phase" << endl;
  
	//Read input informations
  	ReadInput.open(dir+file);
  	if (!ReadInput.is_open()) {
		cerr << "Error: unable to open " << file << endl;
		exit(1);
	}

  	ReadInput >> temp;
  	beta = 1.0/temp;
  	cout << "Temperature = " << temp << endl;

  	ReadInput >> npart;
  	cout << "Number of particles = " << npart << endl;

  	ReadInput >> rho;
  	cout << "Density of particles = " << rho << endl;
  	vol = (double)npart/rho;
  	box = pow(vol,1.0/3.0);
  	cout << "Volume of the simulation box = " << vol << endl;
  	cout << "Edge of the simulation box = " << box << endl;

  	ReadInput >> rcut;
  	cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  	//Tail corrections for potential energy and pressure
  	vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  	ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  	cout << "Tail correction for the potential energy = " << vtail << endl;
  	cout << "Tail correction for the virial           = " << ptail << endl; 

  	ReadInput >> delta;

  	ReadInput >> nblk;

  	ReadInput >> nstep;
  
  	ReadInput >> neqstep;
  
  	ReadInput >> restart;
  	
  	ReadInput >> print;//variable useful to decide whether to print instantaneous values of energy and pressure (print = 1)  or not (print=0)

	ReadInput >> nbins; //variable useful for the computation of te histogram of g(r)
  	cout << "The program perform Metropolis moves with uniform translations" << endl;
  	cout << "Moves parameter = " << delta << endl;
  	cout << "Number of blocks = " << nblk << endl;
  	cout << "Number of steps in one block = " << nstep << endl << endl;
  	cout << "Number of bins = " << nbins << endl;
  	cout << "The program restarts from an old configuration: " << restart << endl;
  	ReadInput.close();
  	
  	iblock = 0;
  	istep = 0;
  	
  	//vector inizializations:
  	histogram.reserve(nbins);
  	g_block.reserve(nbins);
	ave_g.reserve(nbins);
	ave2_g.reserve(nbins);
	g.reserve(nbins);
	g_err.reserve(nbins);
  	
  	
  	//measurement of g(r)  
  	//we see the histrogram as a vector: every element of the vector is a bin 
  	bin_size = (box/2.0)/(double)nbins;
  	
  	
  	if(restart){
  		Restart();
  	}else{ Initialization();};
 
}

void Monte_Carlo_NVT :: Initialization(void){

	//Read initial configuration
	
	ifstream ReadConf;
	
  	cout << "Read initial configuration from file config.0 " << endl << endl;
  	cout << "Equilibrate the system " << endl << endl;
  	ReadConf.open("MonteCarlo_NVT/config/config.0");
  	if (!ReadConf.is_open()) {
		cerr << "Error: unable to open config.0" << endl;
		exit(1);
	}
  	for (int i=0; i<npart; ++i){
    		ReadConf >> x[i] >> y[i] >> z[i];
    		x[i] = Pbc( x[i] * box );//the dimension in the files are in LJ units: we have to multilpy by box
    		y[i] = Pbc( y[i] * box );
    		z[i] = Pbc( z[i] * box );
  	}
  	ReadConf.close();
  	
  	//Equilibration with neqstep
  	for (int i = 0; i < neqstep; i++)
		Move();
  	
  	//after the equilibration we have to reset the counters that have been used in the method: 
  	accepted = 0;
  	attempted = 0;
  	
	//Evaluate potential energy and virial of the initial configuration after the equilibration
  	Measure();

	
//Print initial values for the potential energy and virial
  	cout << "Initial potential energy (with tail corrections) = " << v_blocking/(double)npart + vtail << endl;
  	cout << "Virial                   (with tail corrections) = " << w_blocking/(double)npart + ptail << endl;
  	cout << "Pressure                 (with tail corrections) = " << rho * temp + (w_blocking + (double)npart * ptail) / vol << endl << endl;

}


void Monte_Carlo_NVT :: Restart(void){

	//instead of reading from config.0 read fro an old config and do not execute the equilibration
	ifstream ReadConf;
	
	cout << "Read initial configuration from file config.final " << endl << endl;
	cout << "The system is already equilibrated " << endl << endl;
	ReadConf.open("MonteCarlo_NVT/config/config.final");
	if (!ReadConf.is_open()) {
		cerr << "Error: unable to open config.final" << endl;
		exit(1);
	}
	
	for (int i=0; i<npart; ++i){
    		ReadConf >> x[i] >> y[i] >> z[i];
    		x[i] = Pbc( x[i] * box );//the dimension in the files are in LJ units: we have to multilpy by box
    		y[i] = Pbc( y[i] * box );
    		z[i] = Pbc( z[i] * box );
  	}
	ReadConf.close();
	
	//Evaluate potential energy and virial of the initial configuration without equilibration
  	Measure();

	
//Print initial values for the potential energy and virial
  	cout << "Initial potential energy (with tail corrections) = " << v_blocking/double(npart) + vtail << endl;
  	cout << "Virial                   (with tail corrections) = " << w_blocking/double(npart)  + ptail << endl;
  	cout << "Pressure                 (with tail corrections) = " << rho * temp + (w_blocking + (double)npart * ptail) / vol << endl << endl;

}

void Monte_Carlo_NVT :: Move(void){

	int o;
  	double p, energy_old, energy_new;
  	double xold, yold, zold, xnew, ynew, znew;


  	for(int i=0; i<npart; ++i){
  		//Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
    		o = (int)(Rand.Rannyu()*npart);

  		//Old
    		xold = x[o];
    		yold = y[o];
    		zold = z[o];

    		energy_old = Boltzmann(xold,yold,zold,o);

  		//New
    		xnew = Pbc( x[o] + delta*(Rand.Rannyu() - 0.5) );
    		ynew = Pbc( y[o] + delta*(Rand.Rannyu() - 0.5) );
    		znew = Pbc( z[o] + delta*(Rand.Rannyu() - 0.5) );

    		energy_new = Boltzmann(xnew,ynew,znew,o);

  		//Metropolis test
    		p = exp(beta*(energy_old-energy_new));
    		if(p >= Rand.Rannyu()){
    		
    			//Update
       			x[o] = xnew;
       			y[o] = ynew;
       			z[o] = znew;
    
       			accepted += 1.0;
    		}
    		attempted += 1.0;
  	}
  	
}


double Monte_Carlo_NVT :: Boltzmann(double xx, double yy, double zz, int ip){
  	double ene=0.0;
  	double dx, dy, dz, dr;

  	for (int i=0; i<npart; ++i){
  	
    		if(i != ip){
			// distance ip-i in pbc
      			dx = Pbc(xx - x[i]);
      			dy = Pbc(yy - y[i]);
      			dz = Pbc(zz - z[i]);

      			dr = dx*dx + dy*dy + dz*dz;
      			dr = sqrt(dr);

      			if(dr < rcut){
       				ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      			}
    		}
  	}

  return 4.0*ene;
}		
	
	
void Monte_Carlo_NVT :: Measure(void){
	
 	//reset the histrogram to zero.
 	for (int i = 0; i<nbins; i++){
 		histogram[i] = 0;
 	}
 	//fill(histogram.begin(), histogram.end(), 0);
 	
  	double v = 0.0, w = 0.0;
  	double vij, wij;
  	double dx, dy, dz, dr;

	//instead of using a vector "walker", we use instantaneous variables and execute 
	//reset the hystogram of g(r)
  	/*for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;*/

	//cycle over pairs of particles
  	for (int i=0; i<npart-1; ++i){
    		for (int j=i+1; j<npart; ++j){

			// distance i-j in pbc
     			dx = Pbc(x[i] - x[j]);
     			dy = Pbc(y[i] - y[j]);
     			dz = Pbc(z[i] - z[j]);

     			dr = dx*dx + dy*dy + dz*dz;
     			dr = sqrt(dr);

			//update of the histogram of g(r)
			//we increase the number of elements in the bin of two
			int index = int(dr/box*2*nbins);
						
			if (index < nbins and index >= 0){
				histogram[index] ++;
				histogram[index] ++;
			}	

     			if(dr < rcut){
       				vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       				wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

				// contribution to energy and virial
       				v += vij;
       				w += wij;
     			}
    		}          
 	}

  	//walker[iv] = 4.0 * v;
  	v *= 4.0;
  	w *= 16.0;
  	//walker[iw] = 48.0 * w / 3.0;
  	  	 
  	for(int i = 0; i < nbins; i++){
  		double normalization = 4./3. * pi * rho * npart *  (pow((i+1)*bin_size, 3) - pow(i*bin_size, 3));
  	 	
  		g[i] = histogram[i] / normalization;
  		//cout << " g : " << g[i] << endl;
  		
  	}
  	
  	
  	//compute the block average on fly
  	
  	v_blocking = istep/(double)(istep + 1)*v_blocking + 1./(istep + 1)*v;
  	w_blocking = istep/(double)(istep + 1)*w_blocking + 1./(istep + 1)*w;
  	
  	  	
  	for(int i = 0; i < nbins; i++){
  		g_block[i]= istep/double(istep+1)*g_block[i]+1./double(istep+1)*g[i];
  		//cout << " g_block : " << g_block[i] << endl;
  		 
  	}
  	
  	
  	istep++;
  	
  	if(print){
  	
  		//print the instantaneous values of energy and pressure
  	
  		ofstream outfile_en, outfile_pres;
  	
  		outfile_en.open("data/instant_en.dat", ios::app);
  		if (!outfile_en.is_open()) {
			cerr << "Error: unable to open instant_en.dat" << endl;
			exit(1);
		}
	
		outfile_pres.open("data/instant_pres.dat", ios::app);
  		if (!outfile_pres.is_open()) {
			cerr << "Error: unable to open instant_pres.dat" << endl;
			exit(1);
		}
	
		
		outfile_en << v/double(npart)  + vtail << endl;
		outfile_pres << (w + ptail*npart)/vol + temp*rho << endl;
	
		outfile_en.close();
		outfile_pres.close();
		
	}
}

void Monte_Carlo_NVT :: Accumulate(void){ //Update block averages calculating them on fly

   	ave_v = iblock/(double)(iblock+1)*ave_v + 1./(double)(iblock+1)*v_blocking;
   	ave2_v = iblock/(double)(iblock+1)*ave2_v + 1./(double)(iblock+1)*v_blocking*v_blocking;
    	ave_w = iblock/(double)(iblock+1)*ave_w + 1./(double)(iblock+1)*w_blocking;
   	ave2_w = iblock/(double)(iblock+1)*ave2_w + 1./(double)(iblock+1)*w_blocking*w_blocking;
   	
   	for(int i = 0; i < nbins; i++){
   		ave_g[i] = iblock/(double)(iblock+1)*ave_g[i] + 1./(double)(iblock+1)*g_block[i];
   		//cout << " media g : " << ave_g[i] << endl;
   		ave2_g[i] = iblock/(double)(iblock+1)*ave_g[i] + 1./(double)(iblock+1)*g_block[i]*g_block[i];
   		//cout << " media2 g : " << ave2_g[i] << endl;
   	}
   	
   	iblock += 1;

}

void Monte_Carlo_NVT :: Averages(void){ //Print results for current block

    	
   	double v_err = Error(ave_v, ave2_v, iblock);
   	double w_err = Error(ave_w, ave2_w, iblock);
   	for(int i = 0; i < nbins; i++){
   		g_err[i] = Error(ave_g[i], ave2_g[i], iblock);
   		//cout << " errore g : " << g_err[i] << endl;
   	}
   	
   
   	double pres = (ave_w + ptail*npart)/vol + temp*rho;
   	double pres_err = w_err/vol;
   	
   //double r, gdir;
   	ofstream Gofr, Gave, Epot, Pres;
   	const int wd=12;
    
   	cout << "Block number " << iblock << endl;
    	cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    	Epot.open("data/epot_0.dat",ios::app);
    	if (!Epot.is_open()) {
		cerr << "Error: unable to open epot_0.dat" << endl;
		exit(1);
	}
    	Pres.open("data/pres_0.dat",ios::app);
    	if (!Pres.is_open()) {
		cerr << "Error: unable to open pres_0.dat" << endl;
		exit(1);
	}
	
	Gofr.open("data/gofr_0.dat",ios::app);
    	if (!Gofr.is_open()) {
		cerr << "Error: unable to open gofr_0.dat" << endl;
		exit(1);
	}
	
	Epot << iblock << "\t" << ave_v/double(npart)  + vtail << "\t" << v_err/(double)npart << endl;
	
	Pres << iblock << "\t" << pres << "\t" << pres_err << endl;
	
	Gofr << "___NEW___BLOCK___NUMBER___" << iblock << endl;
	for(int i = 0; i < nbins; i++){
		Gofr << iblock << "\t" << i*bin_size << "\t" << (i+1)*bin_size << "\t" << g_block[i] << endl;
		//cout << " GBLOCK : " << g_block[i] << endl;
	}
	Epot.close();
	Pres.close();
	Gofr.close();
	
	cout << "============================" << endl << endl;
	
	Reset();
	
	//print the results using only the value of the last block
	
	if( iblock == nblk){
		ofstream results;
		results.open("data/results.dat", ios::app);
		if (!results.is_open()) {
			cerr << "Error: unable to open results.dat" << endl;
			exit(1);
		}
		Gave.open("data/gave_0.dat",ios::app);
    		if (!Gave.is_open()) {
			cerr << "Error: unable to open gave.dat" << endl;
			exit(1);
		}
		
		results << ave_v/double(npart)  + vtail << "\t" << setw(wd) << v_err/(double)npart << "\t" << setw(wd) << pres << "\t" << setw(wd) << pres_err << endl;
		//into gave I put g(r) g_err, rmin and rmax
		for (int i = 0; i < nbins; i++){
			
			Gave << ave_g[i] << "\t" << g_err[i] << "\t" << i*bin_size << "\t" << (i+1)*bin_size << endl;
		}
		
		results.close();
		Gave.close();
	}
}
	
   
void Monte_Carlo_NVT :: Reset(void){ //Reset block averages

   	v_blocking = 0;
   	w_blocking = 0;
   	istep = 0;
   	attempted = 0;
   	accepted = 0;
   	//fill(g_block.begin(), g_block.end(), 0);
   	for (int i = 0; i < nbins; i++){
   		g_block[i] = 0;
   	}
}


void Monte_Carlo_NVT :: ConfFinal(void){
  	ofstream WriteConf;

  	cout << "Print final configuration to file config.final " << endl << endl;
  	WriteConf.open("MonteCarlo_NVT/config/config.final");
  	if (!WriteConf.is_open()) {
			cerr << "Error: unable to open config.final" << endl;
			exit(1);
	}
  	for (int i=0; i<npart; ++i){
  	
    		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  	}
  	WriteConf.close();

  	Rand.SaveSeed();
}

void Monte_Carlo_NVT :: ConfXYZ(int nconf){ //Write configuration in .xyz format
  	ofstream WriteXYZ;

  	WriteXYZ.open("MonteCarlo_NVT/frames/config_" + to_string(nconf) + ".xyz");
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

double Monte_Carlo_NVT :: Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Monte_Carlo_NVT :: Error(double ave, double ave2, int iblk)
{
    if( iblk == 1 or iblk == 0) return 0.0;
    else return sqrt(abs(ave2- pow(ave,2))/(double)(iblk-1));
}

int Monte_Carlo_NVT :: GetNumBlocks(){ return nblk; }

int Monte_Carlo_NVT :: GetNumStepsPerBlock(){ return nstep; }


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
