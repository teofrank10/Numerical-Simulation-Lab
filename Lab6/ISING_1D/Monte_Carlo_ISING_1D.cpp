#include "Monte_Carlo_ISING_1D.h"

using namespace std;

Monte_Carlo_ISING_1D :: Monte_Carlo_ISING_1D(): Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in"){ }

Monte_Carlo_ISING_1D :: ~Monte_Carlo_ISING_1D(){ }

void Monte_Carlo_ISING_1D :: Input(void){

	ifstream ReadInput;

	cout << "Classic 1D Ising model             " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
  	cout << "Nearest neighbour interaction      " << endl << endl;
  	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  	cout << "The program uses k_B=1 and mu_B=1 units " << endl;
  
//Read input informations
  	ReadInput.open("ISING_1D/input.dat");
  	if (!ReadInput.is_open()) {
		cerr << "Error: unable to open input.dat" << endl;
		exit(1);
	}

	ReadInput >> temp;
  	
  	beta = 1.0/temp;
  	cout << "Temperature = " << temp << endl;

  	ReadInput >> nspin;
  	cout << "Number of spins = " << nspin << endl;

  	ReadInput >> J;
  	cout << "Exchange interaction = " << J << endl;

  	ReadInput >> h;
  	cout << "External field = " << h << endl << endl;
    
  	ReadInput >> metro; // if=1 Metropolis, else Gibbs

  	ReadInput >> nblk;//num of blocks

  	ReadInput >> nstep;//num of steps for each block
  	
  	ReadInput >> n_eqstep;//num of equilibration steps 
  	
  	ReadInput >> restart;//restart from a previous spin configuration

  	if(metro==1) cout << "The program perform Metropolis moves" << endl;
  	else cout << "The program perform Gibbs moves" << endl;
  	cout << "Number of blocks = " << nblk << endl;
  	cout << "Number of steps in one block = " << nstep << endl << endl;
  	
  	beta = 1./temp;
  	
  	ReadInput.close();

	
  	if(restart){
  		Restart();
  	}else{Inizialization();}
}

void Monte_Carlo_ISING_1D :: Restart(void){
	//read old configuration

	ifstream ReadConf;
	cout << "Reading the initial configuration from config.final" << endl;
	ReadConf.open("config.final");
	if (!ReadConf.is_open()) {
		cerr << "Error: unable to open config.final" << endl;
		exit(1);
	}
  	for (int i=0; i<nspin; ++i){//read the previous configuration
  	
    		ReadConf >> s[i];
    	}
    	
    	ReadConf.close();
    
}

void Monte_Carlo_ISING_1D :: Inizialization(void){

	for (int i=0; i<nspin; ++i){//create the configuration: this corresponds to the case T=inf
    		
    		if(Rand.Rannyu() >= 0.5) s[i] = 1;
    		else s[i] = -1;
  	}
 
  	for (int i = 0; i < n_eqstep; i++){//equilibration
  	
  		Move();
  	}
  	
  	attempted = 0;
  	accepted = 0;
  	
}


void Monte_Carlo_ISING_1D :: Move(){


	if(metro==1){
  		for(int i=0; i<nspin; ++i){
  		//choose a particle between 0 and nspin-1
   			int o = (int)(Rand.Rannyu()*nspin);

    			Metropolis(o);
    		}

    	}else{//Gibbs sampling
    		for(int i=0; i<nspin; ++i){

			Gibbs(i);
    		}
  	}
  	
}

void Monte_Carlo_ISING_1D :: Metropolis(int o){

	attempted++;
	double alpha = fmin(1, exp(-beta * (2 * J * s[o] * (s[Pbc(o-1)] + s[Pbc(o+1)]) + 2 * h * s[o])));
	double r = Rand.Rannyu();
	if(r < alpha){
		s[o] = -s[o];
		accepted++;
	}//otherwise we keep the same configuration
	
}

void Monte_Carlo_ISING_1D :: Gibbs(int i){
	attempted++;
	accepted++;
	double prob = 1./(1.+ exp(- beta * ( 2 * J * (s[Pbc(i-1)] + s[Pbc(i+1)]) + 2 * h)));
	double r = Rand.Rannyu();
	if(r < prob){
		s[i] = 1;
	}else{
		s[i] = -1;
	}
	
}

void Monte_Carlo_ISING_1D :: Measure(){


	double u = 0.0, m = 0.0;//initialise the single values of hamiltonian and magetization

	//cycle over spins
  	for (int i=0; i<nspin; ++i){
     		u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		m += s[i];
  	}
  	
  	//calculate the single values in a single block
  	block_U = step_run/double(step_run+1)*block_U + 1./double(step_run+1)*u;
  	block_U2 = step_run/double(step_run+1)*block_U2 + 1./double(step_run+1)*u*u;
  	block_M = step_run/double(step_run+1)*block_M + 1./double(step_run+1)*m;
  	block_M2 = step_run/double(step_run+1)*block_M2 + 1./double(step_run+1)*m*m;
  	
  	step_run++;	
}
  		

void Monte_Carlo_ISING_1D :: MeasureBlock(void){
	double error_U, error_M, error_C, error_Chi;
	
	//calculate the averages on single blocks 
	ave_U = block_run/double(block_run+1)*ave_U + 1./double(block_run+1)*block_U;
	ave2_U= block_run/double(block_run+1)*ave2_U + 1./double(block_run+1)*block_U*block_U;
	ave_M = block_run/double(block_run+1)*ave_M + 1./double(block_run+1)*block_M;
	ave2_U = block_run/double(block_run+1)*ave2_M + 1./double(block_run+1)*block_M*block_M;	
	
	double C = beta*beta*(block_U2-block_U*block_U);
	double Chi = 0;
	
	if(h == 0){
	
		Chi = beta * block_M2;
	}else{
		Chi = beta * (block_M2-block_M*block_M);
	}
	
	ave_C = block_run/double(block_run+1)*ave_C + 1./double(block_run+1)*C;
	ave2_C = block_run/double(block_run+1)*ave2_C + 1./double(block_run+1)*C*C;
	ave_Chi = block_run/double(block_run+1)*ave_Chi + 1./double(block_run+1)*Chi;
	ave2_Chi = block_run/double(block_run+1)*ave2_Chi + 1./double(block_run+1)*Chi*Chi;
  	
  	error_U = Error(ave_U, ave2_U, block_run);
  	error_M = Error(ave_M, ave2_M, block_run);
  	error_C = Error(ave_C, ave2_C, block_run);
  	error_Chi = Error(ave_Chi, ave2_Chi, block_run);
  	
  	block_run++;
  	
  	ofstream Ustream, Mstream, Cstream, Chistream;
  	
  	Ustream.open("data/ave_U.dat", ios:: app);
  	if (!Ustream.is_open()) {
		cerr << "Error: unable to open ave_U.dat" << endl;
		exit(1);
		}
  	Mstream.open("data/ave_M.dat", ios:: app);
  	if (!Mstream.is_open()) {
		cerr << "Error: unable to open ave_M.dat" << endl;
		exit(1);
		}
  	Cstream.open("data/ave_C.dat", ios:: app);
  	if (!Cstream.is_open()) {
		cerr << "Error: unable to open ave_C.dat" << endl;
		exit(1);
		}
  	Chistream.open("data/ave_Chi.dat", ios:: app);
  	if (!Chistream.is_open()) {
		cerr << "Error: unable to open ave_Chi.dat" << endl;
		exit(1);
		}
	
	Ustream << ave_U/double(nspin) << "\t" << error_U/double(nspin) << endl;
  	Mstream << ave_M/double(nspin) << "\t" << error_M/double(nspin) << endl;
  	Cstream << ave_C/double(nspin) << "\t" << error_C/double(nspin) << endl;
  	Chistream << ave_Chi/double(nspin) << "\t" << error_Chi/double(nspin) << endl;	
  	
  	Ustream.close();
  	Mstream.close();
  	Cstream.close();
  	Chistream.close();
  	
  	//put to zero in order to restart in a new block
  	block_U = 0.;
  	block_U2 = 0.;
  	block_M = 0.;
  	block_M2 = 0.;
  	step_run = 0;
  	
  	if(block_run == nblk){//when the blocks are finished write the results to file
  		const int wd = 12;
  		ofstream outfile;
  		outfile.open("data/results.dat", ios :: app);
  		if (!outfile.is_open()) {
			cerr << "Error: unable to open results.dat" << endl;
			exit(1);
		}
		outfile << temp << "\t" << setw(wd) << h << "\t" << setw(wd) << double(accepted)/double(attempted);
		outfile << "\t" << setw(wd) <<  ave_U/double(nspin) << "\t" << setw(wd) << error_U/double(nspin);
		outfile << "\t" << setw(wd) << ave_M/double(nspin) << "\t" << setw(wd) << error_M/double(nspin);
		outfile << "\t" << setw(wd) << ave_C/double(nspin) << "\t" << setw(wd) << error_C/double(nspin);
		outfile << "\t" << setw(wd) << ave_Chi/double(nspin) << "\t" << setw(wd) << error_Chi/double(nspin) << endl;
		outfile.close();
	}
}
	

void  Monte_Carlo_ISING_1D :: ConfFinal(void){
  	ofstream WriteConf;

  	cout << "Print final configuration to file config.final " << endl << endl;
  	WriteConf.open("config.final");
  	if (!WriteConf.is_open()) {
		cerr << "Error: unable to open config.final" << endl;
		exit(1);
	}
  	for (int i=0; i<nspin; ++i){
    		WriteConf << s[i] << endl;
  	}
  	WriteConf.close();

  	Rand.SaveSeed();
}

int Monte_Carlo_ISING_1D :: GetNumBlocks(){ return nblk; }

int Monte_Carlo_ISING_1D :: GetNumStepsPerBlock(){ return nstep; }


int Monte_Carlo_ISING_1D :: Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Monte_Carlo_ISING_1D :: Error(double sum, double sum2, int iblk)
{//error for the blocking method
    if(iblk==1 or iblk==0) return 0.0;
    else return sqrt(abs((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1)));
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
