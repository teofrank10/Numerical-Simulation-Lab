#include <iostream>
#include <fstream>
#include <random.h>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std; 

//==========================Wave function=========================//
double psi_t(double x, double mu, double sigma){
	double exp1 = -pow((x - mu),2);
	double exp2 = -pow((x + mu),2);
	double denom = 2*pow(sigma,2);
	
	return exp(exp1/denom)+exp(exp2/denom);
}

double psi2_t(double x, double mu, double sigma){//probability density to be sampled
	
	return pow(psi_t(x,mu,sigma),2);
}

double D2_psi_t(double x, double mu, double sigma){//second derivative 

	double exp1 = -pow((x - mu),2);
	double exp2 = -pow((x + mu),2);
	double denom = 2*pow(sigma,2);
	
	return -psi_t(x,mu,sigma)/pow(sigma,2) - exp1/pow(sigma,4)*exp(exp1/denom) - exp2/pow(sigma,4)*exp(exp2/denom);
}
	
//================================================================//
//==============================Energies==========================//
double Potential(double x){

	return pow(x,4) - 5. * x * x/2.;
}

double Kinetic(double x, double mu, double sigma){

	return - D2_psi_t(x,mu,sigma)/(2*psi_t(x,mu,sigma));
}

double error(double ave, double ave2, double n){

	if(n == 0){return 0;}else {
		return sqrt(abs(ave2-ave*ave)/n);
	}
}
//================================================================//


int main (){

	ifstream ReadInput;
	
	ReadInput.open("data/params_optimization.dat");//read the file which contains the values of the energy for different values of the parameters
	
	if (!ReadInput.is_open()) {
		cerr << "Error: unable to open params_optimization.dat" << endl;
		exit(1);
	}

	cout << " Estimation of the Ground State energy " << endl;
	cout << " External potential: V(x) = x^4 - 5/2*x^2 " << endl;
	cout << " Metropolis algorithm with transition probability uniformly distributed " << endl;
	cout << " The variational parameters and the Metropolis delta are already optimized  " << endl;
	
	
	string line;
	vector<double> energies;
	double err_energy[225];
	double val_mu[225];
	double val_sigma[225];
	double val_delta[225];
	
	
	
	
	int i = 0;
	while(getline(ReadInput, line)){//read one line from ReadInput
		
		istringstream iss(line); //access line as a stream 
		
		double en;
		//we need only the third line of the file 
		iss >> val_mu[i] >> val_sigma[i] >> en >> err_energy[i] >> val_delta[i];
		//fill the vector with energies
		energies.push_back(en);
		i++;

	}
	
	ReadInput.close();
	
	//find the index of the minimum of the energy.
	int MinimumIndex = min_element(energies.begin(),energies.end()) - energies.begin();
	
	
	cout << " The program runs with optimized parameters " << endl;
	cout << " OPTIMIZED MU = " << val_mu[MinimumIndex] << endl;
	cout << " OPTIMIZED SIGMA = " << val_sigma[MinimumIndex] << endl;
	cout << " OPTIMIZED DELTA = " << val_delta[MinimumIndex] << endl << endl;
		
	cout << " The minimum energy for the GS is " << energies[MinimumIndex] << " +- " << err_energy[MinimumIndex] << endl;
	
	double x = 0.; //starting point
	double mu = val_mu[MinimumIndex];
	double sigma = val_sigma[MinimumIndex];
	double delta = val_delta[MinimumIndex];
	
	int nsteps = 1e+6;//total number of steps 
	int nblocks = 100;//number of blocks
	int neqsteps = 500;//equilibrations steps
	int L = nsteps / nblocks; //number of throws per block
	
	//Create a random number generator
	Random Rand (SEEDDIR "/Primes", SEEDDIR "/seed.in");
	
	//we remind that  the variational parameters are already optimized.
	//so we can perform the equilibration and the estimatons required in the exercise
		
	double alpha_ave = 0.;
	
	ofstream Equilib;
	Equilib.open("data/equilibration.dat", ios::app);
	if (!Equilib.is_open()){
		cout << "Error: unable to open equilibration.dat" << endl;
		exit(1);
	}
	for (int i = 0; i < neqsteps; i++){
		Rand.Metropolis_Unif_1D(x, delta, psi2_t, mu, sigma);
		alpha_ave = i/double(i+1)*alpha_ave + 1./double(i+1)*Rand.Metropolis_GetAlpha();
		//double x_t = Rand.Metropolis_GetX();
		
		//update the vector pos:
		x = Rand.Metropolis_GetX();
		//double r_t = sqrt(pow(x_t,2)+pow(y_t,2)+pow(z_t,2));
		Equilib << x << endl;
	}
	cout << " ===Equilibration:============== " << endl;
	cout << " ACCEPTANCE RATE = " << alpha_ave << endl;
	
	//estimation of H
	double ave = 0;//initialize the accumulators
	double ave2 = 0;
	
	double H = 0;
	
	ofstream Ham, Positions;
	Ham.open("data/hamiltonian.dat");
	if (!Ham.is_open()){
		cout << "Error: unable to open hamiltonian.dat" << endl;
		exit(1);
	}
	
	Positions.open("data/positions.dat");
	if (!Positions.is_open()){
		cout << "Error: unable to open positions.dat" << endl;
		exit(1);
	}
	double AR = 0.;//initialize the acceptance rate
	for(int i = 0; i < nblocks; i++ ){
		double block_ave = 0;
		
		for (int j = 0; j < L; j++){
			
			Rand.Metropolis_Unif_1D(x, delta, psi2_t, mu, sigma);//use metropolis to get x
			AR = j/double(j+1)*AR + 1./double(j+1)*Rand.Metropolis_GetAlpha();//calculate the average of acceptance rate
			x = Rand.Metropolis_GetX();
			H = Kinetic(x, mu, sigma) + Potential (x);//compute the value of hamiltonian as the sum of kinetic part and potential part
			Positions << x << endl;
			block_ave += H;///accumulate for each block
			
		}
		cout << " ===Estimation:============== " << i <<  endl;
		cout << " ACCEPTANCE RATE = " << AR << endl;
			
		block_ave /= L;//make the average on the number of steps for each block
		ave = i/double(i+1)*ave + 1./double(i+1)*block_ave;//as always: blocking method
		ave2 = i/double(i+1)*ave2 + 1./double(i+1)*block_ave*block_ave;
		double err = error(ave,ave2,i);
		
		Ham << ave << "\t" << err << endl;
		
	}
	
	Positions.close();
	Ham.close();

	return 0;
}
