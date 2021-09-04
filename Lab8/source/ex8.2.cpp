#include <iostream>
#include <fstream>
#include <random.h>
#include <vector>
#include <string>
#include <cmath>

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

double D2(double x, double mu, double sigma){//second derivative 

	/*double exp1 = -pow((x - mu),2);
	double exp2 = -pow((x + mu),2);
	double denom = 2*pow(sigma,2);*/
	
	return -(pow(x,2) + pow(mu,2) - pow(sigma,2) - 2 * mu * x * tanh(mu * x / pow(sigma,2)))/(2*pow(sigma,4));
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
	
	ReadInput.open("input/input.dat");//read the input file
	
	if (!ReadInput.is_open()) {
		cerr << "Error: unable to open input.dat" << endl;
		exit(1);
	}


	cout << " Variational optimization of a single quantum particle Ground State " << endl;
	cout << " Variational Monte Carlo in 1D " << endl;
	cout << " External potential: V(x) = x^4 - 5/2*x^2 " << endl;
	cout << " Metropolis algorithm with transition probability uniformly distributed is used" << endl;
	
	double Input[9];
	double x = 0.;
	double delta = 0.;
	double mu = 0.;
	double sigma = 0.;
	int nsteps = 0;
	int nblocks = 0;
	double neqsteps = 0;
	bool optimization;
	bool paramsOpt;
	
	//fill the vector with Input parameters
	for (int i = 0; i < 9; i++){
	
		ReadInput >> Input[i];
	}
	
	x = Input[0];//initial point
	cout << " The simulation starts in x = " << x << endl;
	delta = Input[1];//delta Metropolis
	mu = Input[2];//variational parameter
	sigma = Input[3];//variational parameter
	nsteps = Input[4];
	nblocks = Input[5];
	neqsteps = Input[6];
	optimization = Input[7];//optimization of delta
	paramsOpt = Input[8];
	int L = nsteps / nblocks; //number of throws per block
	
	ReadInput.close();
	
	string WriteInput []= {"\t", "\tReadInput >> startingPoint;", "\tReadInput >> delta;", "\tReadInput >> mu;", "\tReadInput >> sigma;", "\tReadInput >> nsteps;", "\tReadInput >> nblocks;", "\tReadInput >> neqsteps;", "\tReadInput >> optimization;", "\tReadInput >> paramsOpt;"};
  
  
	//Create a random number generator
	Random Rand (SEEDDIR "/Primes", SEEDDIR "/seed.in");
	
	//code for the optimization of delta, mu and sigma: 
	// delta is optimized in order to reach 50% of AR
	//mu and sigma are optimized in order to obtain the minimum of energy
		
	ofstream Ham, outfile;
	Ham.open("data/params_optimization.dat");
	if (!Ham.is_open()){
		cout << "Error: unable to open params_optimization.dat" << endl;
		exit(1);
	}
	outfile.open("data/acceptance.dat");
	if (!outfile.is_open()){
		cout << "Error: unable to open outfile.dat" << endl;
		exit(1);
	}
				
	if(paramsOpt){
	
		//double energies [50][50]; // mu from 0.5 to 0.8 and sigma from 0.5 to 0.8: steps of 0.02
		for(int i = 0; i < 15; i++){
			
			mu += 0.02;//increment mu
			sigma = 0.5;//inizialize sigma before every cycle
			
			for(int j = 0; j < 15; j++){
				
				sigma += 0.02;//increment sigma
				
				//for each couple (mu,sigma) we execute the program in ex8.1 and take the final estimated energy.
	
				if(optimization){//delta optimization
			
					//cout << " The program performs delta optimization and equilibrates the system " << endl;
					int index = 0;
					double alpha_ave = 0.;
					do{
					
					double correction = 0.00007;
					
					
					//for (int index = 0; index < neqsteps; index++){
						Rand.Metropolis_Unif_1D(x, delta, psi2_t, mu, sigma);
						//calculate average of alpha
						alpha_ave = index/double(index+1)*alpha_ave + 1./double(index+1)*Rand.Metropolis_GetAlpha();
						
						//cout << index << endl;
						
						//correction of alpha
						if(alpha_ave < 0.50){
	
							delta -= correction;
		
						} else if( alpha_ave > 0.50){
		
							delta += correction;
						}
			
						//cout << i << "   " << j << "   " << index << "   " <<  "ACCEPTANCE RATE " << alpha_ave <<  endl;
						index++;
					} while(abs(0.50 - alpha_ave) > 0.005);//set the condition of correction
					
					//cout << delta << endl;
					
					//outfile << delta << "\t" << "ACC RATE: " << alpha_ave << endl;
		
				}	
				
				
					
				//Equilibration
				for (int k = 0; k < neqsteps; k++){
					Rand.Metropolis_Unif_1D(x, delta, psi2_t, mu, sigma);
		
					x = Rand.Metropolis_GetX();
				}
				
				//accumulators
				double ave = 0;
				double ave2 = 0;
	
				double H = 0;
	
				for(int l = 0; l < nblocks; l++ ){//computation of the energy
					double block_ave = 0;
		
					for (int m = 0; m < L; m++){
			
						Rand.Metropolis_Unif_1D(x, delta, psi2_t, mu, sigma);
						x = Rand.Metropolis_GetX();
						H = Kinetic(x, mu, sigma) + Potential (x);
						
		
						block_ave += H;
			
					}
					
					//blocking method
					block_ave /= L;
					ave = l/double(l+1)*ave + 1./double(l+1)*block_ave;
					ave2 = l/double(l+1)*ave2 + 1./double(l+1)*block_ave*block_ave;
					double err = error(ave,ave2,l);
					
					if (l == nblocks - 1){//print to file
						Ham << mu << "\t" << sigma << "\t" << ave << "\t" << err << "\t" << delta << endl;
		
					}
				}
			}
		}
	}
	
	outfile.close();
	
	Ham.close();
		
	return 0;
	
}
