#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

int main ()
{

	const int M_throws = 1000000;//number of throws
	const int N_blocks = 100;//number of blocks 
	int L = M_throws/N_blocks; //throws per block
	
	//Create a random number generator 
	Random Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in");
	
	//Save information about the numbers of our generation on a file
	ofstream outfile("data/info_ex2.1.dat");
	if (!outfile.is_open()) {
		cerr << "Error: unable to open info_ex2.1.dat" << endl;
		exit(1);
	}
	outfile << "Number of blocks: " << N_blocks << endl << "Number of throws: " << M_throws << endl << "Number of throws per block: " << L << endl;
	outfile.close();
	
	//open file to save the generated data 
	outfile.open("data/ex2.1.dat");
	if(!outfile.is_open()){
		cerr << "Error: unable to open ex2.1.dat" << endl;
		exit(1);
	}
	
	for (int i =0; i < N_blocks; i++){
		//accumulators
		double sum = 0., sum_sigma = 0.;
		double sum_imp = 0, sum_sigma_imp = 0.;// where imp stands for "importance"
		
		for (int j = 0; j < L; j++){
		
			//Generate a random number in [0,1)
			//These lines refer to the simple sampling
			double x = Rand.Rannyu();
			double f_x = M_PI/2 * cos(M_PI/2 * x);
			
			//Using the same number generated, we make a distribution following another probability density: IMPORTANCE SAMPLING
			double x_imp = 1 - sqrt(1 - x);
			double f_x_imp = M_PI/2 * cos(M_PI/2 * x_imp)/(2 - 2 * x_imp);
			//accumulate the values of the integrating function
			//uniform sampling 			
			sum += f_x;
			sum_sigma += (f_x - 1) * (f_x -1);
			//importance sampling 
			sum_imp += f_x_imp;
			sum_sigma_imp += (f_x_imp - 1.) * (f_x_imp - 1.);
		}
		//compute the averages and print the results
		double ave = sum/L;
		double ave_sigma= sum_sigma/L;
		double ave_imp = sum_imp/L;
		double ave_sigma_imp = sum_sigma_imp/L;
		
		outfile << ave << "\t" << ave_sigma << "\t" << ave_imp << "\t" << ave_sigma_imp << endl;
		
	}
	
	outfile.close();
	
	Rand.SaveSeed();
	
	return 0;
}
