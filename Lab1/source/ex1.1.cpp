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
#include <cmath>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

int main ()
{

	const int M_throws = 1000000;
	const int N_blocks = 100;
	int L = M_throws/N_blocks; //throws per block
	
	//Create a random number generator 
	Random Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in");
	
	//Save information about the numbers of our generation on a file
	ofstream outfile("data/info_ex1.1.dat");
	if (!outfile.is_open()) {
		cerr << "Error: unable to open info_ex1.1.dat" << endl;
		exit(1);
	}
	outfile << "Number of blocks: " << N_blocks << endl << "Number of throws: " << M_throws << endl << "Number of throws per block: " << L << endl;
	outfile.close();
	
	//open file to save the generated data 
	outfile.open("data/ex1.1.dat");
	if(!outfile.is_open()){
		cerr << "Error: unable to open ex1.1.dat" << endl;
		exit(1);
	}
	
	//for cycle on the number of blocks
	for (int i = 0; i < N_blocks; i++){
		double sum = 0., sum_sigma = 0.;
		//create a vector which represents the division of [0,1) in N_blocks bins
		int n[N_blocks] = {0}; //this vector will contain the number of effective observation of numbers in a certain interval
	
		for (int j = 0; j < L ; j++){
		
			//generate a random number r in [0, 1);
			double r = Rand.Rannyu();
			
			//define the total sum of random number and the total sum of square difference
			sum += r;
			sum_sigma += (r - 0.5) * (r - 0.5);
			
			//increment by one the elements contained in the bin between [i/N_blocks, (i+1)/N_blocks)
			n[(int) (r*N_blocks)]++;
		}
	
		//calculate the average amd the average of square difference
		double ave = sum/L;
		double ave_sigma = sum_sigma/L;
		
		double chi2 = 0.;
		for (int j = 0; j < N_blocks; j++){
			//use the definition of chi square to compare the actual distribution with the expected one
			chi2 += (n[j]-L/N_blocks) * (n[j]-L/N_blocks);
		}
		chi2 *= ((double) N_blocks/L);
		//write data to file 
		outfile << ave << "\t" << ave_sigma << "\t" << chi2 << endl;
	}
	
	outfile.close();
	
	Rand.SaveSeed();
		
	return 0;
}
