#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

int main ()
{
	const int N = 10000; //number of throws
	int N_sums[4] = {1, 2, 10, 100}; //number of dices thrown together 
	const char printing[4] = {'\t', '\t', '\t', '\n'}; //vector of printing spacings 
	
	//Create a random number generator 
	Random Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in");
	
	//Save information about the numbers of our generation nd the number of dices on a file
	ofstream outfile("data/info_ex1.2.dat");
	
	if (!outfile.is_open()) {
		cerr << "Error: unable to open info_ex1.2.dat" << endl;
		exit(1);
	}
	outfile << "Number of throws: " << N << endl;
	outfile << "Number of dices summed up: ";
	for (int i = 0; i < 4; i++)
		outfile << N_sums[i] << printing[i];
	outfile.close();
	
	//open file to save the generated throws of dices
	ofstream outfile0("data/ex1.2UNIF.dat");
	ofstream outfile1("data/ex1.2EXP.dat");
	ofstream outfile2("data/ex1.2LOR.dat");
	if(!outfile0.is_open()){
		cerr << "Error: unable to open ex1.2UNIF.dat" << endl;
		exit(1);
	}
	if(!outfile1.is_open()){
		cerr << "Error: unable to open ex1.2EXP.dat" << endl;
		exit(1);
	}
	if(!outfile2.is_open()){
		cerr << "Error: unable to open ex1.2LOR.dat" << endl;
		exit(1);
	}
	
	for (int i = 0; i < N; i++){
		for(int j = 0; j < 4; j++){ //for each throw we do a cycle for each kind of throw: with 1, 2, 10 or 100 dices
			double sum = 0;
			double sum_exp = 0; 
			double sum_lorentz = 0;
			for (int k = 0; k < N_sums[j]; k++){
				sum += (int) Rand.Rannyu(1., 7.);// max is 7 because we take the integer part and so we cannot obtain 7, but we can obtain 6
				sum_exp += Rand.Exponential(1.);
				sum_lorentz += Rand.Lorentzian(0., 1.);
			}
			outfile0 << sum/N_sums[j] << printing[j];
			outfile1 << sum_exp/N_sums[j] << printing[j];
			outfile2 << sum_lorentz/N_sums[j] << printing[j];
		}
	}
	
	outfile0.close();
	outfile1.close();
	outfile2.close();
	
	Rand.SaveSeed();
	
	return 0;
}

			
			
		
	
