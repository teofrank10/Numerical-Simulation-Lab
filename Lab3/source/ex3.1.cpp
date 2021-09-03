#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

int main ()
{

	const int M_throws = 10000;//number of throws for each block
	const int N_blocks = 100;
	const int L_steps = 100; //number of steps for each throw
	const double S0 = 100; // asset price at t=0
	const double T = 1; //delivery time
	const double K = 100; //strike price 
	const double r = 0.1; //risk-free interest rate
	const double sigma = 0.25; //volatility;	
	
	//Create a random number generator 
	Random Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in");
	
	//Save information about the numbers of our generation on a file
	ofstream outfile("data/info_ex3.1.dat");
	if (!outfile.is_open()) {
		cerr << "Error: unable to open info_ex3.1.dat" << endl;
		exit(1);
	}
	outfile << "Number of blocks: " << N_blocks << endl << "Number of throws per block: " << M_throws << endl << "Number of steps per throw: " << L_steps << endl;
	outfile.close();
	
	//open file to save the generated data 
	outfile.open("data/CALL_ex3.1.dat");
	if(!outfile.is_open()){
		cerr << "Error: unable to open CALL_ex3.1.dat" << endl;
		exit(1);
	}
	
	//============CALL PRICING==============//
	
	for (int i = 0; i < N_blocks; i++){
		double sum_dir= 0.;
		double sum_disc= 0.;
		for (int j = 0; j < M_throws; j++){
			//direct sampling:
			double W = Rand.Gauss(0., T); //extract a random number for W
			double ST_dir = S0 * exp((r - 0.5 * sigma * sigma) * T - sigma * W * sqrt(T)); //use the asset price formula 
			ST_dir -= K; // the profit for a call il max(0, ST - K). We have to work with positive profit:
			if(ST_dir > 0){ 
				sum_dir += ST_dir * exp(-r * T); // in order to find the intial value we have to make a discount with the exponential
			}
			//discrete sampling from path:
			double ST_disc = S0;
			for (int k = 0; k < L_steps; k++){
				double Z = Rand.Gauss(0., T);
				ST_disc *= exp((r - 0.5 * sigma * sigma) * T/L_steps + sigma * Z * sqrt(T/L_steps));
			}
			ST_disc -= K;
			if(ST_disc > 0){ 
				sum_disc += ST_disc * exp(-r * T); // in order to find the intial value we have to make a discount with the exponential
			}
		}
		double ave_dir = sum_dir/M_throws;
		double ave_disc = sum_disc/M_throws;
		
		outfile << ave_dir << "\t" << ave_disc << endl;
	}
	outfile.close();
	
	//open file to save the generated data 
	outfile.open("data/PUT_ex3.1.dat");
	if(!outfile.is_open()){
		cerr << "Error: unable to open PUT_ex3.1.dat" << endl;
		exit(1);
	}
	
	//============PUT PRICING==============//
	
	for (int i = 0; i < N_blocks; i++){
		double sum_dir = 0.;
		double sum_disc = 0.;
		for (int j = 0; j < M_throws; j++){
			//direct sampling
			double W = Rand.Gauss(0., T);
			double ST = S0 * exp((r - 0.5 * sigma * sigma) * T - sigma * W * sqrt(T));
			ST = K - ST; // the profit for a call il max(0, K - ST). We have to work with positive profit:
			
			if(ST > 0){ 
				sum_dir += ST * exp(-r * T); // in order to find the intial value we have to make a discount with the exponential
			}
			
			//discrete sampling from path:
			double ST_disc = S0;
			for (int k = 0; k < L_steps; k++){
				double Z = Rand.Gauss(0., T);
				ST_disc *= exp((r - 0.5 * sigma * sigma) * T/L_steps + sigma * Z * sqrt(T/L_steps));
			}
			ST_disc = K - ST_disc;
			if(ST_disc > 0){ 
				sum_disc += ST_disc * exp(-r * T); // in order to find the intial value we have to make a discount with the exponential
			}
		}
		double ave_dir = sum_dir/M_throws;
		double ave_disc = sum_disc/M_throws;
		
		outfile << ave_dir << "\t" << ave_disc << endl;
	}
	outfile.close();
	
	Rand.SaveSeed();
	
	return 0;
}
		
	
	
