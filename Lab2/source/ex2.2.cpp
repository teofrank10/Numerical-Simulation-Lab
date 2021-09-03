#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;
//function useful to calculate dev std
double error(double ave, double ave2, double n){

	if(n == 0){return 0;}else{
		return sqrt((ave2-ave*ave)/n);
	}
}

int main (){

	const int M_walks = 10000;//number of walks 
	const int N_steps = 100;//number of steps for each walk
	const int N_blocks = 100;
	const int L= M_walks/N_blocks; //number of walks for each block
	
	//Create a random number generator 
	Random Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in");
	
	//Save information about the numbers of our generation on a file
	ofstream outfile("data/info_ex2.2.dat");
	if (!outfile.is_open()) {
		cerr << "Error: unable to open info_ex2.2.dat" << endl;
		exit(1);
	}
	outfile << "Number of blocks: " << N_blocks << endl << "Number of step per walk: " << N_steps << endl << "Number of walks per block: " << L << endl;
	outfile.close();
	
	//open file to save the generated data 
	outfile.open("data/ex2.2.dat");
	if(!outfile.is_open()){
		cerr << "Error: unable to open ex2.2.dat" << endl;
		exit(1);
	}
	
	
//==================================================================//
//this first part of the code does not satisfy the requests of the exercise//	
	for (int i = 0; i < N_blocks; i++){
		//define the accumulator for the square module of r
		double sum_r2 = 0.;
		double sum_r2_cont= 0.;
		
		for (int j = 0; j < L; j++){
		
			//define the 3D vector containing the components of r
			int r_disc[3] = {0};
			double r_cont[3] = {0.}; // is the discrete case: components are integers
			for (int k = 0; k < N_steps; k++){
			
				//Discrete case: select one axis and the verse
				int direction = (int) Rand.Rannyu(0, 3);
				double verse = Rand.Rannyu(-1, 1);
				
				if ( verse < 0.){
					r_disc[direction]--;
				} else { r_disc[direction]++;}
				
				
				//continuum case: spherical coordinates
				
				double phi = Rand.Rannyu(0., 2 * M_PI);
				double x = Rand.Rannyu(-1, 1);
				double theta = acos(x);
				
				r_cont[0] += sin(theta) * cos(phi);
				r_cont[1] += sin(theta) * sin(phi);
				r_cont[2] += cos(theta);
			
				
				
			}
			
			sum_r2 += r_disc[0] * r_disc[0] + r_disc[1] * r_disc[1] + r_disc[2] * r_disc[2];
			sum_r2_cont += r_cont[0] * r_cont[0] + r_cont[1] * r_cont[1] + r_cont[2] * r_cont[2];
		}
		
		double ave_disc = sum_r2/L;
		double ave_cont = sum_r2_cont/L;
		outfile << sqrt(ave_disc) << "\t" << sqrt(ave_cont) << endl;
	}
	
	outfile.close();
//==================================================================//
	
	double r2_disc[N_blocks] = {0.};
	double r4_disc[N_blocks] = {0.};
	double r2_cont[N_blocks] = {0.};
	double r4_cont[N_blocks] = {0.};
	
	
	//inizialize a matrix of the random walks with their coordinates
	double RW_disc[M_walks][3];
	double RW_cont[M_walks][3]; 
	
	//fora each RW set the initial point to the origin
	for (int i = 0; i < M_walks; i++){
		RW_disc[i][0] = 0.;
		RW_disc[i][1] = 0.;
		RW_disc[i][2] = 0.;
		
		RW_cont[i][0] = 0.;
		RW_cont[i][1] = 0.;
		RW_cont[i][2] = 0.;
	}
		
	
	outfile.open("data/ex2.2_new.dat");
	if(!outfile.is_open()){
		cerr << "Error: unable to open ex2.2_new.dat" << endl;
		exit(1);
	}
	
	//write on the file that the origin is the point occupied at zero steps
	outfile << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;
	for (int i = 0; i < N_steps-1; i++){ //nsteps-1 because we already wrote the first one 
		
		//accumulators for each step of the block
		double accu_disc = 0.;
		double accu2_disc = 0.;
		double accu_cont = 0.;
		double accu2_cont = 0.;
		for (int j = 0; j < N_blocks; j++){	
			//define the accumulator for the square module of r
			double sum_r2_disc = 0.;
			double sum_r2_cont= 0.;	
			
			for (int k = L*j; k < L*(j+1); k++){
		
				//Discrete case: select one axis and the verse
				int direction = (int) Rand.Rannyu(0, 3);
				double verse = Rand.Rannyu(-1, 1);
			
				if ( verse < 0.){
					RW_disc[k][direction]--;
				} else { RW_disc[k][direction]++;}
				
				//continuum case: spherical coordinates
				
				double phi = Rand.Rannyu(0., 2 * M_PI);
				double x = Rand.Rannyu(-1, 1);
				double theta = acos(x);
				
				RW_cont[k][0] += sin(theta) * cos(phi);
				RW_cont[k][1] += sin(theta) * sin(phi);
				RW_cont[k][2] += cos(theta);
			
				sum_r2_disc += RW_disc[k][0]*RW_disc[k][0] + RW_disc[k][1]*RW_disc[k][1] + RW_disc[k][2]*RW_disc[k][2];
				
				sum_r2_cont += RW_cont[k][0]*RW_cont[k][0] + RW_cont[k][1]*RW_cont[k][1] + RW_cont[k][2]*RW_cont[k][2];
			
			}
			
			//perform the average on the walks in each block
			sum_r2_disc /= L;
			sum_r2_cont /= L;
			
			//discrete case and continuous: fill the vector with the averages and the squared averages of each block
			r2_disc[j] = sum_r2_disc;
			r4_disc[j] = sum_r2_disc * sum_r2_disc;
			accu_disc += r2_disc[j];
			accu2_disc += r4_disc[j];
			
			r2_cont[j] = sum_r2_cont;
			r4_cont[j] = sum_r2_cont * sum_r2_cont;
			accu_cont += r2_cont[j];
			accu2_cont += r4_cont[j];
		}
		
		//compute the averages of the squared distance from the origin as functions of the step 
		double ave_disc = accu_disc/N_blocks;
		double ave_cont = accu_cont/N_blocks;
		
		double ave2_disc = accu2_disc / N_blocks;
		double ave2_cont = accu2_cont / N_blocks;
		 
		//compute the statistical errors
		double e_disc = error(ave_disc,ave2_disc, N_blocks);
		double e_cont = error(ave_cont,ave2_cont, N_blocks);
		
		outfile << sqrt(ave_disc) << "\t" << e_disc/(2*sqrt(ave_disc))  << "\t" << sqrt(ave_cont) << "\t" << e_cont/(2*sqrt(ave_cont)) << endl;
	}
	
	outfile.close();
	
	Rand.SaveSeed();
	
	
	
	return 0;
}

