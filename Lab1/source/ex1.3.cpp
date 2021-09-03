#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

int main ()
{
	const int N_throws = 10000;
	const int N_blocks = 100;
	const double L = 0.6; //the length of the needle 
	const double d = 1.; //the distance between two lines 
	
	//Create a random number generator 
	Random Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in");
	
	//Save information about the numbers of blocks and the number of throws in info file 
	ofstream outfile("data/info_ex1.3.dat");
	
	if (!outfile.is_open()) {
		cerr << "Error: unable to open info_ex1.3.dat" << endl;
		exit(1);
	}
	
	outfile << "Number of blocks: " << N_blocks << endl;
	outfile << "Number of throws per block: " << N_throws << endl;
	outfile.close();
	
	//open file to save the generated throws
	outfile.open("data/ex1.3.dat");
	if (!outfile.is_open()) {
		cerr << "Error: unable to open ex1.3.dat" << endl;
		exit(1);
	}
	
	for (int i = 0; i < N_blocks; i++){
		//define the counter of the times the needle touch a line
		int N_hit=0;
		
		for (int j =0; j < N_throws; j++){
			//generate one random point as the initial endpoint of the needle
			double a = Rand.Rannyu();
			
			//if the first endpoint fall on the line, increment the counter
			if (a == 0.) N_hit++;
			
			//Now control the inclination of the needle generating a point in an half of circle. y are always >= 0
			double x,y; 
			do{
				x = Rand.Rannyu(-1., 1);
				y = Rand.Rannyu();
			} while (x*x + y*y > 1.);
			
			//Using the generated 2D point, sample a theta angle in (0, pi) which represents the inclination of the needle
			double theta = acos(x/(sqrt(x*x + y*y)));
			//calculate the finale endpoint of the needle
			double b = a + L*sin(theta);
			
			//increment the counter if the final endpoint is beyond or on the line
			if (b > d or b == 1.) N_hit++;
		}
		//estimate pi and write it to file
		double pi_estim = 2 * L * N_throws/(N_hit * d);
		outfile << pi_estim << endl;
		
	}
	
	outfile.close();
	
	Rand.SaveSeed();
	
	return 0;
}
			
			
			
