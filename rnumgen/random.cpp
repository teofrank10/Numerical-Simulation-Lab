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
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "random.h"

using namespace std;

Random :: Random(string file_primes, string file_seed){

	int seed[4];
	int n1, n2;
	ifstream primes(file_primes);
	if (primes.is_open()){
		primes >> n1 >> n2;
	}
	else { 
		cerr << "Error: unable to open " << file_primes << endl;
		exit(1);
	}
	primes.close();
	
	ifstream start(file_seed);
	string property;
	if (start.is_open()){
		while ( !start.eof() ){
			start >> property;
			if (property == "RANDOMSEED"){
				start >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				SetRandom(seed, n1, n2);
			}
		}
		start.close();
	}
	else { 
		cerr << " Error: unable to open " << file_seed << endl;
		exit(1);
	}
}
		

Random :: ~Random(){ }

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}


void Random :: SetRandom(int num){
	int seed[4];
	int n1, n2;
	ifstream primes("../seed/Primes");
	if (primes.is_open()){
		for (int i = 0; i < num + 1; i++){
			primes >> n1 >> n2;
		}
	} else { 
		cerr << "Error: unable to open ../seed/Primes" << endl;
		exit(1);
	}
	primes.close();
	
	ifstream start("../seed/seed.in");
	string property;
	if (start.is_open()){
		while ( !start.eof() ){
			start >> property;
			if (property == "RANDOMSEED"){
				start >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				SetRandom(seed, n1, n2);
			}
		}
		start.close();
	} else { 
		cerr << " Error: unable to open seed.in " << endl;
		exit(1);
	}
}


double Random :: Exponential (double lambda){
	return - (1/lambda) * log(1 - Rannyu());
}

double Random :: Lorentzian (double mean, double gamma){
	return gamma * tan(M_PI * (Rannyu() - 0.5)) + mean;
}
//Metropolis algorithm used in ex5.1
void Random :: Metropolis_Unif(double pos[], double delta, double (*p)(double[], int), int dim){
	
	//set the coordinates of Metropolis equal to the actual position
	x = pos[0];
	y = pos[1];
	z = pos[2];
		
	double pos_fin[3];
	for (int i = 0; i < 3; i++){
		pos_fin[i] = Rannyu(pos[i] - delta, pos[i] + delta);
	}
	
	alpha = fmin(1, p(pos_fin, dim) / p(pos, dim));
	double r = Rannyu();
	if ( r < alpha) {
		x = pos_fin[0];
		y = pos_fin[1];
		z = pos_fin[2];
		return;
	} //otherwise we use the actual point
}
//Metropolis algorithm used in ex5.1
void Random :: Metropolis_Gauss(double pos[], double sigma, double (*p)(double[], int), int dim){

	//set the coordinates of Metropolis equal to the actual position
	x = pos[0];
	y = pos[1];
	z = pos[2];

	double pos_fin[3];
	for (int i = 0; i < 3; i++){
		pos_fin[i] = Gauss(pos[i], sigma);
	}
	
	alpha = fmin(1, p(pos_fin, dim) / p(pos, dim));
	double r = Rannyu();
	if ( r < alpha) {
		x = pos_fin[0];
		y = pos_fin[1];
		z = pos_fin[2];
		return;
	} //otherwise we use the actual point
}


void Random :: Metropolis_Unif_1D(double pos, double delta, double (*p)(double, double, double), double mu, double sigma){//useful for ex6.1 _ Ising Model
	
	//set the coordinates of Metropolis equal to the actual position
	x = pos;
		
	double pos_fin;
	
	pos_fin = Rannyu(pos - delta, pos + delta);
	
	
	alpha = fmin(1, p(pos_fin, mu, sigma) / p(pos, mu, sigma));
	double r = Rannyu();
	if ( r < alpha) {
		x = pos_fin;
		return;
	} //otherwise we use the actual point
}

double Random :: Metropolis_GetX(){ return x;}
double Random :: Metropolis_GetY(){ return y;}
double Random :: Metropolis_GetZ(){ return z;}
double Random :: Metropolis_GetAlpha(){ return alpha;}
		
		
	

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
