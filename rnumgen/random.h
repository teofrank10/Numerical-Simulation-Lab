/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

#include <string>
#include <cmath>

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;
  double x, y, z, alpha; //coordinates for Metropolis algorithm and acceptance probability.

protected:

public:
  // constructors
  Random(std::string primes, std::string seed);
  // destructor 
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SetRandom(int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
  double Exponential(double lambda);
  double Lorentzian(double mean, double gamma);
  void Metropolis_Unif(double pos[], double delta, double (*p)(double[], int),int dim);
  void Metropolis_Gauss(double pos[], double sigma, double (*p)(double[], int), int dim);
  void Metropolis_Unif_1D(double pos, double delta, double (*p)(double, double, double), double mu, double sigma);
  double Metropolis_GetX();
  double Metropolis_GetY();
  double Metropolis_GetZ();
  double Metropolis_GetAlpha();
};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
