/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Monte_Carlo_NVT__
#define __Monte_Carlo_NVT__

#include <cmath>
#include <string>
#include <ostream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
//Random numbers
#include "../../rnumgen/random.h"

using namespace std;

class Monte_Carlo_NVT{

private:

	Random Rand;
	
	double vtail,ptail;
	double v_blocking, w_blocking; 
	double ave_v, ave2_v;
	double ave_w, ave2_w;
	double bin_size,nbins;

	// averages
	double accepted,attempted;

	// thermodynamical state
	int npart;
	double beta,temp,vol,rho,box,rcut;
	
	//configuration
	static const int m_part=108;
	double x[m_part],y[m_part],z[m_part];
	
	//vector for the calculation of g(r)
	vector<int> histogram;
	vector<double> g_block;
	vector<double> ave_g;
	vector<double> ave2_g;
	vector<double> g;
	vector<double> g_err;
	

	// simulation
	int nstep, neqstep, nblk , istep, iblock;
	double delta;
	bool restart;
	bool print;//boolean variable which allows to control the printing of instantaneous variable energy and pressure.

	//pigreco
	const double pi=3.1415927;

protected:

public:

	//Constructor
	Monte_Carlo_NVT();
	//desctructor
	~Monte_Carlo_NVT();
	//methods
	void Input(void);
	void Restart(void);
	void Initialization(void);
	void Move(void);
	double Boltzmann(double, double, double, int);
	void Measure(void);
	void Accumulate(void);
	void Averages(void);
	void Reset(void);
	void ConfFinal(void);
	void ConfXYZ(int);
	double Pbc(double);
	double Error(double,double,int);
	int GetNumBlocks();
	int GetNumStepsPerBlock();
	
};

#endif //__Monte_Carlo_NVT__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
