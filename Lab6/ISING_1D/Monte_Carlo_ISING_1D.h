/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Monte_Carlo_Ising_1D__
#define __Monte_Carlo_Ising_1D__

#include <cmath>
#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
//Random numbers
#include "../../rnumgen/random.h"

class Monte_Carlo_ISING_1D{

private:
	Random Rand;
	// thermodynamical state
	int nspin;
	double beta,J,h;
	// simulation
	int nstep, n_eqstep, nblk, metro;
	bool restart;
	
	double accepted, attempted;
	
	double alpha;
	
	int step_run;
	int block_run;
	
	//configuration
	static const int m_spin=50;
	double s[m_spin];
	
	// averages
	double block_U, block_M;
	double block_U2, block_M2;
	double ave_U, ave_C, ave_Chi, ave_M;
	double ave2_U, ave2_C, ave2_Chi, ave2_M;
	
protected:

public:
	//temperature
	double temp;
	//constructor
	Monte_Carlo_ISING_1D();
	//destructor
	~Monte_Carlo_ISING_1D();
	//methods
	void Input(void);
	void Restart(void);
	void Inizialization(void);
	void Metropolis(int o);
	void Gibbs(int i);
	void Move(void);
	void ConfFinal(void);
	void Measure(void);
	void MeasureBlock(void);
	int GetNumBlocks();
	int GetNumStepsPerBlock();
	int Pbc(int);
	double Error(double,double,int);
	

};

#endif //__Monte_Carlo_Ising_1D__


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
