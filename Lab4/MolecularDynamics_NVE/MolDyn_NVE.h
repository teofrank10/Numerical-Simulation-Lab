/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __MolDyn_NVE__
#define __MolDyn_NVE__

#include <cstdlib>
#include <vector>
#include "../../rnumgen/random.h"

using namespace std;

class MolDyn_NVE{

private:
	//parameters, observables
	double temp, vol, rho, box, rcut, delta;
	int npart, nstep, iprint;
	int iblock, n_blocks; //useful for ex4.2 and 4.3 for the blocking method. in the first exercise are unused
	bool restart, print;
	int nbins= 100; //only useful in the ex7.3: computation of the radial distribution function
	double bin_size;// ""
	double vtail,ptail; //"" : tail corrections due to the cutoff
	
	//pigreco
	const double pi=3.1415927;
	
	//configuration
	static const int m_part=108;
	double x[m_part],y[m_part],z[m_part],x_old[m_part],y_old[m_part],z_old[m_part];
	double vx[m_part],vy[m_part],vz[m_part];
	
	double blockEpot, blockEkin, blockEtot, blockTemp, blockVir;//useful for ex4.2 and 4.3 for the blocking method. in the first exercise are unused
	double aveEpot, aveEkin, aveEtot, aveTemp, aveVir; //" "
	double ave2Epot, ave2Ekin, ave2Etot, ave2Temp, ave2Vir; //" "
	
	void Inizialization();
	void Restart();
	double Pbc(double);
	double Force(int, int);
	
	//vector for the calculation of g(r)
	vector<int> histogram;
	vector<double> g_block;
	vector<double> ave_g;
	vector<double> ave2_g;
	vector<double> g;
	vector<double> g_err;
	//only useful in the ex7.3: computation of the radial distribution function
	
protected:

public:

	//constructor
	MolDyn_NVE();
	//destructor
	~MolDyn_NVE();
	//methods: see the .cpp file for explanations
	void Input(void);
	void Move(void);
	void ConfFinal(void);
	void SaveOldConfig(void);
	void ConfXYZ(int);
	void Measure(void);
	int Steps(void);
	int Prints(void);
	void Blocking(void);//useful for ex4.2 and 4.3 for the blocking method. in the first exercise is unused
	
};

#endif //__MolDyn_NVE__


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
