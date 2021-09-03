#ifndef __Simulated_Annealing__
#define __Simulated_Annealing__


#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <utility>
#include <numeric>
//Random numbers
#include "../rnumgen/random.h"

using namespace std;

class Simulated_Annealing{

private:
	
	Random Rand;
	int n_cities;
	string distribution;//configuration of the cities
	vector<double> route;
	vector<double> best_route;
	double best_length;
	vector<vector<double>> pos; //2d vector containing x and y coordinates
	
	//simulation
	int SA_steps = 0;
	int n_temp_steps, n_metrop_steps;
	double temp, temp_decay_rate;
	int attempted, accepted;
	double beta;
	double alpha;


protected:

public:

	Simulated_Annealing();
	~Simulated_Annealing();
	//methods
	void Input(void);
	void Cities(void);
	vector<double> Create_Route(int);
	vector<double>Annealation(vector<double> route);
	void Metropolis(void);
	void Reset(void);
	void Find_best_length(void);
	double Length(vector<double> route);
	void Print_best_length(int step);
	void Print_best_route(void);
	int Get_Accepted(void){ return accepted;}
	int Get_Attempted(void){ return attempted;}
	int Get_T_Steps(void){ return n_temp_steps;}
	int Get_Metrop_Steps(void){ return n_metrop_steps;}
	void Set_N_Metrop_Steps(int n){ n_metrop_steps += n;}
	double Get_Temp(void){ return temp;}
	void Set_Temp(void){ temp *= temp_decay_rate;}
	double Get_Best_length(void){ return best_length;}
	
};










#endif //__Simulated_Annealing__
