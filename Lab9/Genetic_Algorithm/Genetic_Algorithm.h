
#ifndef __Genetic_Algorithm__
#define __Genetic_Algorithm__


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

class Genetic_Algorithm{

private:
	
	Random Rand;
	int n_cities, n_routes, n_generations;
	string distribution;//configuration of the cities
	
	vector<vector<double>> pos; //2d vector containing x and y coordinates
	
	//genetic parameters
	double p_crossover, p_mutation;
	double power;
	
	


protected:

public:

	Genetic_Algorithm();
	~Genetic_Algorithm();
	//matrix of all the routes
	vector<vector<double>> population;
	//methods
	void Input(void);
	void Cities(void);
	void Create_Population();
	void My_Sort(void);
	vector<double> Create_Route(int);
	void New_Generation(double, double); 
	vector<vector<double>> Crossover(int mum, int dad);
	int Selection(void);
	double Metric(vector<double> city_A, vector<double> city_B);
	double Length(vector<double> route);
	void Print_best_length(int step);
	void Print_best_route(void);
	int Get_n_generations(void){ return n_generations; };
	int Get_n_cities(void){ return n_cities; };
	int Get_n_routes(void){ return n_routes; };
	double Get_p_crossover(void){ return p_crossover; };
	double Get_p_mutation(void){ return p_mutation; };
	void Assign_Length(void);
};










#endif //__Genetic_Algorithm__
