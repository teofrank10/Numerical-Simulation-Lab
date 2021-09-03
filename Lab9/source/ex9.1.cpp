#include "Genetic_Algorithm.h"

using namespace std;

int main(){

	Genetic_Algorithm GA;
	
	GA.Input();
	
	//create the population of routes
	GA.Create_Population();
	
	
	//we create the population of n_population routes composed by n_cities cities
	//we sort the population looking to the length of its routes: the population vector is filled from the shortest to the longest ruote.
	
	GA.Print_best_length(0);
	
	for(int i = 0; i < GA.Get_n_generations(); i++){
	
		GA.New_Generation(GA.Get_p_crossover(), GA.Get_p_mutation());

		if((i+1)%100 == 0){
			cout << endl << " Generation number " << i+1 << "    min length: " << GA.population[0][GA.Get_n_cities()] << endl;
			
		}
		//save to file the route with the minimum length
		GA.Print_best_length(i+1);
		
		}
		//save to file the coordinates of the cities visited in the best route
		GA.Print_best_route();
	
	
	return 0;
	
}
		
		
