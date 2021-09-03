#include "Genetic_Algorithm.h"

using namespace std;

int main(int argc, char *argv[]){

	

	MPI_Init(&argc, &argv);//initialize the MPI 
	
	Genetic_Algorithm GA;
	
	
	MPI_Comm_size(MPI_COMM_WORLD, &GA.Size);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &GA.Rank);
	
	MPI_Barrier(MPI_COMM_WORLD);//synchronization
	//start the chrono	
	double t_start = MPI_Wtime();
	
	GA.Input();
	
	
	//create the population of routes
	GA.Create_Population();
	
	
	//we create the population of n_population routes composed by n_cities cities
	//we sort the population looking to the length of its routes: the population vector is filled from the shortest to the longest ruote.

	
	GA.Print_best_length(0);
	
	
	
	for(int i = 0; i < GA.Get_n_generations(); i++){
	
		GA.New_Generation(GA.Get_p_crossover(), GA.Get_p_mutation());

		if((i+1)%1100== 0){//print the status of the program
			cout << endl << "Node number: " << GA.Rank << "  is calculating...   Generation number " << i+1 << "    min length: " << GA.population[0][GA.Get_n_cities()] << endl;	
		}
		
		if((i+1)%GA.Get_n_generations() == 0){//print the end
			cout << endl << "Node number: " << GA.Rank << "  finished the calculations.   Generation number " << i+1 << "    min length: " << GA.population[0][GA.Get_n_cities()] << endl;	
		}
		
		if((i+1)%GA.Get_n_migrations() == 0){//migrations
		
			GA.Migration();

		}
		//save to file the route with the minimum length
		GA.Print_best_length(i+1);
		
	}
		//save to file the coordinates of the cities visited in the best route
	GA.Print_best_route();
		
	double t_final = MPI_Wtime();//stop time 
	double delta_t = t_final-t_start;
	//print the time 
	cout << endl << "Rank " << GA.Rank << ".  Time " << delta_t << endl;
	
	MPI_Finalize();
	
	
	return 0;
	
}
		
		
