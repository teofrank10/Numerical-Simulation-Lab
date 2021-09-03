#include "Simulated_Annealing.h"

using namespace std;

int main(){

	Simulated_Annealing SA;
	
	SA.Input();
	
	
	SA.Print_best_length(0);
	
	//we create the population of n_population routes composed by n_cities cities
	//we sort the population looking to the length of its routes: the population vector is filled from the shortest to the longest ruote.
	
	for (int i = 0; i < SA.Get_T_Steps(); i++){
		SA.Reset();//reset the accumulators and increase the number of Simulated annealing steps 
		
		for (int j = 0; j < SA.Get_Metrop_Steps(); j++){
		
			SA.Metropolis();//use metropolis algorithm
		}
		
		SA.Find_best_length();
		SA.Print_best_length(i+1);
		
		if((i+1)%30==0){
			cout << "Temperature step number " << i+1 << ", T = " << SA.Get_Temp() << endl << " min length = " << SA.Get_Best_length() << endl;
			
			double AR = 0;
			AR = (double) SA.Get_Accepted()/ (double) SA.Get_Attempted();// compute acceptance rate
			
			cout << endl << "Acceptance rate = " << AR << endl;
		}
		SA.Set_Temp();//use the decay rate
		SA.Set_N_Metrop_Steps(100);//increase the number of metropolis steps
	}
		//save to file the coordinates of the cities visited in the best route
		SA.Print_best_route();
	
	
	return 0;
	
}
		
		
