/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     
#include <iostream>     
#include <fstream>      
#include <cmath>        
#include "../MolecularDynamics_NVE/MolDyn_NVE.h"

using namespace std;

int main(){ 

	MolDyn_NVE NVE;//Construct the NVE object which has all the methods
	
	NVE.Input(); //Inizialization of the system or restarting of the simulation 
	
	int nconf = 1; //put the number of configurations equal to 1
	
	
	//Equilibrate system 
	
	for (int i = 0; i < NVE.Steps(); i++) {//this is a cycle on the steps of integration of the Hamilton equations	
		if(i%NVE.Prints() == 0) {
     			cout << "Number of time-steps: " << i << endl;//print command to show the point reached during the run
     		}
     			
     		if(i%10 == 0){//every 10 steps we take a measure of the quantities which are E_pot E_kin E_tot Temp
        		NVE.Measure(); //Properties measurement
				
			nconf += 1; //increment the number of configurations. This number is useful in order to calculate the mean of the different observables.. To do the means we have to know the number of configurations on which we did the measure
		}
		
		NVE.Move();
	}
	
	NVE.SaveOldConfig();
	NVE.ConfFinal();         /*Write final configuration to restart. Prints the final configuration not in the xyz format but in the format used by this code, i.e. the 3 coordinates for each d.o.f.*/
  	/*In both configurations, initial and final, the unity of measure of the lenght is the lenght of the side called "box". This is useful because independently on the density of the system we can use whatever generated configuration, because starting from the density of the system I can find the other observables multiplying the density by the lenght of the side I'm using.*/	
	

  	return 0;
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
