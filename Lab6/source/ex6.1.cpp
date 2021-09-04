#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){
 	
 	Monte_Carlo_ISING_1D MC_ISING;
 	
	MC_ISING.Input(); //Inizialization
	
  	for(int iblk=1; iblk <= MC_ISING.GetNumBlocks(); iblk++){ //Simulation
    		for(int istep=1; istep <= MC_ISING.GetNumStepsPerBlock(); istep++){
      				MC_ISING.Move();
      				MC_ISING.Measure();
      				
    		}
    		MC_ISING.MeasureBlock();  //Print results for current block
  		}
  		MC_ISING.ConfFinal(); //Write final configuration
  	
  	return 0;

}
