#include "Monte_Carlo_NVT.h"

using namespace std;

int main(){
 	
 	Monte_Carlo_NVT NVT;
 	
	NVT.Input(); //Inizialization
	
  	for(int iblk=0; iblk < NVT.GetNumBlocks(); ++iblk){ //Simulation
    		for(int istep=0; istep < NVT.GetNumStepsPerBlock(); ++istep){
      			NVT.Move();
      			NVT.Measure();
      				
    		}
    		NVT.Accumulate();   
    		NVT.Averages(); //Print results for current block
  	}
  	NVT.ConfFinal(); //Write final configuration
  	
  	return 0;
  	
}
