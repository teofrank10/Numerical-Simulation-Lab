#include "Simulated_Annealing.h"

using namespace std;

Simulated_Annealing :: Simulated_Annealing(): Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in"){ }

Simulated_Annealing :: ~Simulated_Annealing(){ }
 

void Simulated_Annealing :: Input(void){//initialization

	ifstream ReadInput;
	ReadInput.open("Simulated_Annealing/input.dat");
	if (!ReadInput.is_open()) {
		cerr << "Error: unable to open input.dat" << endl;
		exit(1);
	}

	cout << " Travelling Salesman Problem: solution with the Simulated Annealing algorithm"<< endl;
	
	ReadInput >> n_cities;
	ReadInput >> distribution;
	
	cout << "Number of cities : " << n_cities << endl;
	cout << "Configuration: " << distribution << endl;

	//set the configuration of the cities and print them
	Cities();
	
	ReadInput >> temp;
	ReadInput >> temp_decay_rate;
	ReadInput >> n_temp_steps;
	ReadInput >> n_metrop_steps;
	
	cout << "Initial temperature = " << temp << endl;
	cout << "Temperature decay rate = " << temp_decay_rate << endl;
	cout << "Number of temperature steps = " << n_temp_steps << endl;
	cout << "Number of Metropolis steps = " << n_metrop_steps << endl;
	
	ReadInput.close();
	
	//we set only the final value of number of steps. So, in order to have the initial number of steps, we do:
	double denom = n_metrop_steps - 100*(n_temp_steps -1);
	denom = n_metrop_steps/denom;
	n_metrop_steps = int(n_metrop_steps/denom);
	
	route = Create_Route(n_cities);  //set the route to the initial random route
	best_route = route;   //set the best route
	best_length = Length(route);  //and the best length
		
	
}

void Simulated_Annealing :: Cities(void){

	//set the configuration and print the result
	
	if(distribution == "circumference"){
		//generate a distribution of cities on a circumference of radius 10
		for(int i = 0; i < n_cities; i++){
		
			double phi = Rand.Rannyu() * 2 * M_PI;
			double x = 10 * cos(phi);
			double y = 10 * sin(phi); 
			pos.push_back({x,y});
		}
	}
	
	if(distribution == "square"){
		//generate a distribution of cities in a square of edge 10
		for(int i = 0; i < n_cities; i++){
		
			double x = 10 * Rand.Rannyu();
			double y = 10 * Rand.Rannyu();
			pos.push_back({x,y});
		}
	}
	
	
	//save the positions into file
	ofstream outfile("data/pos_cities.dat");
	
	for (int i = 0; i < n_cities; i++){
		outfile << pos[i][0] << "\t" << pos[i][1] << endl;
	}
	//the last point has to be the intial point 
	outfile << pos[0][0] << "\t" << pos[0][1] << endl;
	
	outfile.close();
	
}

vector<double> Simulated_Annealing :: Create_Route(int n_cities){

	vector<double> route;
	
	route.push_back(0.0);//we fix the first city to be the 0th for each route
	int r = (int) Rand.Rannyu(0,n_cities-1);
	route.push_back(1+r);//first visited city
	
	//condition for which all the cities are visited once 
	while( (int) route.size() != n_cities){
		int rr = 1 + (int) Rand.Rannyu(0, n_cities-1);
		//check if such value exists in the vector otherwise push it back
		if(find(route.begin(), route.end(), rr) == route.end()){
			
			route.push_back(rr);
		}
	}
	
	route.push_back(0.0); //here we save the length of the route: as the 33rd element of the vector

	return route;
}

vector<double> Simulated_Annealing :: Annealation(vector<double> route){

	vector<double> new_route = route;
	
	int index = (int) Rand.Rannyu(0,3);
	//here we use the same code of the mutations in the genetic algorithm
	
	if(index == 0){//permutation of a couple of cities 
		int n1 = (int)Rand.Rannyu(1,n_cities);
		int n2 = (int)Rand.Rannyu(1,n_cities);
		
		swap(new_route[n1], new_route[n2]);
	}
	
	if(index == 1){//Inversion of the vector
		int n1 = (int)Rand.Rannyu(1,n_cities);
		int n2 = (int)Rand.Rannyu(1,n_cities);
		
		if(n2 > n1){
			reverse(new_route.begin() + n1, new_route.begin() + n2);
				//cout << " mutation2.1 AVVENUTO " << endl;
		} else {
			reverse(new_route.begin() + n2, new_route.begin() + n1);
				//cout << " mutation2.2 AVVENUTO " << endl;
		}
	}
	
	if(index == 2){//permutation of adjacent
				
		vector<int> s;
		//we want to use rotate, so we have to random generate the first, the last and the middle iterators. Rotate will act rotating the vector until "middle" becomes the first element of the vector.
		s.push_back((int)Rand.Rannyu(1, n_cities));
		s.push_back((int)Rand.Rannyu(1, n_cities));	
		s.push_back((int)Rand.Rannyu(1, n_cities));
		sort(s.begin(), s.end());// we sort the vector s in order to have first, middle and last ordered
		
		rotate(new_route.begin()+s[0], new_route.begin()+s[1], new_route.begin()+s[2]);
	}
	
	return new_route;
			
}

void Simulated_Annealing :: Metropolis(void){//Metropolis algorithm

	vector<double> route_1 = Annealation(route);//new route obtained with annealing
	
	double delta = Length(route_1) - Length(route);
	double p = exp(-beta*delta);
	alpha = min(1.,p);
	
	double random = Rand.Rannyu();
	if (random <= alpha){
		route = route_1;
		accepted++;
	}
	attempted++;
}

void Simulated_Annealing :: Reset(void){//reset the accumulators and increase the steps of Simulated Annealing

	accepted = 0;
	attempted = 0;
	beta = 1./temp;
	SA_steps++;
}

void Simulated_Annealing :: Find_best_length(void){
	double length = Length(route);
	if(length < best_length){
		best_route = route;
		best_length = length;
	}
}

double Simulated_Annealing :: Length(vector<double> route){

	double sum = 0; 
	for (int i = 0; i < n_cities; i++){
		double metric = 0;
		for(int j = 0; j < 2; j++){
		
			metric += (pos[(int)route[i]][j] - pos[(int)route[(i+1)%n_cities]][j]) * (pos[(int)route[i]][j] - pos[(int)route[(i+1)%n_cities]][j]);
		}
		
		sum += sqrt(metric);
		//here we put the module in order to make the route anular
 	}
 	
 	return sum;
}

void Simulated_Annealing :: Print_best_length(int SA_step){

	ofstream outfile;
	outfile.open("data/Best_lengths.dat", ios::app);
	//we print the length of the first route in the vector population and the mean of the best lengths.
	outfile << SA_step << "\t" << temp << "\t" << best_length << endl;
	
	outfile.close();
	
}
	

void Simulated_Annealing :: Print_best_route(void){
	
	ofstream outfile;
	outfile.open("data/Best_route.dat");
	
	for( int i = 0; i < n_cities; i++){
	
		outfile << pos[(int) best_route[i]][0] << "\t" << pos[(int) best_route[i]][1] << endl;
	}
	
	outfile << pos[(int) best_route[0]][0] << "\t" << pos[(int) best_route[0]][1] << endl;
	
	outfile.close();
	
	Rand.SaveSeed();
}

