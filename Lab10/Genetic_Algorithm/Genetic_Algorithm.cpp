#include "Genetic_Algorithm.h"

using namespace std;

Genetic_Algorithm :: Genetic_Algorithm(): Rand(SEEDDIR "/Primes", SEEDDIR "/seed.in"){ }

Genetic_Algorithm :: ~Genetic_Algorithm(){ }
 

void Genetic_Algorithm :: Input(void){//initialization

	ifstream ReadInput;
	ReadInput.open("Genetic_Algorithm/input.dat");
	if (!ReadInput.is_open()) {
		cerr << "Error: unable to open input.dat" << endl;
		exit(1);
	}
	
	//the cout is only for Rank == 0, while information about parameters is common
	if (Rank == 0) {
		cout << " Travelling Salesman Problem: " << "solution with the Genetic Algorithm and parallel computing"<< endl;
	
		cout << endl << "Number of nodes used = " << Size << endl;
	}
	
	ReadInput >> n_cities;
	ReadInput >> distribution;
	
	if (Rank == 0) {
		cout << "Number of cities : " << n_cities << endl;
		cout << "Configuration: " << distribution << endl;
	}
	
	ReadInput >> n_routes; 
	ReadInput >> n_generations;
	
	if (Rank == 0) {
		cout << "Number of routes = " << n_routes << endl;
		cout << "Number of generations = " << n_generations << endl;
	}
	
	ReadInput >> p_crossover;
	ReadInput >> p_mutation;
	ReadInput >> power;
	
	if (Rank == 0) {
		cout << "Probability of crossing over = " << p_crossover*100 << "%" << endl;
		cout << "Probability of mutation = " << p_mutation*100 << "%" << endl;
		cout << "Crossover power = " << power << endl;
	}
	
	ReadInput >> n_migrations;
	
	if (Rank == 0) {
	
		cout << "Every Migration takes place after " << n_migrations << " generations " << endl;
	}
	

	//set the configuration of the cities and print them
	Cities();
	
	ReadInput.close();
	
	Rand.SetRandom(Rank);
	
	
}

void Genetic_Algorithm :: Cities(void){

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
	ofstream outfile("data/data_10.2/pos_cities.dat");
	
	for (int i = 0; i < n_cities; i++){
		outfile << pos[i][0] << "\t" << pos[i][1] << endl;
	}
	//the last point has to be the intial point 
	outfile << pos[0][0] << "\t" << pos[0][1] << endl;
	
	outfile.close();
	
}

void Genetic_Algorithm :: Create_Population(void){

	for(int i=0; i < n_routes; i++){
		
		//create a population
		vector<double> starting_route = Create_Route(n_cities);
		population.push_back(starting_route);
		
	}
	
	Assign_Length();
	My_Sort();
	
}

void Genetic_Algorithm :: My_Sort(void){

	//order the population
	sort(population.begin(), population.end(), [](const vector<double>& a, const vector<double>& b) { return (a[a.size()-1]) < (b[b.size()-1]); });
	
}


vector<double> Genetic_Algorithm :: Create_Route(int n_cities){

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
	//Assing_Length() method is supposed to do so, here we only create the space for that value.

	return route;
}

void Genetic_Algorithm :: New_Generation(double p_crossover, double p_mutation){

	vector<vector<double>> new_population;
	
	for(int i = 0; i < n_routes/2; i++){// n_routes/2 and create 2 sons from two parents: every generation has to have the same number of individuals
	
		int mum, dad;
		vector<vector<double>> children;
		
		mum = Selection();
		dad = Selection();
		
		children = {population[mum], population[dad]};
		
		if(Rand.Rannyu() < p_crossover){
			children = Crossover(mum,dad);
			//cout << " CROSSOVER AVVENUTO " << endl;
		} else {
		
			//cout << " CROSSOVER NON AVVENUTO " << endl;
		}
		//With the crossover we create the two children
		//now we perform the mutation on the children
		
		
		//first son
		
		//first mutation
		if(Rand.Rannyu() < p_mutation){ //permutation of a couple of cities 
			int n1 = (int)Rand.Rannyu(1,n_cities);
			int n2 = (int)Rand.Rannyu(1,n_cities);
			
			swap(children[0][n1], children[0][n2]);
		
			//cout << " mutation1 AVVENUTO " << endl;
		} else {
			//cout << " mutation1 NON AVVENUTO " << endl;
		}
		
		//second mutation
		if(Rand.Rannyu() < p_mutation){ //Inversion of the vector
		
			int n1 = (int)Rand.Rannyu(1,n_cities);
			int n2 = (int)Rand.Rannyu(1,n_cities);
			
			if(n2 > n1){
				reverse(children[0].begin() + n1, children[0].begin() + n2);
				//cout << " mutation2.1 AVVENUTO " << endl;
			} else {
				reverse(children[0].begin() + n2, children[0].begin() + n1);
				//cout << " mutation2.2 AVVENUTO " << endl;
			}
		} else{ 
			//cout << " mutation2 NON AVVENUTO " << endl;
		}
		
		//third mutation	
		if(Rand.Rannyu() < p_mutation){//permutation of adjacent
				
			vector<int> s;
			//we want to use rotate, so we have to random generate the first, the last and the middle iterators. Rotate will act rotating the vector until "middle" becomes the first element of the vector.
			
			s.push_back((int)Rand.Rannyu(1, n_cities));
			s.push_back((int)Rand.Rannyu(1, n_cities));
			s.push_back((int)Rand.Rannyu(1, n_cities));
			sort(s.begin(), s.end());// we sort the vector s in order to have first, middle and last ordered
			
			rotate(children[0].begin()+s[0], children[0].begin()+s[1], children[0].begin()+s[2]);
			
		} else{
		
			//cout << " mutation 3 NON AVVENUTO " << endl;
		}
		
		//second son
		
		//first mutation
		if(Rand.Rannyu() < p_mutation){ //permutation of a couple of cities 
			int n1 = (int)Rand.Rannyu(1,n_cities);
			int n2 = (int)Rand.Rannyu(1,n_cities);
			swap(children[1][n1], children[1][n2]);
			
			//cout << " mutation1 AVVENUTO " << endl;
		} else {
			//cout << " mutation1 NON AVVENUTO " << endl;
		}
		
		//second mutation
		if(Rand.Rannyu() < p_mutation){ //Inversion of the vector
		
			int n1 = (int)Rand.Rannyu(1,n_cities);
			int n2 = (int)Rand.Rannyu(1,n_cities);
			
			if(n2 > n1){
				reverse(children[1].begin() + n1, children[1].begin() + n2);
				//cout << " mutation2.1 AVVENUTO " << endl;
			} else {
				reverse(children[1].begin() + n2, children[1].begin() + n1);
				//cout << " mutation2.2 AVVENUTO " << endl;
			}
		} else { 
			//cout << " mutation2 NON AVVENUTO " << endl;
		}
			
		if(Rand.Rannyu() < p_mutation){//permutation of adjacent
				
			vector<int> s;
			s.push_back((int)Rand.Rannyu(1, n_cities));
			s.push_back((int)Rand.Rannyu(1, n_cities));
			s.push_back((int)Rand.Rannyu(1, n_cities));
			sort(s.begin(), s.end());
				
			rotate(children[1].begin()+s[0], children[1].begin()+s[1], children[1].begin()+s[2]);
			//cout << " mutation 3 AVVENUTO " << endl;
		} else{
		
			//cout << " mutation 3 NON AVVENUTO " << endl;
		}
		
		//cout << endl << "FINE DELLE MUTATION " << endl;
		
		new_population.push_back(children[0]);
		new_population.push_back(children[1]);
			
	}
	
	for (int i = 0; i < n_routes; i++){
		fill(population[i].begin(), population[i].end(), 0.0);
	}
	
	population = new_population;
	
	Assign_Length();
	
	My_Sort();
			
}

vector<vector<double>> Genetic_Algorithm :: Crossover(int mum, int dad){

	vector<double> route_mum = population[mum];
	vector<double> route_dad = population[dad];
	vector<double> son_1;
	vector<double> son_2;
	
	son_1.reserve(route_mum.size());
	copy(route_mum.begin(), route_mum.end() -1, back_inserter(son_1));
	son_1.push_back(0.0); //element relative to the length of the route
	
	son_2.reserve(route_dad.size());
	copy(route_dad.begin(), route_dad.end()-1, back_inserter(son_2));
	son_2.push_back(0.0); //element relative to the length of the route
	
	//random create a cutting point
	int alt = (int)Rand.Rannyu(n_cities/2 - 3, n_cities/2 + 3);
	
	int index = alt; 
	
	for(int i = 0; i < n_cities; i++){
		
		int single_city = route_dad[i];
		if (find(route_mum.begin(), route_mum.begin() + alt, single_city) == route_mum.begin() + alt){
			// we change the index element because find runs in [begin, last)
			son_1[index] = single_city;
			index++;
		}
	}
	
	index = alt;
	for(int i = 0; i < n_cities; i++){
	
		int single_city = route_mum[i];
		if (find(route_dad.begin(), route_dad.begin() + alt, single_city) == route_dad.begin() + alt){
			son_2[index] = single_city;
			index++;
		}
	}
	//cout << "FINE CROSSOVER" << endl;
	return {son_1, son_2};

}
//Once we sorted the population we choose mum and dad with selection(), resorting to the advices
//the greater is the power and the more demanding is the selection
int Genetic_Algorithm :: Selection(void){

	double r = Rand.Rannyu();
	return (int) (n_routes * pow(r, power));
	//we use such value for selection because if we put "+1" then we go beyond the vector limits
	
}

double Genetic_Algorithm :: Length(vector<double> route){

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

void Genetic_Algorithm :: Print_best_length(int step){

	vector<double> lengths;
	int n_best_routes = (int)(n_routes/2);
	//we take only the first half of the vector since it is sorted
	
	for(int i = 0; i < n_best_routes; i++){
		lengths.push_back(population[i][n_cities]);
	}
	
	ofstream outfile;
	outfile.open("data/data_10.2/rank" + to_string(Rank) +"/ Best_lengths.dat", ios::app);
	//we print the length of the first route in the vector population and the mean of the best lengths.
	outfile << step << "\t" << *(lengths.begin()) << "\t" << accumulate(lengths.begin(), lengths.end(), 0.0)/n_best_routes << endl;
	
	outfile.close();
	
}

void Genetic_Algorithm :: Print_best_route(void){

	vector<double> best_route;
	vector<vector<double>> best_pos;
	
	ofstream outfile;
	outfile.open("data/data_10.2/rank" +to_string(Rank) +"/Best_route.dat");
	
	best_route = population[0];
	for( int i = 0; i < n_cities; i++){
	
		outfile << pos[(int) best_route[i]][0] << "\t" << pos[(int) best_route[i]][1] << endl;
	}
	
	outfile << pos[(int) best_route[0]][0] << "\t" << pos[(int) best_route[0]][1] << endl;
	
	outfile.close();
	
	Rand.SaveSeed();
}

void Genetic_Algorithm :: Assign_Length(void){
	
	for(int i = 0; i< n_routes; i++){
		population[i][n_cities] = Length(population[i]);
	}

}

void Genetic_Algorithm :: Migration(void){//different nodes communicate their best routeto to each other

	//create the sending and the receiving vector and fill teh sending one with the population
	vector<double> migrating_routes;
	migrating_routes.resize(n_cities+1);
	
	
	for ( int i = 0; i < n_cities+1; i++){
		migrating_routes[i] = population[0][i];
	}
	

	vector<double> receiving_routes;
	receiving_routes.resize(n_cities+1);
	
	//initialize the number of continent. They will indicate the sending and the receiving nodes
	int continent1 = -1;
	int continent2 = -1;
	int continent3 = -1;
	
	//take a random continent
	if (Rank == 0){
		continent2 = 0;
		do{
			continent1 = int(Rand.Rannyu()*Size);
		} while (continent1 == continent2);
		
		continent3 = continent1;
	}
	//give information about the random continent to all the ranks
	MPI_Bcast(&continent1, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (Rank == continent1){
		continent2 = 0;
		continent3 = continent1;
	}
	//set all the possibilities
	if (Rank != continent1 and Rank != 0){
		if (continent1 == 1){
			continent2 = 2;
			continent3 = 3;
		} else if (continent1 == 2){
			continent2 = 1;
			continent3 = 3;
		} else if (continent1 == 3){
			continent2 = 1;
			continent3 = 2;
		}
	}
	
	//define status and tags for sending and receiving methods of MPI
	MPI_Status status1, status2;
	int itag1 = 1;
	int itag2 = 2;	
		
	MPI_Barrier(MPI_COMM_WORLD);//wait for all the nodes 
	
	if (Rank == continent2){
	
			MPI_Send(migrating_routes.data(), n_cities+1, MPI_DOUBLE, continent3, itag1, MPI_COMM_WORLD); 
			MPI_Recv(receiving_routes.data(), n_cities+1, MPI_DOUBLE, continent3, itag2, MPI_COMM_WORLD, &status2);
			//cout << endl << "messaggio da c2 a c3 mandato e da c3 a c2 ricevuto" << endl;
		
	} 
	else if (Rank == continent3){
			MPI_Send(migrating_routes.data(), n_cities+1, MPI_DOUBLE, continent2, itag2, MPI_COMM_WORLD);
			MPI_Recv(receiving_routes.data(), n_cities+1, MPI_DOUBLE, continent2, itag1, MPI_COMM_WORLD, &status1);
	}
		
	
//============= unused alternative with Isend and Irecv	================
	
	/*int n = 100;
	vector<int> imesg1;// = new int[n]; 
	vector<int> imesg2;// = new int[n];
	imesg1.resize(n);
	imesg2.resize(n);
	for(int i = 0; i < n; i++){
		imesg1[i] = Rank;
		imesg2[i] = Rank +1;
	}
	
	if (Rank == 1){
		MPI_Send(imesg1.data(), n, MPI_INTEGER, 0, itag1, MPI_COMM_WORLD);
		MPI_Recv(imesg2.data(), n, MPI_INTEGER, 0, itag2, MPI_COMM_WORLD, &status2);
		cout << endl << "messaggio = " << imesg2[0] << endl;
	}
	else if (Rank == 0) {
		MPI_Send(imesg2.data(), n, MPI_INTEGER, 1, itag2, MPI_COMM_WORLD);
		MPI_Recv(imesg1.data(), n, MPI_INTEGER, 1, itag1, MPI_COMM_WORLD, &status1);
		cout << endl << " messaggio = " << imesg1[0] << endl;
	}*/
	
	MPI_Barrier(MPI_COMM_WORLD);	
		
	/*MPI_Request request;
	MPI_Status status;
	int request_complete = 0;
	
	MPI_Request request1;
	MPI_Status status1;
	
	MPI_Request request2;
	MPI_Status status2;
	int request_complete2 = 0;
	
	MPI_Request request3;
	MPI_Status status3;
	
	MPI_Win win;
	MPI_Errhandler errhandler;
	MPI_Win_set_errhandler(win, errhandler);
	
	if (Rank == continent2){
		MPI_Isend(&migrating_routes, n_cities+1, MPI_DOUBLE, continent3, 1, MPI_COMM_WORLD, &request);
		
		if (!request_complete)
			MPI_Test(&request, &request_complete, &status);
	
		if (!request_complete)
			MPI_Wait(&request, &status);
			
		MPI_Irecv(&receiving_routes, n_cities+1, MPI_DOUBLE, continent3, 2, MPI_COMM_WORLD, &request1);
		
		MPI_Wait(&request1, &status1);
	} else if (Rank == continent3){
		MPI_Isend(&migrating_routes, n_cities+11, MPI_DOUBLE, continent2, 2, MPI_COMM_WORLD, &request);
		
		if (!request_complete2)
			MPI_Test(&request2, &request_complete2, &status2);
	
		if (!request_complete2)
			MPI_Wait(&request2, &status2); 
		MPI_Irecv(&receiving_routes, n_cities+1, MPI_DOUBLE, continent2, 1, MPI_COMM_WORLD, &request3);
		
		MPI_Wait(&request3, &status3);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);*/
//==================================================================
	//set the population to the received one
	for (int j = 0; j < n_cities+1; j++){
		population[0][j] = receiving_routes[j];
	}
	
	My_Sort();

	
}

/*void Genetic_Algorithm :: Broadcasting(void){//unused alternative for sending and receive messages

	MPI_Bcast(&n_cities, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	//since we want to make things easy: tranform the "distribution" string t a int:
	int config = 0;
	if( distribution == "square"){
		config = 1;
	}
	
	MPI_Bcast(&config, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_routes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_generations, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&p_crossover, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&p_mutation, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&power, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_migrations, 1, MPI_INT, 0, MPI_COMM_WORLD);
}*/

	


