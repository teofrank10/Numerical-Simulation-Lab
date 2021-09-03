#include <iostream>
#include <fstream>
#include <random.h>
#include <vector>

using namespace std; 

//==========================Wave functions=========================//
double psi_100(double pos[], int dim){
	double r = 0;
	for (int i = 0; i < dim; i++){
		r += pos[i]*pos[i];
	}
	return exp(-2.* sqrt(r));
}

double psi_210(double pos[], int dim){
	double r = 0;
	for (int i = 0; i < dim; i++){
		r += pos[i]*pos[i];
	}
	return pow(pos[dim-1],2)*exp(-sqrt(r));
}
//=================================================================//

//========================Blocking functions=======================//
//same functions used in previous notebooks and written in python.
//here we implemented them in c++ code
double error(double ave, double ave2, double n){

	if(n == 0){return 0;}else {
		return sqrt((ave2-ave*ave)/n);
	}
}

double** observable(const vector<double>& v){

	double sum = 0.; 
	double sum2 = 0.;
	double ave = 0.;
	double ave2 = 0.;
	int size = static_cast<int>(v.size());
	double** val = 0;
	val = new double*[2];
	for (int i = 0; i < 2; i++){
		val[i] = new double[size];
	}
		
	for (int i = 0; i < size; i++){
	
		sum += v[i];
		sum2 += v[i]*v[i];
		ave = sum / (i+1);
		ave2 = sum2 / (i+1);
		val[0][i] = ave;
		val[1][i] = error(ave, ave2, i);
	}
	
	return val;
}


double partial_ave(const vector<double>& v, int start, int interval){//make average on a certain slice of vector


	double sum = 0.;
	for (int i = start; i < start + interval; i++){
	
		sum += v[i];
	}
	
	return sum / interval;
	
}

//=================================================================//


int main (){

	const int N_blocks  = 100;
	const int M_throws = 1000000;
	const int L = M_throws/ N_blocks; //throws per block
	const int M_equilib = 1510;//equilibration step
	
	//create vectors
	vector<double> psi_100_r(M_throws);
	vector<double> psi_210_r(M_throws);
	vector<double> psi_100_obs(N_blocks);
	vector<double> psi_210_obs(N_blocks);
	vector<double> psi_100_ave(N_blocks);
	vector<double> psi_210_ave(N_blocks);
	vector<double> psi_100_error(N_blocks);
	vector<double> psi_210_error(N_blocks);
	
	const double sigma_100 = 0.759;//set parameters for the Gaussian sampling
	const double sigma_210 = 1.874;
	
	double alpha_100_ave = 0;//set to zero the mean of the acceptance
	double alpha_210_ave = 0;
	
	//Create a random number generator
	Random Rand (SEEDDIR "/Primes", SEEDDIR "/seed.in");
	
	//Save information about the numbers of our generation on a file
	ofstream outfile("data/Gauss/info_ex5.1.dat");
	if (!outfile.is_open()) {
		cerr << "Error: unable to open info_ex5.1.dat" << endl;
		exit(1);
	}
	outfile << "Number of blocks: " << N_blocks << endl << "Number of throws: " << M_throws << endl << "Number of throws per block: " << L << endl;
	outfile.close();
	
	//Equilibration for psi_100
	
	ofstream outfile_points;
	ofstream outfile_radii;
	
	outfile_points.open("data/Gauss/eq_positions100.dat");
	if (!outfile_points.is_open()){
		cout << "Error: unable to open eq_positions100.dat" << endl;
		exit(1);
	}
	
	outfile_radii.open("data/Gauss/eq_radii100.dat");
	if (!outfile_radii.is_open()){
		cout << "Error: unable to open eq_radii100.dat" << endl;
		exit(1);
	}
	
//===================Uniform sampling=====================//
	//starting point for psi_100: {0.2, 0.2, 0.2}
	//code for the equilibration starting from a far point:
	/*double pos_100[3];
	pos_100[0] = 40;
	pos_100[1] = 30;
	pos_100[2] = 2;*/
	 
	double pos_100[3];
	for (int i = 0; i < 3; i++){
		pos_100[i] = 0.2;
	}
	
	for (int i = 0; i < M_equilib; i++){//use metropolis to equilibrate the system
		Rand.Metropolis_Gauss(pos_100, sigma_100, psi_100, 3);
		
		double x_100 = Rand.Metropolis_GetX();
		double y_100 = Rand.Metropolis_GetY();
		double z_100 = Rand.Metropolis_GetZ();
		//update the vector pos:
		pos_100[0] = Rand.Metropolis_GetX();
		pos_100[1] = Rand.Metropolis_GetY();
		pos_100[2] = Rand.Metropolis_GetZ();
		double r_100 = sqrt(pow(x_100,2)+pow(y_100,2)+pow(z_100,2));//calculate radius
		
		alpha_100_ave = i/double(i+1)*alpha_100_ave + 1./double(i+1)*Rand.Metropolis_GetAlpha();//calculate the average of acc. rate
		
		outfile_points << x_100 << "\t" << y_100 << "\t" << z_100 << endl;
		outfile_radii << r_100 << endl;
	}
	cout << "alpha eq 100: " << Rand.Metropolis_GetAlpha() << endl;//instantaeous acc. rate
	cout << "alpha eq 100AVE: " << alpha_100_ave << endl;//acc. rate to be set to 0.5
	outfile_points.close();
	outfile_radii.close();
	
	//Equilibration for psi_210
	
	outfile_points.open("data/Gauss/eq_positions210.dat");
	if (!outfile_points.is_open()){
		cout << "Error: unable to open eq_positions210.dat" << endl;
		exit(1);
	}
	
	outfile_radii.open("data/Gauss/eq_radii210.dat");
	if (!outfile_radii.is_open()){
		cout << "Error: unable to open eq_radii210.dat" << endl;
		exit(1);
	}
	
	
	//starting point for psi_210: {1.2, 1.2, 1.2}
	
	//code for the equilibration starting from a far point:
	/*double pos_210[3];
	pos_210[0] = 40;
	pos_210[1] = 30;
	pos_210[2] = 2;*/
	
	
	double pos_210[3];
	for (int i = 0; i < 3; i++){
		pos_210[i] = 1.2;
	}
	
	for (int i = 0; i < M_equilib; i++){
		Rand.Metropolis_Gauss(pos_210, sigma_210, psi_210, 3);
		
		double x_210 = Rand.Metropolis_GetX();
		double y_210 = Rand.Metropolis_GetY();
		double z_210 = Rand.Metropolis_GetZ();
		//update the vector pos:
		pos_210[0] = Rand.Metropolis_GetX();
		pos_210[1] = Rand.Metropolis_GetY();
		pos_210[2] = Rand.Metropolis_GetZ();
		double r_210 = sqrt(pow(x_210,2)+pow(y_210,2)+pow(z_210,2));
		
		alpha_210_ave = i/double(i+1)*alpha_210_ave + 1./double(i+1)*Rand.Metropolis_GetAlpha();
		
		outfile_points << x_210 << "\t" << y_210 << "\t" << z_210 << endl;
		outfile_radii << r_210 << endl;
	}
	cout << "alpha eq 210: " << Rand.Metropolis_GetAlpha() << endl;
	cout << "alpha eq 210AVE: " << alpha_210_ave << endl;
	outfile_points.close();
	outfile_radii.close();
	
	
	//sampling 
	ofstream outfile_100;
	ofstream outfile_210;
	
	outfile_100.open("data/Gauss/positions_100.dat");
	if (!outfile_100.is_open()){
		cout << "Error: unable to open positions_100.dat" << endl;
		exit(1);
	}
	
	outfile_210.open("data/Gauss/positions_210.dat");
	if (!outfile_210.is_open()){
		cout << "Error: unable to open positions_210.dat" << endl;
		exit(1);
	}
	
	for (int i = 0; i < M_throws; i++){//simulation
		Rand.Metropolis_Gauss(pos_100, sigma_100, psi_100, 3);
		
		double x_100 = Rand.Metropolis_GetX();
		double y_100 = Rand.Metropolis_GetY();
		double z_100 = Rand.Metropolis_GetZ();
		//update the vector pos:
		pos_100[0] = Rand.Metropolis_GetX();
		pos_100[1] = Rand.Metropolis_GetY();
		pos_100[2] = Rand.Metropolis_GetZ();
		double r_100 = sqrt(pow(x_100,2)+pow(y_100,2)+pow(z_100,2));
		
		psi_100_r[i] = r_100;//save the radius in the vector
		
		alpha_100_ave = i/double(i+1)*alpha_100_ave + 1./double(i+1)*Rand.Metropolis_GetAlpha();
		
		if(i%200 == 0 ){//every 200 throws print the position
		outfile_100 << x_100 << "\t" << y_100 << "\t" << z_100 << endl;
		
		}
	}
	cout << "alpha sampling 100: " << Rand.Metropolis_GetAlpha() << endl;
	cout << "alpha sampling 100AVE: " << alpha_100_ave << endl;
	outfile_100.close();
	
	for (int i = 0; i < M_throws; i++){
		Rand.Metropolis_Gauss(pos_210, sigma_210, psi_210, 3);
		
		double x_210 = Rand.Metropolis_GetX();
		double y_210 = Rand.Metropolis_GetY();
		double z_210 = Rand.Metropolis_GetZ();
		//update the vector pos:
		pos_210[0] = Rand.Metropolis_GetX();
		pos_210[1] = Rand.Metropolis_GetY();
		pos_210[2] = Rand.Metropolis_GetZ();
		double r_210 = sqrt(pow(x_210,2)+pow(y_210,2)+pow(z_210,2));
		
		psi_210_r[i] = r_210;
		
		alpha_210_ave = i/double(i+1)*alpha_210_ave + 1./double(i+1)*Rand.Metropolis_GetAlpha();
		
		if(i%200 == 0 ){//print 5000 positions
		outfile_210 << x_210 << "\t" << y_210 << "\t" << z_210 << endl;
		
		}
	}
	cout << "alpha sampling 210: " << Rand.Metropolis_GetAlpha() << endl;
	cout << "alpha sampling 210AVE: " << alpha_210_ave << endl;
	outfile_210.close();
	
	//calculate the N_blocks observables:
	
	for (int i = 0; i < N_blocks; i++){//make average with partial_Ave
	
		psi_100_obs[i] = partial_ave(psi_100_r, i*L, L);
		psi_210_obs[i] = partial_ave(psi_210_r, i*L, L);
		
	}
	
	
	//calculate the average and the error with the "observable" function
	
	double** vett_100 = observable(psi_100_obs);
	for ( int i = 0; i < N_blocks; i++){
		psi_100_ave[i] = vett_100[0][i];
		psi_100_error[i] = vett_100[1][i];
	}
	
	double** vett_210 = observable(psi_210_obs);
	for ( int i = 0; i < N_blocks; i++){
		psi_210_ave[i] = vett_210[0][i];
		psi_210_error[i] = vett_210[1][i];
	} 
	
	outfile_100.open("data/Gauss/r_100.dat");
	if (!outfile_100.is_open()){
		cout << "Error: unable to open positions_100.dat" << endl;
		exit(1);
	}
	
	outfile_210.open("data/Gauss/r_210.dat");
	if (!outfile_210.is_open()){
		cout << "Error: unable to open positions_210.dat" << endl;
		exit(1);
	}
	
	for (int i = 0; i < N_blocks; i++){
		
		outfile_100 << psi_100_ave[i] << "\t" << psi_100_error[i] << endl;
		outfile_210 << psi_210_ave[i] << "\t" << psi_210_error[i] << endl;
	}
//===================end of the uniform sampling============//	

		
	Rand.SaveSeed();
	return 0;
}


