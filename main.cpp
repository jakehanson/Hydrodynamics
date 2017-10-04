
//#include "header.h"
#include <fstream>

#include <vector>
#include <iostream>
#include <random>
#include <stdio.h>
#include <map>


/* Main Function to Evolve Ising Model and Return list of States */
int main(int argc, char** argv)
{
	// Sim Params
	int nx = 100; // size of domain
	double dx = 1.; // grid resolution
	int n_steps = 1000; // number of time steps
	double t = 0.; // time
	double dt; // timestep
	double temp; // holds velocity for calculating min timestep
	double cfl = 0.5; // courant number
	int error_flag = 0;

	double flux_left;
	double flux_right;

	// Arrays
	std::vector<double> x(nx+1,0); // number of cells (with 2 ghost cells)
	std::vector<double> q1(nx+1,0.); // array to hold density
	std::vector<double> q1_new(nx+1,0.);  // array to hold updated q1
	std::vector<double> q2(nx+1,0.); // array to hold density*velocity
	std::vector<double> q2_new(nx+1,0.);  // array to hold updated q2
	std::vector<double> u(nx+2,0.);  // array to hold updated q2

	// Pressure
	double gamma = 7./5.;
	double e = 1.;
	double cs = sqrt(gamma*(gamma-1)*e); // should be (gamma-1)*e

	// Init x_array
	for(int i=0;i<x.size();i++){
		x[i] = 100*double(i)/100.;
	}

	// Init q1_array with gaussian centered at midpoint
	for(int i=0;i<q1.size();i++){
		q1[i] = 1.0+0.3*exp(-pow((x[i]-double(nx)/2.)/10,2));
	}

	// Init velocity at interfaces
	for(int i=1;i<u.size()-1;i++){
		u[i] = 0.5*(q2[i-1]/q1[i-1]+q2[i]/q1[i]);
	}
	//BC's
	u[0] = 0;
	u[u.size()-1] = 0;

	// Open file to store output
	std::ofstream data("data.txt");
	data << "time\tx\tq1\tq2\n";
	for(int i=0;i<q1.size();i++){
		data << t << "\t" << x[i] << "\t" << q1[i] << "\t" << q2[i] << std::endl;
	}	

	// Run Sim
	for(int i=1;i<n_steps;i++){
		if(error_flag != 0){
			break;
		}
		// Get dt from max velocity
		for(int j=0;j<q1.size();j++){
			if(q2[i] == 0 and q1[i] == 0){
				temp = cfl*dx/cs;
			}else{
				temp = cfl*dx/(cs+std::abs(q2[i]/q1[i]));
			}
			if(isnan(temp)){
				std::cout << "ERROR! NAN VALUE FOR TIME\t" << q2[i] << "\t" << q1[i] << std::endl;
				error_flag = 1;
				break;
			}
			if(j==0){
				dt = temp;
			}
			if(temp < dt){
				dt = temp;
			}
			if(dt < 0){
				std::cout << "ERROR! NEGATIVE TIME STEP" << std::endl;
				error_flag = 1;
				break;
			}
		}
		// Update time
		t = t + dt;


		/* HALF STEP - OPERATOR SPLITTING PART 1 */
		// Update q1
		for(int i=1;i<q1.size()-1;i++){
			if(u[i+1]>0){
				flux_right = q1[i]*u[i+1];
			}else{
				flux_right = q1[i+1]*u[i+1];
			}
			if(u[i]>0){
				flux_left = q1[i-1]*u[i];
			}else{
				flux_left = q1[i]*u[i];
			}
			q1_new[i] = q1[i] - dt/dx*(flux_right-flux_left);
		}
		//BC's
		if(u[1] > 0){
			q1_new[0] = q1[0]-dt/dx*(q1[0]*u[1]);
		}else{
			q1_new[0] = q1[0]-dt/dx*(q1[1]*u[1]);
		}	
		if(u[u.size()-2] > 0){
			q1_new[q1.size()-1] = q1[q1.size()-1]+dt/dx*(q1[q1.size()-2]*u[u.size()-2]);
		}else{
			q1_new[q1.size()-1] = q1[q1.size()-1]+dt/dx*(q1[q1.size()-1]*u[u.size()-2]);
		}	
		// Update q2
		for(int i=1;i<q1.size()-1;i++){
			if(u[i+1]>0){
				flux_right = q2[i]*u[i+1];
			}else{
				flux_right = q2[i+1]*u[i+1];
			}
			if(u[i]>0){
				flux_left = q2[i-1]*u[i];
			}else{
				flux_left = q2[i]*u[i];
			}
			q2_new[i] = q2[i] - dt/dx*(flux_right-flux_left);
		}
		//BC's
		if(u[1] > 0){
			q2_new[0] = q2[0]-dt/dx*(q2[0]*u[1]);
		}else{
			q2_new[0] = q2[0]-dt/dx*(q2[1]*u[1]);
		}	
		if(u[u.size()-2] > 0){
			q2_new[q1.size()-1] = q2[q1.size()-1]+dt/dx*(q2[q1.size()-2]*u[u.size()-2]);
		}else{
			q2_new[q1.size()-1] = q2[q1.size()-1]+dt/dx*(q2[q1.size()-1]*u[u.size()-2]);
		}	


		/* FULL STEP - OPERATOR SPLITTING PART 2 */
		// Update q1
		for(int i=0;i<q1.size();i++){
			q1[i] = q1_new[i];
			if(q1[i]<0){
				std::cout << "ERROR! NEGATIVE DENSITY" << std::endl;
				error_flag = 1;
				break;
			}
		}
		// Update q2
		for(int i=1;i<q2.size()-1;i++){
			q2[i] = q2_new[i] - pow(cs,2)/(2*dx)*(q1_new[i+1]-q1_new[i-1]);
		}
		// BC's
		q2[0] = q2_new[0] - pow(cs,2)/dx*(q1_new[1]-q1_new[0]);
		q2[q2.size()-1] = q2_new[q2.size()-1] - pow(cs,2)/dx*(q1_new[q2.size()-1]-q1_new[q2.size()-2]);
		// Update velocity
		for(int i=1;i<u.size()-1;i++){
			u[i] = 0.5*(q2[i-1]/q1[i-1]+q2[i]/q1[i]);
		}


		/* Write output */
		for(int i=0;i<q1.size();i++){
			data << t << "\t" << x[i] << "\t" << q1[i] << "\t" << q2[i] << std::endl;
		}	
	}

	data.close();
	std::cout << "Done.\n";

}