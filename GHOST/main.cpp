
#include "header.h"
#include <fstream>


/* Main Function to Evolve Ising Model and Return list of States */
int main(int argc, char** argv)
{
	// Sim Params
	int nx = 104; // number of cells (including ghost cells)
	int n_ghost = 2;
	int BC_flag = 1; // 1 = Periodic BC's, 2 = Mirror BC's
	int n_steps = 1000;
	double x0 = 0.;
	double x1 = 100; // 
	double xmid = 0.5*(x1+x0); // used for gaussian density perturbation
	double t = 0;
	double dt = 0.25;
	double dx = 1.0;
	double cfl = 0.25; // courant number
	double gamma = 7/5.;
	double e = 1.0; // what does this mean if this is 1?
	double cs = sqrt((gamma-1)*e); // sound speed
	double dg = 0.1*(x1-x0);
	double dummy; // used to calculate dt_min
	std::vector<double> x(nx,0); // number of cells + 2 ghost cells for BC's
	std::vector<double> xi(nx+1,0); // location of cell boundaries
	std::vector<double> rho(nx,0); // density
	std::vector<double> rhou(nx,0); // momentum
	std::vector<double> ui(nx+1,0); // velocity at cell boundaries

	// Initialize x
	for(int i=0;i<x.size();i++){
		x[i] = i-1;
		//std::cout << "x[i] at i = " << i << "\t" << x[i] << std::endl;
	}

	// Init xi array (location of boundaries)
	for(int i=0;i<nx+1;i++){
		xi[i] = i-1.5;
		//std::cout << "xi[i] at i = " << i << "\t" << xi[i] << std::endl;
	}

	// Init rho
	for(int i=0;i<rho.size();i++){
		rho[i] = 1+0.3*exp(-pow(x[i]-xmid,2)/pow(dg,2)); // gaussian in center
	}

	// Open file to store output
	std::ofstream data("data.txt");
	data << "time\tx\trho\trhou\n";
	for(int i=0;i<x.size();i++){
		data << t << "\t" << x[i] << "\t" << rho[i] << "\t" << rhou[i] << std::endl;
	}

	// Run sim
	for(int nt=1;nt<n_steps;nt++){
		// Get dt
		for(int i=0;i<rho.size();i++){
			if(rho[i] != 0){
				dummy = cfl*dx/(cs+std::abs(rhou[i]/rho[i]));
			}else{
				dummy = cfl*dx/cs;
			}
			if(i==0){
				dt = dummy;
			}
			if(dummy < dt){
				dt = dummy;
			}
		}
		t = t + dt;

		// Impose BC's
		boundary(rho,rhou,BC_flag);

		// Compute ui
		double term_1;
		double term_2;
		for(int i=1;i<nx+1;i++){
			if(rho[i] == 0){
				term_1 = 0;
			}else{
				term_1 = rhou[i]/rho[i];
			}
			if(rho[i-1] == 0){
				term_2 = 0;
			}else{
				term_2 = rhou[i-1]/rho[i-1];
			}
			ui[i]=0.5*(term_1+term_2);
		}

		// Update rho and rhou
		try{
			advect(x,xi,rho,ui,dt);
			advect(x,xi,rhou,ui,dt);
		}
		catch(std::runtime_error &e){
			std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
			return 1;
		}

		// reimpose BC's
		boundary(rho,rhou,BC_flag);

		// Get Pressure 
		std::vector<double> pressure(nx,0);
		for(int i=0;i<nx;i++){
			pressure[i] = (gamma-1.)*rho[i]*e; // rho is updated, Pressure uses values at half step
		}

		// Update rhou
		for(int i=2;i<nx-2;i++){
			rhou[i] = rhou[i]-dt*(pressure[i+1]-pressure[i-1])/(x[i+1]-x[i-1]); // delta x not delta xi
		}

		// Reimpose BC's before writing to file
		boundary(rho,rhou,BC_flag);

		// Write to file
		for(int i=0;i<x.size();i++){
			data << t << "\t" << x[i] << "\t" << rho[i] << "\t" << rhou[i] << std::endl;
		}

	}	
	data.close();
	std::cout << "DONE." << std::endl;
}

