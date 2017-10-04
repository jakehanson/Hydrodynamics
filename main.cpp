
#include "header.h"
#include <fstream>


/* Main Function to Evolve Ising Model and Return list of States */
int main(int argc, char** argv)
{
	// Sim Params
	int nx = 100; // size of domain
	int n_steps = 1000;
	double x0 = 0.;
	double x1 = 100.;
	double xmid = 0.5*(x1+x0);
	double t = 0;
	double dt = 0.25;
	double cfl = 0.5;
	double gamma = 7/5.;
	double e = 1.0;
	double cs = sqrt(gamma*(gamma-1)*e); // sound speed
	double dg = 0.1*(x1-x0);
	double dummy; // used to calculate dt_min
	std::vector<double> x(nx,0); // number of cells
	std::vector<double> rho(nx,0); // density
	std::vector<double> rhou(nx,0); // momentum
	std::vector<double> xi(nx+1,0); // location of cell boundaries
	std::vector<double> dx(nx,0); // cell size

	// Initialize x
	for(int i=0;i<x.size();i++){
		x[i] = x0+(x1-x0)*i/(nx-1);
	}

	// Init rho with gaussian
	for(int i=0;i<rho.size();i++){
		rho[i] = 1.0+0.3*exp(-pow(x[i]-xmid,2)/pow(dg,2));
	}

	// Init xi array (location of boundaries)
	for(int i=1;i<xi.size()-1;i++){
		xi[i] = 0.5*(x[i]+x[i-1]);
	}
	xi[0] = 2*xi[1]-xi[2];
	xi[nx] = 2*xi[nx-1]-xi[nx-2];

	// Init dx array
	for(int i=1;i<nx+1;i++){
		dx[i-1] = xi[i]-xi[i-1];
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
			dummy = dx[i]/(cs+std::abs(rhou[i]/rho[i]));
			if(i==0){
				dt = dummy;
			}
			if(dummy < dt){
				dt = dummy;
			}
		}
		t = t + dt;

		// Update rho and rhou
		hydroiso_cen(x,xi,rho,rhou,e,gamma,dt);

		for(int i=0;i<x.size();i++){
			data << t << "\t" << x[i] << "\t" << rho[i] << "\t" << rhou[i] << std::endl;
		}

	}	
	data.close();
	std::cout << "DONE." << std::endl;
}

