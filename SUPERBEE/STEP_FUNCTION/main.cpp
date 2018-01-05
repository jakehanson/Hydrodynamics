#include "header.h"

int main(int argc, char** argv)
{
	int nx = 100; // size of domain
	double u = 1.0; // velocity
	double dt = 0.1; // timestep
	int N_steps = 300; // number of steps
	double dx = 1.0; // grid size

	std::vector<double> q(nx,0); // quantity to advect
	std::vector<double> sigma(nx,0); // slope at each cell
	double slope1,slope2; // used for slope limiter

	// Open file to store results
	std::ofstream data("data.txt");

	// Initialize Step Function
	for(int i=0;i<30;i++){
		q.at(i) = 1.0;
	}

	std::cout << "STARTING.." << std::endl;
	data << q << std::endl;

	// Evolve q through N steps
	for(int i=0;i<N_steps;i++){

		// Get slope at each cell (except boundary where slope = 0 and we use simple donor cell)
		for(int j=1;j<sigma.size()-1;j++){
			slope1 = (q.at(j+1)-q.at(j))/dx;
			slope2 = (q.at(j)-q.at(j-1))/dx;
			sigma.at(j)=maxmod(minmod(slope1,2*slope2),minmod(2*slope1,slope2)); // SUPERBEE
		}

		// Evolve all cells but boundary (assumes u > 0 since it only uses j,j-1)
		for(int j=1;j<q.size();j++){
			q.at(j) = q.at(j) - u*dt/dx*(q.at(j)-q.at(j-1))-u*dt/(2.*dx)*(sigma.at(j)-sigma.at(j-1))*(dx-u*dt);
		}

		data << q << std::endl;

	}
	
	data.close();
	std::cout << "COMPLETE" << std::endl;

}
