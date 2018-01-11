#include "header.h"

int main(int argc, char** argv)
{
	int nx = 100; // size of domain
	int nb = 2; // number of ghost cells on either side of domain
	nx = nx+2*nb;
	double u = 1.0; // velocity
	double dt = 0.1; // timestep
	int N_steps = 300; // number of steps
	//int N_steps = 1; // number of steps
	double dx = 1.0; // grid size

	std::vector<double> q(nx,0); // quantity to advect
	std::vector<double> r(nx,0); // r value at each cell
	double flux_left,flux_right; // used in advection
	double phi_left, phi_right; // flux limiters
	double r_left, r_right; // slope differences

	// Open file to store results
	std::ofstream data("data.txt");

	// Initialize Step Function
	q.at(0) = 0.;
	q.at(1) = 0.;
	for(int i=nb;i<=30+nb;i++){
		q.at(i) = 1.0;
	}

	std::cout << "STARTING.." << std::endl;
	data << q << std::endl;

	// Evolve q through N steps
	for(int i=0;i<N_steps;i++){

		// Evolve all cells but boundary (assumes u > 0 since it only uses j,j-1)
		for(int j=nb;j<q.size()-nb;j++){

			r_left = (q.at(j-1)-q.at(j-2))/(q.at(j)-q.at(j-1));
			r_right = (q.at(j)-q.at(j-1))/(q.at(j+1)-q.at(j));

			/* Superbee */
			// phi_left = std::max(0.,std::max(std::min(1.,2.*r_left),std::min(2.,r_left)));
			// phi_right = std::max(0.,std::max(std::min(1.,2.*r_right),std::min(2.,r_right)));

			/* Donor Cell */
			// phi_left = 0.;
			// phi_right = 0.;

			/* Lax-Wendroff (Downwind) */
			// phi_left = 1.;
			// phi_right = 1.;

			/* Monotonized Central-Difference (MC) */
			// phi_left = std::max(0.,std::min((1.+r_left)/2.,std::min(2.,2.*r_left)));
			// phi_right = std::max(0.,std::min((1.+r_right)/2.,std::min(2.,2.*r_right)));

			/* van Leer */
			phi_left = (r_left+std::abs(r_left))/(1.+std::abs(r_left));
			phi_right = (r_right+std::abs(r_right))/(1.+std::abs(r_right));

			flux_left = u*q.at(j-1)+u/2.*(1-u*dt/dx)*phi_left*(q.at(j)-q.at(j-1));	
			flux_right = u*q.at(j)+u/2.*(1-u*dt/dx)*phi_right*(q.at(j+1)-q.at(j));	

			q.at(j) = q.at(j) + dt/dx*(flux_left-flux_right);
		}

		data << q << std::endl;

	}
	
	data.close();
	std::cout << "COMPLETE" << std::endl;

}
