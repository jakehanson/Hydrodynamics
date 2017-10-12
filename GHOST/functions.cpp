#include "header.h"
#include <cmath>

/* Function to advect a given quantity using advection velocity*/
void advect(std::vector<double> x,std::vector<double> xi,std::vector<double> &q,std::vector<double> ui,double dt){
	
	int nx = x.size();

	if(xi.size() != nx+1){
		throw std::runtime_error("XI ARRAY IS WRONG SIZE");
	}
	if(q.size() != nx){
		throw std::runtime_error("Q ARRAY IS WRONG SIZE");
	}
	if(ui.size() != nx+1){
		throw std::runtime_error("UI ARRAY IS WRONG SIZE");
	}
	if(dt <= 0.){
		throw std::runtime_error("TIMESTEP IS NOT POSITIVE");
	}

	// Determine r_{i-1/2} for flux limiter
	std::vector<double> r(nx+1,0);
	double dq;
	for(int i = 2;i<nx-1;i++){
		dq = q[i]-q[i-1];
		if(std::abs(dq) > 0.){
			if(ui[i] > 0.){
				r[i] = (q[i-1]-q[i-2])/dq;
			}else{
				r[i] = (q[i+1]-q[i])/dq;
			}
		}
	}

	// Determine flux limiter -- SUPERBEE (ROE, 1986)
	std::vector<double> phi(nx+1,0);
	double a;
	double b;
	for(int i=1;i<nx;i++){
		a = std::min(1.0,2.0*r[i]);
		b = std::min(2.0,r[i]);
		phi[i] = std::max({0.0,a,b});
	}

	// Construct flux	
	std::vector<double> flux(nx+1,0);
	for(int i=1;i<nx;i++){
		// Donor Cell First
		if(ui[i] > 0.0){
			flux[i] = ui[i]*q[i-1];
		}else{
			flux[i] = ui[i]*q[i];
		}
		// Then Higher Order Correction
		flux[i] = flux[i] + 0.5*std::abs(ui[i])*(1-std::abs(ui[i]*dt/(x[i]-x[i-1])))*phi[i]*(q[i]-q[i-1]);
	}

	// Update cells
	for(int i=1;i<nx-1;i++){
		q[i] = q[i] - dt*(flux[i+1]-flux[i])/(xi[i+1]-xi[i]);
	}
}

/* Function to impose BC's */
void boundary(std::vector<double> &rho, std::vector<double> &rhou,int BC_flag){
	int nx = rho.size();

	// Periodic case
	if(BC_flag == 1){
		rho[0] = rho[nx-4];
		rho[1] = rho[nx-3];
		rho[nx-2] = rho[2];
		rho[nx-1] = rho[3];
		rhou[0] = rhou[nx-4];
		rhou[1] = rhou[nx-3];
		rhou[nx-2] = rhou[2];
		rhou[nx-1] = rhou[3];
	}

	// Mirror case
	if(BC_flag == 2){
		rho[0] = rho[3];
		rho[1] = rho[2];
		rho[nx-2] = rho[nx-3];
		rho[nx-1] = rho[nx-4];
		rhou[0] = rhou[3];
		rhou[1] = rhou[2];
		rhou[nx-2] = rhou[nx-3];
		rhou[nx-1] = rhou[nx-4];
	}	
}