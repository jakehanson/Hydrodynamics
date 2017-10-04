#include "header.h"
#include <cmath>
#include <stdexcept>


void hydroiso_cen(std::vector<double> x,std::vector<double> xi,std::vector<double> &rho,std::vector<double> &rhou,double e,double gamma,double dt){

	int nx = x.size(); // number of cells

	// Initialize ui
	std::vector<double> ui(nx+1,0); // velocities at cell interfaces
	for(int i=1;i<nx;i++){
		ui[i]=0.5*(rhou[i]/rho[i]+rhou[i-1]/rho[i-1]);
	}

	// Compute flux for rho
	std::vector<double> fluxrho(nx+1,0); // rho flux at cell interfaces
	for(int i=1;i<nx;i++){
		if(ui[i] > 0.){
			fluxrho[i] = rho[i-1]*ui[i];
		}else{
			fluxrho[i] = rho[i]*ui[i];
		}
	}

	// Update density (uses delta xi not delta x)
	for(int i=0;i<nx;i++){
		rho[i] = rho[i] - dt/(xi[i+1]-xi[i])*(fluxrho[i+1]-fluxrho[i]);
	}

	// Get flux for rhou -- HERE IS A DISCREPENCY
	std::vector<double> fluxrhou(nx+1,0); // rhou flux at cell interfaces
	for(int i=1;i<nx;i++){
		if(ui[i] > 0){
			fluxrhou[i] = pow(rhou[i-1],2)/rho[i-1]; // rho[i] has already been updated??
		}else{
			fluxrhou[i] = pow(rhou[i],2)/rho[i];
		}
	}

	// Update momentum 
	for(int i=0;i<nx;i++){
		rhou[i] = rhou[i]-dt/(xi[i+1]-xi[i])*(fluxrhou[i+1]-fluxrhou[i]);
	}

	// Compute pressure -- DENSITY IS UPDATED)
	std::vector<double> pressure(nx,0); // rhou flux at cell interfaces
	for(int i=0;i<nx;i++){
		pressure[i] = (gamma-1.)*rho[i]*e;
	}

	// Now add pressure
	for(int i=1;i<nx-1;i++){
		rhou[i] = rhou[i]-dt*(pressure[i+1]-pressure[i-1])/(x[i+1]-x[i-1]); // delta x not delta xi
	}
	rhou[0] = rhou[0] - 0.5*dt*(pressure[1]-pressure[0])/(x[1]-x[0]);
	rhou[nx-1] = rhou[nx-1] - 0.5*dt*(pressure[nx-1]-pressure[nx-2])/(x[nx-1]-x[nx-2]);

}
