#pragma once
#include <vector>
#include <iostream>
#include <random>
#include <stdio.h>
#include <map>
#include <stdexcept>

void advect(std::vector<double> x,std::vector<double> xi,std::vector<double> &q,std::vector<double> ui,double dt);

void boundary(std::vector<double> &rho, std::vector<double> &rhou,int BC_flag);
