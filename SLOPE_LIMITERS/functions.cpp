#include "header.h"

// Function to write q array to ostream
std::ostream &operator<<(std::ostream &out,std::vector<double> q){
	for (int i=0; i<q.size()-1;i++){
		out << q.at(i) << "\t";
	}
	out << q.at(q.size()-1);
	//out << std::endl;
	return out;
}

// minmod function
double minmod(double a, double b){
	if(a*b<=0){
		return 0.;
	}else{
		if(std::abs(a)<std::abs(b)){
			return a;
		}else{
			return b;
		}
	}
}

// maxmod function
double maxmod(double a, double b){
	if(a*b<=0){
		return 0.;
	}else{
		if(std::abs(a)>std::abs(b)){
			return a;
		}else{
			return b;
		}
	}
}