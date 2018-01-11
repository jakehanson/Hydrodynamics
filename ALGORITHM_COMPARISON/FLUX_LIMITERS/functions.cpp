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