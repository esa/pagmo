#include"discrepancy.h"
#include<iostream>
int main() {	
	pagmo::util::discrepancy::simplex rng(3,1);
	for (size_t i=1;i<40;++i){
		std::vector<double> tmp = rng();
		std::cout << tmp[0] << " " << tmp[1] << " " << tmp[2] << std::endl;
	}
}
