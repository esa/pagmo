#include"discrepancy.h"
#include<iostream>
int main() {
	
	pagmo::util::discrepancy::faure rng(2,1);
	for (size_t i=1;i<13;++i){
		std::vector<double> tmp = rng();
		std::cout << tmp[0] << " " << tmp[1] << std::endl;
	}
}
