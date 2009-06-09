#include <iostream>

int main()
{
	int n = 0;
	std::cout << __sync_fetch_and_add(&n,1) << std::endl;
	return 0;
}
