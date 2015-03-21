#include<iostream>
#include<random>
#include"Distributions.h"
#include"Sampler.hpp"
#include<vector>
#include<algorithm>
#include<fstream>
#include<string>
#include<numeric>
#include<array>
#include<limits>
//#include<gsl/gsl_sf_psi.h>
template<int m,int n>
using Matrix = std::array < std::array<double, n>, m > ;
int test() 
{
	std::vector<int> test(100);
	std::iota(test.begin(), test.end(), 1);
	std::for_each(test.begin(), test.end(), [](int n){std::cout << n << std::endl; });
	auto test2 = test;
	test2.insert(test2.end(), test.begin(), test.end());
	std::for_each(test2.begin(), test2.end(), [](int n){std::cout << n << std::endl; });
	int ans = std::accumulate(test2.begin(), test2.end(), 0, [](int s, int now)->int{
		return s + now;
	});
	std::cout << ans << std::endl;
	std::cout << std::numeric_limits<double>::min() << std::endl;
//	std::cout << gsl_sf_psi(10) << std::endl;
//	std::cout << gsl_sf_psi_1(100) << std::endl;
	return 0;
}