#include<vector>
#include<iostream>
#include<numeric>
template<int n>
class test{
public:
	std::vector<double> a;
	test() :a(n){
		std::iota(a.begin(), a.end(), 0);
	}
	void print();
};


template<int n>
void test<n>::print(){
	for (auto item : a)
		std::cout << item << std::endl;
}
