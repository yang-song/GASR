//#include"Header.h"
#include<iostream>
template<int n>
void test<n>::print(){
	for (auto item : a)
		std::cout << item << std::endl;
}