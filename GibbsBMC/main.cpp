#include "EMSampler.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <ctime>
#include <algorithm>
int main(int argc,char* argv[]){
	const int N = 100;
	std::stringstream name;
	name << "diary" << N << ".txt";
	std::ofstream diary(R"(..\)"+name.str());

//	std::ofstream diary("Jester2-200iter-Rank40.txt");
	//Movielens 100k
	/*
	const int m = 943, n = 1682, r = 100, iter = 100, averageNumber = 40;
	std::string trainFiles("ua.base");
	std::string testFiles("ua.test");
	*/
	//complete case
	const int m = 500, n = 200, r = 200, iter = 200, averageNumber = 50;
	name.str("");
	name << "X" << N << ".txt";
	std::string trainFiles(R"(..\)"+name.str());
	name.str("");
	name << "Z" << N << ".txt";
	std::string testFiles(R"(..\)"+name.str());
	//Movielens 1M
	/*
	const int m = 6040, n = 3952, r = 30, iter = 100, averageNumber = 40;
	std::string trainFiles("train.dat");
	std::string testFiles("test.dat");
	*/
	//Jester 1
	/*
	const int m = 24983, n = 100, r = 100, iter = 200, averageNumber = 30;
	std::string trainFiles("train1.dat");
	std::string testFiles("test1.dat");
	*/
	//Jester 2
	/*
	const int m = 23500, n = 100, r = 40, iter = 200, averageNumber = 30;
	std::string trainFiles("train2.dat");
	std::string testFiles("test2.dat");
	*/
	//Jester 3
	/*
	const int m = 24938, n = 100, r = 40, iter = 100, averageNumber = 50;
	std::string trainFiles("train3.dat");
	std::string testFiles("test3.dat");
	*/

	EMSampler<m, n, r> learner(std::max(iter - averageNumber + 1, 0), iter);
	learner.read(trainFiles);
	learner.read(testFiles, true);
	learner.init(1);
	
	auto t1 = clock();
	for (int i = 0; i < iter; i++){
		auto begin = clock();
		learner.sample();
		auto d = learner.d;
		std::sort(d.begin(), d.end());
		std::for_each(d.begin(), d.end(), [&](double n){
			std::cout << n << std::endl;
			diary << n << std::endl;
		});
		std::cout << "a = " << learner.a
			<< " b = " << learner.b << std::endl;
		std::cout << "alpha = " << learner.alpha
			<< " beta = " << learner.beta << std::endl;
		std::cout << "lambda = " << learner.lambda << std::endl;
		diary << "lambda = " << learner.lambda << std::endl;
		auto error = learner.NMAE();
		//auto error = learner.err();
		std::cout << "Iteration: " << i << ", Accuracy: " << error << std::endl;
		diary << "Iteration: " << i << ", Accuracy: " << error << std::endl;
		
		std::cout << "Elapsed time: " << (clock() - begin) / 1000.0 << "s" << std::endl;
		diary << "Elapsed time: " << (clock() - begin) / 1000.0 << "s" << std::endl;		
	}
	auto error = learner.averageNMAE();
	//auto error = learner.averageErr();
	std::cout << "Averaged Accuracy for " << averageNumber << " iterations: " <<
		error << std::endl;
	diary << "Averaged Accuracy for " << averageNumber << " iterations: " <<
		error << std::endl;
	auto d = learner.sumD;
	std::sort(d.begin(), d.end());
	std::for_each(d.begin(), d.end(), [&](double n){
		diary << n << std::endl;
	});
	diary << "Averaged Elapsed Time: " << (clock() - t1) / static_cast<double>(iter) / 1000.0 << "s"
		<< std::endl;
	return 0;
}
