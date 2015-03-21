#include "Distributions.h"
#include <cmath>

std::random_device Distributions::rd;
std::mt19937 Distributions::eng(Distributions::rd());
const Distributions::_uniform Distributions::uniform(0, 1.0);
double Distributions::tnorm(double mu, double var, double a, double b){
	/*
	Uses direct reject sampling with uniform proposal distribution
	Specially adjusted for our task!
	*/
	double mode = (mu <= a) ? (a - mu)*(a - mu) : ((mu >= b) ? (b - mu)*(b - mu) : 0.0);
	double r = 0.0;	
	while (true){
		r = a + (b - a)*uniform(eng);
		if (uniform(eng) <= std::exp(1.0 / 2.0 / var*(mode - (r - mu)*(r - mu))))
			break;
	}
	return r;
}
double Distributions::tailNorm(double mu, double var){
	/* Especially suitable for tail norm in extreme conditions.
	* The algorithm is reject sampling based on exponential distritbuions.
	*/
	double lambda = -mu / var;
	double r = 0.0;	
	while (true){
		r = -1 / lambda*std::log(1 - uniform(eng));
		if (uniform(eng) <= std::exp(-1.0 / 2.0 / var*(r - mu)*(r - mu) + lambda*r + mu*mu / 2.0 / var))
			break;
	}
	return r;
}
double Distributions::halfNorm(double mu, double var){
	if (mu < 0){
		struct exp:public std::exception{
			const char* what() const override{
				return "Error using halfNorm!";
			}
		}error;
		throw error;
	}
	double ans = 0.0;
	auto normal = _normal(mu, std::sqrt(var));
	while (true){
		auto u = normal(eng);
		if (u >= 0){
			ans = u;
			break;
		}
	}
	if (ans == 0.0){
		struct exp :public std::exception{
			const char* what() const override{
				return "ans==0.0";
			}
		}error;
		throw error;
	}
	return ans;
}
