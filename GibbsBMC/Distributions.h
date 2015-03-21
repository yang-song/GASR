#ifndef GIBBSBMC_DISTRIBUTIONS
#define GIBBSBMC_DISTRIBUTIONS
#include<random>
class Distributions{
public:
	static std::random_device rd;
	static std::mt19937 eng;
	
	using _uniform = std::uniform_real_distribution < > ;
	using _gamma = std::gamma_distribution < > ;
	using _normal = std::normal_distribution < > ;

	static const _uniform uniform;
public:
	static double gamrnd(double a, double b){
		return _gamma(a, 1.0 / b)(eng);
	}
	static double tnorm(double mu, double var, double a, double b);
	static double tailNorm(double mu, double var);
	static double halfNorm(double mu, double var);
};
#endif