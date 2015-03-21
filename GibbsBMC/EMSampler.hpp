#ifndef GIBBSBMC_EMSAMPLER
#define GIBBSBMC_EMSAMPLER
#include "Sampler.hpp"
#include <vector>
//#include <gsl/gsl_sf_psi.h>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <numeric>

template<int m,int n,int r>
class EMSampler :public Sampler < m, n, r > {
protected:
	template<typename C,typename F>
	double sum(C& cont, F func){
		return std::accumulate(cont.begin(), cont.end(), 0.0, [&](double s, double v){
			return s + func(v);
		});
	}
	template<typename C>
	double sum(C& cont){
		return std::accumulate(cont.begin(), cont.end(), 0.0, [](double s, double v){
			return s + v;
		});
	}

	int EMCounter, EMNumber;
	std::vector<double> S, SL, SLLambda, SLambda;
public:
	EMSampler(int from, int end, double a = 100, double b = 100, int EMNumber=5):
		Sampler<m, n, r>(from, end, a, b), S(EMNumber), SL(EMNumber),
		SLLambda(EMNumber), SLambda(EMNumber),EMCounter(0), EMNumber(EMNumber){}

	void sampleLambda() override{
		double sum = 0.0;
		for (auto pair : omega){
			int i = pair.first;
			int j = pair.second;
			sum += (X[i][j] - Z[i][j]) * (X[i][j] - Z[i][j]);
		}
		sum /= 2.0;
		lambda = Distributions::gamrnd(alpha + observations / 2.0, sum + beta);

		SLambda[EMCounter - 1] = lambda;
		SLLambda[EMCounter - 1] = std::log(lambda);
	}

	void sampleGammas() override{
		for (int i = 0; i < r; i++)
			gamma[i] = Distributions::gamrnd(a + 1, b + d[i]);
		/*
		S[EMCounter - 1] = std::accumulate(gamma.begin(), gamma.end(),0.0,
		[](double s, double v){
		return s + v;
		});

		SL[EMCounter - 1] = std::accumulate(gamma.begin(), gamma.end(),0.0,
		[](double s, double v){
		return s + std::log(v);
		});
		*/
		S[EMCounter - 1] = sum(gamma);
		SL[EMCounter - 1] = sum(gamma, [](double n){
			return std::log(n);
		});
	}
	void EM(){
		auto s = sum(S) / EMNumber;
		auto sl = sum(SL) / EMNumber;
		auto slambda = sum(SLambda) / EMNumber;
		auto sllambda = sum(SLLambda) / EMNumber;
		double ta = a + 1.0;

		while (std::abs(ta - a) > 1e-3){
			ta = a;
			//a -= (gsl_sf_psi(a) - std::log(r*a / s) - sl / r) / (gsl_sf_psi_1(a) - 1 / a);
			a -= (boost::math::digamma(a) - std::log(r*a / s) - sl / r) / (boost::math::trigamma(a) - 1 / a);
			if (a < 0)	a = 0.1;
		}
		b = r*a / s;

		
		alpha = 0.1;
		auto talpha = alpha + 1.0;
		while (std::abs(talpha - alpha) > 1e-3){
			talpha = alpha;
			//alpha -= (std::log(alpha) - gsl_sf_psi(alpha) + sllambda - std::log(slambda)) / (1 / alpha -
			//	gsl_sf_psi_1(alpha));
			alpha -= (std::log(alpha) - boost::math::digamma(alpha) + sllambda - std::log(slambda)) / (1 / alpha -
					boost::math::trigamma(alpha));
			if (alpha < 0)	alpha = 0.1;
			if (alpha > 100000) alpha = 100000;
		}
		beta = alpha / slambda;

		/*
		alpha = 12.112589496;
		beta = alpha / slambda;
		*/
	}
	
	void sample() override{
		sampleCounter++;
		EMCounter++;

		sampleLambda();
		sampleGammas();
		sampleds();
		sampleU();
		sampleV();
		updateZ();

		if (EMCounter >= EMNumber){
			EM();
			EMCounter = 0;
		}
	}
};
#endif