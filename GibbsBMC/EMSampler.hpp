#ifndef GIBBSBMC_EMSAMPLER
#define GIBBSBMC_EMSAMPLER
#include "Sampler.hpp"
#include <vector>
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
	std::vector<double> S, SL, SLLambda, SLambda,se;
public:	
	using Sampler < m, n, r > ::lambda;
	using Sampler<m, n, r>::observations;
	using Sampler<m, n, r>::sampleCounter;
	using Sampler<m, n, r>::a;
	using Sampler<m, n, r>::b;
	using Sampler<m, n, r>::X;
	using Sampler<m, n, r>::sumZ;
	using Sampler<m, n, r>::U;
	using Sampler<m, n, r>::V;
	using Sampler<m, n, r>::d;
	using Sampler<m, n, r>::gamma;
	EMSampler(int from, int end, double a = 100, double b = 100, int EMNumber=5):
		Sampler<m, n, r>(from, end, a, b), S(EMNumber), SL(EMNumber),
		SLLambda(EMNumber), SLambda(EMNumber),EMCounter(0), EMNumber(EMNumber),se(EMNumber){}

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
//		auto slambda = sum(SLambda) / EMNumber;
//		auto sllambda = sum(SLLambda) / EMNumber;
		double ta = a + 1.0;

		while (std::abs(ta - a) > 1e-3){
			ta = a;
			a -= (boost::math::digamma(a) - std::log(r*a / s) - sl / r) / (boost::math::trigamma(a) - 1 / a);
			if (a < 0)	a = 0.1;
		}
		b = r*a / s;

		double mse = sum(se) / EMNumber;
		lambda = observations / mse;
	}
	void updateZ(){
	//Compute the latent matrix Z
		for (auto &line : Z)
			line.fill(0.0);

		auto needed = omega;
		needed.insert(omegaT.begin(), omegaT.end());
		for (auto& pair : needed){
			int i = pair.first;
			int j = pair.second;
			for (int k = 0; k < r; k++)
				Z[i][j] += d[k] * U[i][k] * V[j][k];
				if (sampleCounter <= end && sampleCounter >= from)
					sumZ[i][j] += Z[i][j];
			}
		
		if (sampleCounter <= end && sampleCounter >= from)
			for (int i = 0; i < r; i++)	sumD[i] += d[i];
		se[EMCounter - 1] = 0.0;
		for (auto & item: omega){
			int i = item.first;
			int j = item.second;
			se[EMCounter - 1] += (X[i][j] - Z[i][j])*(X[i][j] - Z[i][j]);
		}
	}
	void sample() override{
		sampleCounter++;
		EMCounter++;

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