#include "Distributions.h"
#include <random>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <utility>
#include <limits>

template<int m,int n,int r>
void Sampler<m ,n ,r>::init()
{
	std::random_device rd;
	std::mt19937 eng(rd());
	std::uniform_real_distribution<> dist(0, 1.0);

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			U[i][j] = dist(eng);
			V[i][j] = dist(eng);
			sumU[j] += U[i][j] * U[i][j];
			sumV[j] += V[i][j] * V[i][j];
		}
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			U[i][j] /= sqrt(sumU[j]) / initlen;
			V[i][j] /= sqrt(sumV[j]) / initlen;
		}
	sumU.fill(initlen*initlen);
	sumV.fill(initlen*initlen);
}

template<int m,int n,int r>
void Sampler<m, n, r>::read(std::string file, bool test = false)
{
	int user, movie;
	double rate;
	ifstream fin(file);

	if (!test)
	{
		fin >> user >> movie >> rate;
		user--;
		movie--;
		X[user][movie] = rate;
		omega.push_back(std::make_pair(user, movie));
		if (row.find(user) == row.end())
			row.emplace(std::make_pair(user,std::unordered_set<unsigned int> tmp{ movie }));
		else
			row[user].emplace(movie);
		if (column.find(movie) == column.end())
			column.emplace(std::make_pair(movie, std::unordered_set<unsigned int> tmp{ user }));
		else
			column[movie].emplace(user);
	}
	else
	{
		fin >> user >> movie >> rate;
		if (real.find(std::make_pair(user, movie)) == real.end())
			real.emplace(std::make_pair(std::make_pair(user, movie), rate));
		else{
			struct exp :public std::exception{
				const char* what() const override{
					return "Error! Multiple test key-value pairs!";
				}
			}error;

			throw error;
		}
	}
	observations = omega.size();
}

template<int m,int n,int r>
void Sampler<m,n,r>::init(int iter) {  }

template<int m,int n,int r>
void Sampler<m, n, r>::sampleLambda(){
	double sum = 0.0;
	for (auto pair : omega){
		int i = pair.first;
		int j = pair.second;
		sum += (X[i][j] - Z[i][j]) * (X[i][j] - Z[i][j]);
	}
	sum /= 2.0;
	lambda = Distributions::gamrnd(alpha + observations / 2.0, sum + beta);
}

template<int m,int n,int r>
void Sampler<m, n, r>::sampleGammas(){
	for (int i = 0; i < r; i++)
		gamma[i] = Distributions::gamrnd(a + 1, b + d[i]);
}

template<int m,int n,int r>
void Sampler<m, n, r>::sampleds(){
	for (int a = 0; a < r; a++){
		double A = std::accumulate(omega.begin(), omega.end(), 0.0, 
			[](double s, std::pair<unsigned int, unsigned int> pair){
			return s + U[pair.first][a] * U[pair.first][a] * V[pair.second][a] * V[pair.second][a];
		});
		double B = std::accumulate(omega.begin(), omega.end(),0.0,
			[](double s, std::pair<unsigned int, unsigned int> pair){
			double sum = 0.0;
			for (k = 0; k < r;k++)
				if (k != a){
					sum += d[k] * U[pair.first][a] * U[pair.first][k] * 
						V[pair.second][a] * V[pair.second][k];
				}
			return s + sum - X[pair.first][pair.second] * U[pair.first][a] * V[pair.second][a];
		}) + gamma[a] / lambda;

		if (-B / A < 0)
			d[a] = Distributions::tailNorm(-B / A, 1 / A / lambda);
		else
			d[a] = Distributions::halfNorm(-B / A, 1 / A / lambda);
	}
}

template<int m,int n,int r>
void Sampler<m, n, r>::sampleU(){
	/*
	!!!!!!!!!!!!!!!!!!
	Mind the danger that some b never exists in the training data!
	Then U obey the prior distribution, e.g. truncated uniform distribution.
	*/
	for (int a = 0; a < r; a++)
		for (int b = 0; b < m; b++){
			sumU[a] -= U[b][a] * U[b][a];
			auto bound = std::sqrt(1 - sumU[a]);
			if (row.find(b) == row.end())
				U[b][a] = Distributions::_uniform(-bound, bound)(Distributions::eng);
			else{
				auto &set = row[b];
				auto C = std::accumulate(set.begin(), set.end(), 0.0, [](double s, unsigned int j){
					return s + d[a] * d[a] * V[j][a] * V[j][a];
				});
				auto D = std::accumulate(set.begin(), set.end(), 0.0, [](double s, unsigned int j){
					double sum = 0.0;
					for (k = 0; k != r; k++)
						if (k != a)
							sum += d[a] * d[k] * U[b][k] * V[j][a] * V[j][k];
					return s + sum - d[a] * X[b][j] * V[j][a];
				});
				U[b][a] = tnorm(-D / C, 1 / C / lambda, -bound, bound);
			}
			sumU[a] += U[b][a] * U[b][a];
		}
}

template<int m,int n,int r>
void Sampler<m, n, r>::sampleV(){
	for (int a = 0; a < r; a++)
		for (int b = 0; b < n; b++){
			sumV[a] -= V[b][a] * V[b][a];
			auto bound = std::sqrt(1 - sumV[a]);
			if (column.find(b) == column.end())
				V[b][a] = Distributions::_uniform(-bound, bound)(Distributions::eng);
			else{
				auto &set = column[b];
				auto E = std::accumulate(set.begin(), set.end(), 0.0, [](double s, unsigned int i){
					return s + d[a] * d[a] * U[i][a] * U[i][a];
				});
				auto F = std::accumulate(set.begin(), set.end(), 0.0, [](double s, unsigned int i){
					double sum = 0.0;
					for (k = 0; k != r; k++)
						if (k != a)
							sum += d[a] * d[k] * V[b][k] * U[i][a] * U[i][k];
					return s + sum - d[a] * X[i][b] * U[i][a];
				});
				V[b][a] = tnorm(-F / E, 1 / E / lambda, -bound, bound);
			}
			sumV[a] += V[b][a] * V[b][a];
		}
}

template<int m,int n,int r>
void Sampler<m, n, r>::updateZ(){
	//Compute the latent matrix Z
	for (auto &line : Z)
		line.fill(0.0);
	auto needed = omega;
	needed.insert(needed.end(), omegaT.begin(), omegaT.end());
	for (auto& pair : needed){
		int i = pair.first;
		int j = pair.second;
		for (k = 0; k < r; k++)
			Z[i][j] += d[k] * U[i][k] * V[j][k];
		if (sampleCounter <= end && sampleCounter >= from)
			sumZ[i][j] += Z[i][j];
	}
	if (sampleCounter <= end && sampleCounter >= from)
		for (i = 0; i < r; i++)	sumD[i] += d[i];
}

template<int m,int n,int r>
double Sampler<m, n, r>::_RMSE(Sampler<m,n,r>::Matrix<m,n>& Z){
	//must run updateZ() first.
	//Calculate the predict accuracy.
	double sum = 0.0;
	int N = 0;
	for (auto &record : real){
		int i = record.first.first;
		int j = record.first.second;
		double rate = record.second;
		double dis = Z[i][j] - rate;
		sum += dis*dis;
	}
	return std::sqrt(sum / N);
}

template<int m,int n,int r>
double Sampler<m, n, r>::averageRMSE(){
	for (i = 0; i < r; i++)	sumD[i] /= (end - from + 1);
	Matrix<m, n> ans = sumZ;
	for (auto &line : ans)
		for (auto &item : line)
			item /= double(end - from + 1);
	return _RMSE(ans);
}

template<int m,int n,int r>
double Sampler<m, n, r>::_NMAE(Sampler<m, n, r>::Matrix<m, n> & Z){
	int N = 0;
	double sum = 0.0;
	double maxi = -std::numeric_limits<double>::min();
	double mini = std::numeric_limits<double>::max();
	for (auto &record : real){
		int i = record.first.first;
		int j = record.first.second;
		double rate = record.second;
		N++;
		sum += std::abs(rate - Z[i][j]);
		maxi = std::max(rate, maxi);
		mini = std::min(rate, mini);
	}
	return sum / N / (maxi - mini);
}

template<int m,int n,int r>
double Sampler<m, n, r>::averageNMAE(){
	for (auto& num : sumD)
		num /= (end - from + 1);
	Matrix<m, n> ans = sumZ;
	for (auto& line : ans)
		for (auto& item : line)
			item /= (end - from + 1);
	return _NMAE(ans);
}

template<int m,int n,int r>
void Sampler<m, n, r>::sample(){
	sampleCounter++;
	sampleLambda();
	sampleGammas();
	sampleds();
	sampleU();
	sampleV();
	updateZ();
}