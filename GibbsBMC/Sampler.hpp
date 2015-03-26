#ifndef GIBBSBMC_SAMPLER
#define GIBBSBMC_SAMPLER
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <utility>
#include <string>
#include "Distributions.h"
#include <random>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <stdexcept>

template<int m,int n,int r>
class Sampler
{
public:
	template<int m, int n>
	using Matrix = std::array < std::array<double , n>, m > ;

	template<int len>
	using Vector = std::array < double, len > ;

	class hash{
	public:
		double operator()(std::pair<int, int> in)
		{
			return std::hash<int>()(in.first) ^ std::hash<int>()(in.second);
		}
	};
	// Count the No. of samples
	int sampleCounter;
	//The start and end of averaged indices
	int from, end;
	//parameters
	double a, b, alpha, beta;
	double lambda;
	//The length of the initial Us and Vs
	const double initlen;
	//The observed matrix
	Matrix<m, n> X;
	//The predicted completed matrix and matrix needed to compute averages
	Matrix<m, n> Z, sumZ;
	//The latent matrices
	Matrix<m, r> U;
	Matrix<n, r> V;
	//The squared norms of each singular vector
	Vector<r> sumU;
	Vector<r> sumV;
	//The singular value vectors
	Vector<r> d;
	Vector<r> sumD;
	//The gamma values
	Vector<r> gamma;
	//The maps to store the Omega_ij
	std::unordered_map<unsigned int, std::unordered_set<unsigned int>> row, column;
	//Train and test indices
	//std::vector<std::pair<unsigned int, unsigned int>> omega, omegaT;
	std::unordered_set <std::pair<unsigned int, unsigned int>, hash> omega, omegaT;
	int observations;
	// The true results
	std::unordered_map<std::pair<unsigned int, unsigned int>, double, hash> real;
public:
	void init(){
		std::random_device rd;
		std::mt19937 eng(rd());
		std::uniform_real_distribution<> dist(0, 1.0);
		Vector<r> initSumU{ { 0 } };
		Vector<r> initSumV{ { 0 } };
		for (int i = 0; i < m; i++)
			for (int j = 0; j < r; j++)
			{
				U[i][j] = dist(eng);
				initSumU[j] += U[i][j] * U[i][j];
			}
		for (int i = 0; i < n; i++)
			for (int j = 0; j < r; j++)
			{
				V[i][j] = dist(eng);
				initSumV[j] += V[i][j] * V[i][j];
			}
		for (int i = 0; i < m; i++)
			for (int j = 0; j < r; j++)
				U[i][j] /= sqrt(initSumU[j]) / initlen;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < r; j++)
				V[i][j] /= sqrt(initSumV[j]) / initlen;

		sumU.fill(initlen*initlen);
		sumV.fill(initlen*initlen);
	}
	void init(int iter){
		Vector<r> w = { 0 };
		std::array<std::array<bool, n>, m> flag = { { false } };
		for (auto & pair : omega)	flag[pair.first][pair.second] = true;
		Eigen::MatrixXd matX(m, n), matZ(m, n);
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
			{
				matX(i, j) = X[i][j];
				matZ(i, j) = Z[i][j];
			}
		lambda = 1.0;
		for (int c = 0; c < iter; c++){
			for (int i = 0; i < r; i++)	w[i] = static_cast<double>(a + 1.0) / (b + d[i]);
			for (int i = 0; i < m; i++)
				for (int j = 0; j < n; j++)
					if (flag[i][j])
						matZ(i, j) = matX(i, j);

			Eigen::JacobiSVD<Eigen::MatrixXd> svd(matZ,Eigen::ComputeThinU | Eigen::ComputeThinV);
			auto U = svd.matrixU();
			auto V = svd.matrixV();
			auto sv = svd.singularValues();
			for (int i = 0; i < r; i++)
				d[i] = sv(i) = std::max(0.0, sv(i) - 1 / lambda*w[i]);
			matZ = U*sv.asDiagonal()*V.transpose();
		}

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(matZ, Eigen::ComputeFullU | Eigen::ComputeFullV);
		auto svdU = svd.matrixU();
		auto svdV = svd.matrixV();
		auto svdD = svd.singularValues();
		for (int i = 0; i < m; i++)
			for (int k = 0; k < r; k++)
				U[i][k] = svdU(i, k)*initlen;

		for (int i = 0; i < n; i++)
			for (int k = 0; k < r; k++)
				V[i][k] = svdV(i, k)*initlen;

		for (int i = 0; i < r; i++)	d[i] = svdD(i);
	}
	
	Sampler(int from, int end, double a = 100, double b = 100) :
		initlen(0.9), lambda(1.0), observations(0), sampleCounter(0),
		from(from), end(end), a(a), b(b), X{ { { 0 } } }, Z{ { { 0 } } }, 
		sumZ{ { { 0 } } }, U{ { { 0 } } }, V{ { { 0 } } }, sumU{ { 0 } }, sumV{ {0} },
		d{ { 0 } }, sumD{ { 0 } }, gamma{ { 0 } }, alpha(0.0),beta(0.0){
		init();
		updateZ();
		/*
		!!!!!!!!!!!!!!!!!!!
			I REALLY NEED THIS UPDATEZ!!!
		!!!!!!!!!!!!!!!!!!!
		*/
	}
	/*
	void clear(Matrix<m, n> &m){
		for (auto& line : m)
			for (auto& item : line)
				item = 0;
	}
	template<typename T>
	void clear(T& vector){
		vector.fill(0);
	}

	Sampler(int from, int end, double a = 100, double b = 100) :
		initlen(0.9), lambda(1.0), observations(0), sampleCounter(0),
		from(from), end(end), a(a), b(b),alpha(0.0),beta(0.0){
		clear(X); clear(Z); clear(sumZ); clear(U); clear(V); clear(sumU);
		clear(sumV); clear(d); clear(sumD); clear(gamma);
		init();
	}
	*/	
	void read(std::string file, bool test = false){
		unsigned int user, movie;
		double rate;
		std::ifstream fin(file);

		if (!test)
		{
			while (fin >> user >> movie >> rate){
				fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				user--;
				movie--;
				X[user][movie] = rate;
				//omega.push_back(std::make_pair(user, movie));
				omega.emplace(std::make_pair(user, movie));
				if (row.find(user) == row.end()){
					std::unordered_set<unsigned int> tmp{ movie };
					row.emplace(std::make_pair(user, tmp));
				}
				else
					row[user].emplace(movie);
				if (column.find(movie) == column.end()){
					std::unordered_set < unsigned int > tmp{ user };
					column.emplace(std::make_pair(movie, tmp));
				}
				else
					column[movie].emplace(user);
			}
		}
		else
		{
			while (fin >> user >> movie >> rate){
				fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				user--;
				movie--;
				if (real.find(std::make_pair(user, movie)) == real.end()){
					real.emplace(std::make_pair(std::make_pair(user, movie), rate));
					//omegaT.push_back(std::make_pair(user, movie));
					omegaT.emplace(std::make_pair(user, movie));
				}
				else
					throw std::runtime_error("Error! Multiple test key-value pairs!");

			}
		}
		observations = omega.size();
	}
	virtual void sampleLambda(){
		double sum = 0.0;
		for (auto pair : omega){
			int i = pair.first;
			int j = pair.second;
			sum += (X[i][j] - Z[i][j]) * (X[i][j] - Z[i][j]);
		}
		sum /= 2.0;
		lambda = Distributions::gamrnd(alpha + observations / 2.0, sum + beta);
	}
	virtual void sampleGammas(){
		for (int i = 0; i < r; i++)
			gamma[i] = Distributions::gamrnd(a + 1, b + d[i]);
	}
	void sampleds(){
		for (int a = 0; a < r; a++){
			double A = std::accumulate(omega.begin(), omega.end(), 0.0,
				[&](double s, std::pair<unsigned int, unsigned int> pair){
				return s + U[pair.first][a] * U[pair.first][a] * V[pair.second][a] * V[pair.second][a];
			});
			double B = std::accumulate(omega.begin(), omega.end(), 0.0,
				[&](double s, std::pair<unsigned int, unsigned int> pair){
				double sum = 0.0;
				for (int k = 0; k < r; k++)
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
	void sampleU(){
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
					auto C = std::accumulate(set.begin(), set.end(), 0.0, [&](double s, unsigned int j){
						return s + d[a] * d[a] * V[j][a] * V[j][a];
					});
					auto D = std::accumulate(set.begin(), set.end(), 0.0, [&](double s, unsigned int j){
						double sum = 0.0;
						for (int k = 0; k != r; k++)
							if (k != a)
								sum += d[a] * d[k] * U[b][k] * V[j][a] * V[j][k];
						return s + sum - d[a] * X[b][j] * V[j][a];
					});
					U[b][a] = Distributions::tnorm(-D / C, 1 / C / lambda, -bound, bound);
				}
				sumU[a] += U[b][a] * U[b][a];
			}
	}
	void sampleV(){
		for (int a = 0; a < r; a++)
			for (int b = 0; b < n; b++){
				sumV[a] -= V[b][a] * V[b][a];
				auto bound = std::sqrt(1 - sumV[a]);
				if (column.find(b) == column.end())
					V[b][a] = Distributions::_uniform(-bound, bound)(Distributions::eng);
				else{
					auto &set = column[b];
					auto E = std::accumulate(set.begin(), set.end(), 0.0, [&](double s, unsigned int i){
						return s + d[a] * d[a] * U[i][a] * U[i][a];
					});
					auto F = std::accumulate(set.begin(), set.end(), 0.0, [&](double s, unsigned int i){
						double sum = 0.0;
						for (int k = 0; k != r; k++)
							if (k != a)
								sum += d[a] * d[k] * V[b][k] * U[i][a] * U[i][k];
						return s + sum - d[a] * X[i][b] * U[i][a];
					});
					V[b][a] = Distributions::tnorm(-F / E, 1 / E / lambda, -bound, bound);
				}
				sumV[a] += V[b][a] * V[b][a];
			}
	}
	virtual void updateZ(){
		//Compute the latent matrix Z
		for (auto &line : Z)
			line.fill(0.0);
		/*
		//Why is this wrong! Repeatitions!!!!
		auto needed = omega;
		needed.insert(needed.end(), omegaT.begin(), omegaT.end());
		for (auto& pair : needed){
			int i = pair.first;
			int j = pair.second;
			for (int k = 0; k < r; k++)
				Z[i][j] += d[k] * U[i][k] * V[j][k];
			if (sampleCounter <= end && sampleCounter >= from)
				sumZ[i][j] += Z[i][j];
		}
		*/
		
		// This proves to be right!!
		/*
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < r; k++)
					Z[i][j] += d[k] * U[i][k] * V[j][k];
				if (sampleCounter <= end && sampleCounter >= from)
					sumZ[i][j] += Z[i][j];
			}
		*/	
		
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
	}
	virtual void sample(){
		sampleCounter++;
		sampleLambda();
		sampleGammas();
		sampleds();
		sampleU();
		sampleV();
		updateZ();
	}
	double _RMSE(Matrix<m, n> &Z){
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
	double RMSE(){ 
		return _RMSE(Z); 
	}
	double averageRMSE(){
		for (int i = 0; i < r; i++)	sumD[i] /= (end - from + 1);
		auto ans = sumZ;
		for (auto &line : ans)
			for (auto &item : line)
				item /= double(end - from + 1);
		return _RMSE(ans);
	}
	double _err(Matrix<m, n> &Z){
		double sum = 0.0;
		double tmp = 0.0;
		for (auto &record : real){
			int i = record.first.first;
			int j = record.first.second;
			double rate = record.second;
			auto dis = Z[i][j] - rate;
			sum += dis*dis;
			tmp += rate*rate;
		}
		return sum / tmp;
	}
	double err() { 
		return _err(Z); 
	}
	double averageErr(){
		for (auto& num : sumD)
			num /= (end - from + 1);
		auto ans = sumZ;
		for (auto& line : ans)
			for (auto& item : line)
				item /= (end - from + 1);
		return _err(ans);
	}
	double _NMAE(Matrix<m, n> &Z){
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
	double NMAE() { 
		return _NMAE(Z); 
	}
	double averageNMAE(){
		for (auto& num : sumD)
			num /= (end - from + 1);
		Matrix<m, n> ans = sumZ;
		for (auto& line : ans)
			for (auto& item : line)
				item /= (end - from + 1);
		return _NMAE(ans);
	}
};

#endif