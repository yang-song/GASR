template<int m,int n,int r>
void EMSampler<m, n, r>::sampleLambda(){
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

template<int m, int n, int r>
void EMSampler<m, n, r>::sampleGammas(){
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

template<int m,int n,int r>
void EMSampler<m, n, r>::EM(){
	auto s = sum(S) / EMNumber;
	auto sl = sum(SL) / EMNumber;
	auto slambda = sum(SLambda) / EMNumber;
	auto sllambda = sum(SLLambda) / EMNumber;
	double ta = a + 1.0;

	while (std::abs(ta - a) > 1e-3){
		ta = a;
		a -= (gsl_sf_psi(a) - std::log(r*a / s) - sl / r) / (gsl_sf_psi_1(a) - 1 / a);
		if (a < 0)	a = 0.1;
	}
	b = r*a / s;
	
	alpha = 0.1;
	auto talpha = alpha + 1.0;
	while (std::abs(talpha - alpha) > 1e-3){
		talpha = alpha;
		alpha -= (std::log(alpha) - gsl_sf_psi(alpha) + sllambda - std::log(slambda)) / (1 / alpha -
			gsl_sf_psi_1(alpha));
		if (alpha < 0)	alpha = 0.1;
		if (alpha > 400000) alpha = 400000;
	}
	beta = alpha / slambda;
}

template<int m, int n, int r>
void EMSampler<m, n, r>::sample(){
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