#define _USE_MATH_DEFINES
#include "signal.h"
#include <future>

using namespace std;
using dcvec = vector<complex<double>>;

vector<complex<double>> Signal::Frequency::FFT(const vector<complex<double>>& X)
{
    //This thing gave me a headache. If you don't want to be cursed, don't touch it.
    //This one function alone is half of the reason I had to compile in GCC instead of MSVC.

	using namespace std::complex_literals;
	size_t N = X.size();

	if (N == 0) { throw runtime_error("Empty data"); }

	if (N == 1) { return X; }

	else
	{
		if ((N % 2) != 0)
		{
			return Signal::Frequency::CZT(X);
		}
		else
		{
			vector<complex<double>> y1(N / 2), y2(N / 2);
			for (size_t n = 0; n < (N / 2); n++)
			{
				y1[n] = X[2 * n];
				y2[n] = X[2 * n + 1];
			}
			y1 = Signal::Frequency::FFT(y1);
			y2 = Signal::Frequency::FFT(y2);

			vector<complex<double>> term1(N / 2), term2(N / 2);
			complex<double> W_N = 1.0;
			complex<double> fact = exp(-2.0i * M_PI / double(N));

			term1[0] = y1[0] + y2[0];
            term2[0] = y1[0] - y2[0];

			for (size_t n = 1; n < N / 2; n++)
			{
				W_N = W_N * fact;
				term1[n] = y1[n] + W_N * y2[n];
				term2[n] = y1[n] + (-W_N) * y2[n];
			}
			term1.insert(term1.end(), term2.begin(), term2.end());		//add vector name //I forgot why I wrote this comment.
																	   //If there is any error, this line is probably the reason.

			return term1;
		}
	}
}

//Same process as GNU Octave's "czt.m"
vector<complex<double>> Signal::Frequency::CZT(const vector<complex<double>>& X)
{
	//Don't ask me what Chirp-z transform is. Just Google it. I have absolutely zero idea.
	//I don't know what it is. I don't know how it is. I don't know why it is.
	//Just don't touch it. DON'T DARE TO TOUCH THIS CODE. Don't even put your cursor here.
	//Because this is all magic and YOU DON'T MESS WITH MAGIC.

	using namespace std::complex_literals;
	size_t n = X.size();
	size_t N = 2 * n - 1;
	complex<double> w = exp(-2.0i * M_PI / double(n));
	double a = 1.0;
	double m = 1 - double(n);

	vector<complex<double>> chirp(N);
	for (size_t k = 0; k < N; k++, m++) { chirp[k] = pow(w, (pow(m, 2) / 2.0)); }

	size_t N2 = pow(2, ceil(log2(N)));

	vector<complex<double>> xp(n);
	for (size_t k = 0; k < n; k++) { xp[k] = X[k] * pow(a, -int(k)) * chirp[k + n - 1]; }
	xp.resize(N2, complex<double>(0, 0));

	vector<complex<double>> ichirp(N);
	transform(chirp.begin(), chirp.end(), ichirp.begin(), [](const complex<double>& k)->complex<double> {return (1.0 / k); });
	ichirp.resize(N2, complex<double>(0, 0));

	auto fut_xp = async(launch::deferred, (dcvec(*) (const dcvec&)) &Signal::Frequency::FFT, xp);
	auto fut_ichirp = async(launch::deferred, (dcvec(*) (const dcvec&)) &Signal::Frequency::FFT, ichirp);
	xp = fut_xp.get();
	ichirp = fut_ichirp.get();
	for (size_t k = 0; k < N2; k++) { xp[k] = xp[k] * ichirp[k]; }
	xp = Signal::Frequency::iFFT(xp);

	vector<complex<double>> ret(n);
	for (size_t k = 0; k < n; k++) { ret[k] = xp[k + n - 1] * chirp[k + n - 1]; }

	return ret;
}

vector<complex<double>> Signal::Frequency::FFT(const vector<double>& X)
{
	size_t n = X.size();
	vector<complex<double>> Y(n);
	for (size_t i = 0; i < n; i++)
	{
		Y[i] = complex<double>(X[i], 0);
	}

	return Signal::Frequency::FFT(Y);
}

vector<complex<double>> Signal::Frequency::iFFT(const vector<complex<double>>& X)
{
	size_t N = X.size();
	double dN = double(N);
	vector<complex<double>> ret = X;
	for (auto& n : ret) { n = conj(n); }
	ret = Signal::Frequency::FFT(ret);
	for (auto& n : ret) { n = conj(n) / dN; }

	return ret;
}

vector<complex<double>> Signal::Frequency::rFFT(const vector<double>& X)
{
	vector<complex<double>> ret = Signal::Frequency::FFT(X);
	size_t N = ret.size();

	if (N % 2 == 0) { ret.resize(N / 2 + 1); }
	else { ret.resize((N + 1) / 2); }

	return ret;
}

vector<double> Signal::Frequency::rFFT_abs(const vector<double>& X)
{
	vector<complex<double>> rfft_complex = Signal::Frequency::rFFT(X);
	size_t n = rfft_complex.size();
	vector<double> ret(n);

	for (size_t i = 0; i < n; i++) { ret[i] = abs(rfft_complex[i]); }

	return ret;
}
