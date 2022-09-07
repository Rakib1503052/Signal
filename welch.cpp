#include "signal.h"

using namespace std;
using dvec = vector<double>;

double cumul(const vector<double>& X)
{
	double sum = 0.0;
	for (auto& n : X) { sum += n; }

	return sum;
}

vector<vector<double>> fft_helper(const vector<double> &X, const vector<double> &win, size_t nperseg, size_t noverlap)
{
	//Helper function. Should not be called outside.	//BTW, the warning is from scipy.
	size_t seg_num = floor((X.size() - nperseg) / (nperseg - noverlap)) + 1;
	vector<vector<double>> segments;
	segments.reserve(seg_num);
	size_t i = 0;

	//Make segments.
	while ((i + nperseg) <= X.size())
	{
		segments.push_back(vector<double>(X.begin() + i, X.begin() + i + nperseg));
		i = i + nperseg - noverlap;
	}

	//Apply window.
	for (auto& n : segments)
	{
		for (size_t i = 0; i < n.size(); i++)
			n[i] = n[i] * win[i];
	}

	for (auto& n : segments)
	{
		n = Signal::Frequency::rFFT_abs(n);

		for (auto& m : n) { m = pow(m, 2); }
	}

	return segments;
}

pair<dvec, vector<dvec>> spectral_helper(const vector<double> &X, const vector<double> &Y, const vector<double> &win, double fs = 1.0,
                                           size_t nperseg = NULL, size_t noverlap = NULL)
{
	//Private function. Used by welch. Should not be called from outside of this library.

	if (nperseg == NULL)
	{
		nperseg = win.size();
	}
	if (nperseg != win.size())
	{
		throw runtime_error("Error: given window size does not match 'nperseg' value.");
	}

	//Window and nperseg set.

	if (noverlap == NULL) { noverlap = floor(nperseg / 2); }
	if (noverlap >= nperseg) { throw runtime_error("Are you dumb? You want bigger overlap than segment size? SERIOUSLY!!!"); }  //Should change that.
																																  //On second thought , Nah!!!
//Determine scaling.
	double sum = 0;
	for (const auto& n : win) { sum += (n * n); }
	double scale = 1.0 / (fs * sum);

	//Determine sampling frequencies.
	double val = fs / nperseg;
	size_t N = size_t(nperseg / 2) + 1;
	vector<double> freqs(N, 0);
	for (size_t i = 0; i < N; i++) { freqs[i] = i * val; }

	//FFT section
	vector<vector<double>> result = fft_helper(X, win, nperseg, noverlap);
	for (auto& n : result) { for (auto& m : n) { m *= scale; } }

	if ((nperseg % 2) == 1) {
		for (auto& n : result) {
			for (size_t i = 1; i < n.size(); i++) { n[i] *= 2; }
		}
	}
	else {
		for (auto& n : result) {
			for (size_t i = 1; i < (n.size() - 1); i++) { n[i] *= 2; }
		}
	}

	//return make_tuple(freqs, time, switch_axis(result));
	return make_pair(freqs, matalg::mat_T(result));
}

map<double, double> Signal::Frequency::welch(const vector<double> &X, const vector<double> &window, double fs, size_t nperseg, size_t noverlap)
{
	map<double, double> result;
	auto [freq, Pxy] = spectral_helper(X, X, window, fs, nperseg, noverlap);

	for (size_t i = 0; i < freq.size(); i++)
	{
		result[freq[i]] = cumul(Pxy[i]) / double(Pxy[i].size());
	}

	return result;
}

map<double, double> Signal::Frequency::welch(const vector<double> &X, string window, double fs, size_t nperseg, size_t noverlap)
{
	if (nperseg == NULL) { nperseg = 256; }
	vector<double> win = Signal::Frequency::get_window(window, nperseg);

	return Signal::Frequency::welch(X, win, fs, nperseg, noverlap);
}
