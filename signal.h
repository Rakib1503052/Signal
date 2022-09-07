#ifndef SIGNAL_H
#define SIGNAL_H

#include <vector>
#include <array>
#include <algorithm>
#include <exception>
#include <cmath>
#include <complex>
#include <string>
#include <map>
#include <tuple>
#include "./libs/matalg.h"

namespace Signal
{
namespace DSF
{

std::vector<double> padding(std::vector<double>, size_t);

std::vector<double> padRemover(std::vector<double>, size_t);

std::vector<double> lfilter(std::vector<double> B,
                           std::vector<double> A,
                           const std::vector<double> &X,
                           std::vector<double> &Zi);

std::vector<double> filtfilt(std::vector<double> B,
                             std::vector<double> A,
                             const std::vector<double> &X);

std::vector<double> lfilter_zi(std::vector<double>, std::vector<double>);

}

namespace Frequency
{
std::map<double,double> welch(const std::vector<double> &X, const std::vector<double> &window, double fs = 1.0, size_t nperseg = NULL, size_t noverlap = NULL);
std::map<double,double> welch(const std::vector<double> &X, std::string window = "hann", double fs = 1.0, size_t nperseg = NULL, size_t noverlap = NULL);
std::vector<double> get_window(std::string, size_t);

//FFT
std::vector<std::complex<double>> FFT(const std::vector<double>&);
std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>>&);
std::vector<std::complex<double>> CZT(const std::vector<std::complex<double>>&);
std::vector<std::complex<double>> rFFT(const std::vector<double>&);
std::vector<double> rFFT_abs(const std::vector<double>&);

//invert FFT
std::vector<std::complex<double>> iFFT(const std::vector<std::complex<double>>&);
}
}
#endif // SIGNAL_H
