#ifndef SIGNAL_DATA_H
#define SIGNAL_DATA_H
#include <cmath>
#include <cstdint>
#include <vector>
#include "../Eigen/Dense"
#include <numeric>
#include <algorithm>
#include <random>

class SignalData {
public:
    SignalData();
    SignalData(std::vector<Eigen::VectorXcd> InputSignal);
    std::vector<Eigen::VectorXcd> Signals;
    std::vector<std::complex<double>> Means;
    std::vector<std::complex<double>> StdDeviation;
    void CalculateMeans();
    void CalculateDeviation();
    void PrintMeanAndDeviation();
    void GenerateRandomSignals(size_t numSignals, size_t numSamples, std::complex<double> mean, double stddev);
    Eigen::VectorXcd RuningSum(const Eigen::VectorXcd& Input);
    Eigen::VectorXcd FirstDifference(const Eigen::VectorXcd& Input);
    Eigen::VectorXcd FFT(Eigen::VectorXcd &a);
    Eigen::VectorXcd IFFT(Eigen::VectorXcd &a); 
    Eigen::VectorXcd ZeroPadPowerTwo(const Eigen::VectorXcd &a);
    Eigen::VectorXcd MovingAverage(const Eigen::VectorXcd &a,int  window_size);
private:
    Eigen::VectorXcd FFTAux(Eigen::VectorXcd a, bool invert = false);
};


#endif
