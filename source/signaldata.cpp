#include "signaldata.hpp"
#include <complex>
#include <cstddef>
#include <iostream>
#include <ostream>
SignalData::SignalData() : Signals{}, Means{}, StdDeviation{} {}

SignalData::SignalData(std::vector<Eigen::VectorXcd> InputSignal)
    : Signals{InputSignal} {}

void SignalData::CalculateMeans() {
    Means.clear();
    for (const auto& signal : Signals) {
        if (signal.size() > 0) {
            std::complex<double> sum = std::accumulate(signal.begin(), signal.end(), std::complex<double>(0.0, 0.0));
            Means.push_back(sum / static_cast<double>(signal.size()));
        } else {
            Means.push_back(0.0);
        }
    }
}

void SignalData::CalculateDeviation() {
    StdDeviation.clear();
    CalculateMeans(); // Asegurarse de que los medios estén actualizados
    for (size_t i = 0; i < Signals.size(); ++i) {
        const auto& signal = Signals[i];
        auto mean = Means[i];
        if (signal.size() > 0) {
            double variance = std::accumulate(signal.begin(), signal.end(), 0.0, 
                                              [mean](double acc, const std::complex<double>& value) {
                                                  auto diff = value - mean;
                                                  return acc + std::norm(diff);
                                              }) / static_cast<double>(signal.size());
            StdDeviation.push_back(std::sqrt(variance));
        } else {
            StdDeviation.push_back(0.0);
        }
    }
}

void SignalData::PrintMeanAndDeviation(){

  if(!(Signals.size() == Means.size() == StdDeviation.size())){
    std::cout<< "Mean or standar deviation not propertly initialized" <<std::endl;
  }
  int Size = Signals.size();
  for( int i = 0 ; i < Size ; i++ ){
    std::cout<< "signal number "<< i << " mean = "<< Means[i] << " standart deviation = " << StdDeviation[i] << std::endl;
  }


}


void SignalData::GenerateRandomSignals(size_t numSignals, size_t numSamples, std::complex<double> mean, double stddev) {
    Signals.clear();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist_real(mean.real(), stddev);
    std::normal_distribution<> dist_imag(mean.imag(), stddev);

    for (size_t i = 0; i < numSignals; ++i) {
        Eigen::VectorXcd signal(numSamples);
        for (size_t j = 0; j < numSamples; ++j) {
            signal[j] = std::complex<double>(dist_real(gen), dist_imag(gen));
        }
        Signals.push_back(signal);
    }
}

std::vector<std::complex<double>> SignalData::FirstDifference(const Eigen::VectorXcd& Input){
  size_t Size = Input.size();

  std::vector<std::complex<double>> Output(Size);

  if(Input.size()<0) return Output;

  Output[0] = Input[0];

  for(int i = 1;i < Size ; i++ ){
    Output[i] = Input[i] - Input[i-1];
  }

  return Output;
}

std::vector<std::complex<double>> SignalData::RuningSum(const Eigen::VectorXcd& Input){
  size_t Size = Input.size();

  std::vector<std::complex<double>> Output(Size);

  if(Input.size()< 0) return Output;
  Output[0] = Input[0];
  for (size_t i = 1; i < Size; ++i) {
      Output[i] = Output[i - 1] + Input[i];
  }
  return Output;
}

