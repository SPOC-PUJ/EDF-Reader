#include "signaldata.hpp"
#include <algorithm>
#include <iostream>
#include <ostream>
SignalData::SignalData() : Signals{}, Means{}, StdDeviation{} {}

SignalData::SignalData(std::vector<std::vector<int16_t>> InputSignal)
    : Signals{InputSignal} {}

void SignalData::CalculateMeans() {
    Means.clear();
    for (const auto& signal : Signals) {
        if (!signal.empty()) {
            float mean = std::accumulate(signal.begin(), signal.end(), 0.0) / signal.size();
            Means.push_back(mean);
        } else {
            Means.push_back(0.0f); // Manejo de señales vacías
        }
    }
}

void SignalData::CalculateDeviation() {
    StdDeviation.clear();
    CalculateMeans(); // Asegurarse de que los medios estén actualizados
    for (size_t i = 0; i < Signals.size(); ++i) {
        const auto& signal = Signals[i];
        float mean = Means[i];
        if (!signal.empty()) {
            float variance = std::accumulate(signal.begin(), signal.end(), 0.0,[mean](float acc, int16_t value) {
                 return acc + (value - mean) * (value - mean);
             }) / signal.size();
            StdDeviation.push_back(std::sqrt(variance));
        } else {
            StdDeviation.push_back(0.0f); // Manejo de señales vacías
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


void SignalData::GenerateRandomSignals(size_t numSignals, size_t numSamples, float mean, float stddev) {
    Signals.clear();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(mean, stddev);

    for (size_t i = 0; i < numSignals; ++i) {
        std::vector<int16_t> signal;
        for (size_t j = 0; j < numSamples; ++j) {
            signal.push_back(static_cast<int16_t>(dist(gen)));
        }
        Signals.push_back(signal);
    }
}
