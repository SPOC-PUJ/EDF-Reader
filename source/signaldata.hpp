#ifndef SIGNAL_DATA_H
#define SIGNAL_DATA_H
#include <cmath>
#include <cstdint>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>

class SignalData {
public:
    SignalData();
    SignalData(std::vector<std::vector<int16_t>> InputSignal);
    std::vector<std::vector<int16_t>> Signals;
    std::vector<float> Means;
    std::vector<float> StdDeviation;
    void CalculateMeans();
    void CalculateDeviation();
    void PrintMeanAndDeviation();
    void GenerateRandomSignals(size_t numSignals, size_t numSamples, float mean, float stddev);
    std::vector<int16_t> RuningSum(const std::vector<int16_t>& Input);
    std::vector<int16_t> FirstDifference(const std::vector<int16_t>& Input);
private:

};


#endif
