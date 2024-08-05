#ifndef SIGNAL_DATA_H
#define SIGNAL_DATA_H
#include <cmath>
#include <cstdint>
#include <vector>
#include <numeric>
#include <algorithm>

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
private:

};


#endif
