#include <emscripten/bind.h>
#include "edfwasm.hpp"
#include "signaldata.hpp"
#include "../Eigen/Dense"
#include <complex>

double get_real(const std::complex<double>& c) {
    return c.real();
}

double get_imag(const std::complex<double>& c) {
    return c.imag();
}
std::complex<double> getElement(const Eigen::VectorXcd& vec, int index) {
    return vec(index);
}

EMSCRIPTEN_BINDINGS(edf_module) {
    emscripten::register_vector<Eigen::VectorXcd>("VectorXcdVector");
    // emscripten::class_<std::complex<double>>("Complex");
    emscripten::class_<edf>("EDF")
        .constructor<const std::string&>()
        .function("PrintHeaderRecords", &edf::PrintHeaderRecords)
        .function("PrintDataRecords", &edf::PrintDataRecords)
        .function("PrintSizeSignals", &edf::PrintSizeSignals)
        .function("PrintTopValues", &edf::PrintTopValues)
        .property("Signals", &edf::Signals);

        emscripten::class_<SignalData>("SignalData")
        .constructor<>()
        .function("CalculateMeans", &SignalData::CalculateMeans)
        .function("CalculateDeviation", &SignalData::CalculateDeviation)
        .function("PrintMeanAndDeviation", &SignalData::PrintMeanAndDeviation)
        .property("signals", &SignalData::Signals);

        emscripten::class_<Eigen::VectorXcd>("EigenVector")
        .property("size", &Eigen::VectorXcd::size)
        .function("get", &getElement);

        emscripten::class_<std::complex<double>>("Complex")
        .function("real", &get_real)
        .function("imag", &get_imag);
        
}
