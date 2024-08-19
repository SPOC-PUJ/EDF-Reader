#include <emscripten/bind.h>
#include "edfwasm.hpp"
#include "signaldata.hpp"
#include "../Eigen/Dense"
EMSCRIPTEN_BINDINGS(edf_module) {
    emscripten::register_vector<Eigen::VectorXcd>("VectorXcdVector");
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

        emscripten::class_<Eigen::VectorXcd>("EigenVector");
        
}
