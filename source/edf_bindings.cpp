#include <emscripten/bind.h>
#include "edfwasm.hpp"
#include <vector>
EMSCRIPTEN_BINDINGS(edf_module) {
    
    emscripten::class_<edf>("EDF")
        .constructor<const std::string&>()
        .function("PrintHeaderRecords", &edf::PrintHeaderRecords)
        .function("PrintDataRecords", &edf::PrintDataRecords)
        .function("PrintSizeSignals", &edf::PrintSizeSignals)
        .function("PrintTopValues", &edf::PrintTopValues);
}
