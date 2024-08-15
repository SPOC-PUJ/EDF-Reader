// edf.hpp

#ifndef EDF_H
#define EDF_H
#include "signaldata.hpp"
#include <string>
#include <vector>

struct DataRecords {
    std::string Label;
    std::string Transducer;
    std::string PhyisicalDimension;
    std::string PhysicalMinimum;
    std::string PhysicalMaximum;
    std::string DigitalMinimum;
    std::string DigitalMaximum;
    std::string Prefiltering;
    int Nr;
    std::string ReservedDos;
};

class edf {
  public:
    edf( const std::string& fileData); // Modify constructor
    SignalData Signals;
    void PrintHeaderRecords();
    void PrintDataRecords();
    void PrintSizeSignals();
    void PrintTopValues(int i);

  private:
    void parseEDFData(const std::string& data); // New method to parse data

    std::string Version;
    std::string PatientId;
    std::string RecordingId;
    std::string StartDate;
    std::string StartTime;
    std::string SizeHeader;
    std::string Reserved;
    int NumDataRecords;
    std::string DurationDataRecords;
    int NumberSignals;
    std::vector<DataRecords> SignalsInfo;
};

#endif
