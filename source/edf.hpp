#ifndef EDF_H
#define EDF_H
#include "signaldata.hpp"
#include <string>
#include <vector>


struct DataRecords{
    std::string Label;
    std::string Transducer; // specify the applied senso
    std::string PhyisicalDimension; 
    std::string PhysicalMinimum;
    std::string PhysicalMaximum; // PMAX > PMIN, In case of a negative amplifier gain the corresponding PMAX < PMIN
    std::string DigitalMinimum;
    std::string DigitalMaximum; // PMIN != PMAX avoid division by 0
    std::string Prefiltering; // HP:0.1Hz HighPass LP:75Hz LowPass N:50Hz Notch
    int Nr; // sample in each data record
    std::string ReservedDos;
};


class edf{

  public:
    edf(const std::string path);
    SignalData Signals;
    void PrintHeaderRecords();
    void PrintDataRecords();
    void PrintSizeSignals();
    void WriteRawCsv(const std::string filename);
    void PrintTopValues(int i);
  private:
  // 61440 Max size for EDF file
  // Header Records
    std::string Version;
    std::string PatientId; //Structure: code, sex, birthdate, name
    std::string RecordingId; //Structure: text = 'Startdate' ,StartDate ,hospital administration , code for investigator, code for equipment,  
    std::string StartDate;
    std::string StartTime;
    std::string SizeHeader;
    std::string Reserved; // EDF+C continous , EDF+D discontinous
    int NumDataRecords; // number can be -1 during recording
    std::string DurationDataRecords;
    int NumberSignals; // ns
    // Data Records
    std::vector<DataRecords> SignalsInfo;
};


#endif
