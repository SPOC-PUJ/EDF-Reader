#include <string>
#include <vector>

class edf{

  public:
    edf(std::string path);
    std::vector<std::vector<int>> DataRecord;

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
    std::string NumDataRecords; // number can be -1 during recording
    std::string DurationDataRecords;
    std::string NumberSignals; // ns
  // Data Records
    std::string Label;
    std::string Transducer; // specify the applied senso
    std::string PhyisicalDimension; 
    std::string PhysicalMinimum;
    std::string PhysicalMaximum; // PMAX > PMIN, In case of a negative amplifier gain the corresponding PMAX < PMIN
    std::string DigitalMinimum;
    std::string DigitalMaximum; // PMIN != PMAX avoid division by 0
    std::string Prefiltering; // HP:0.1Hz HighPass LP:75Hz LowPass N:50Hz Notch
    std::string Nr; // sample in each data record
    std::string ReservedDos;
};
