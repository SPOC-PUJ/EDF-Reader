#include "edf.hpp"
#include <complex>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>



edf::edf(const std::string path){
  std::ifstream is(path , std::ios::binary);
  if(!is.is_open()){
    std::cerr << "Error: Could not open the file " << path << std::endl;  
    throw std::runtime_error("Faile to open file");
  }

  std::vector<char> buffer(80); // Single buffer for reading

  // Read Version
  is.read(buffer.data(), 8);
  if (!is) {
      std::cerr << "Error: Could not read the file (Version) " << path << std::endl;
      throw std::runtime_error("Failed to read file");
  }
  Version = std::string(buffer.data(), 8);

  // Read PatientId
  is.read(buffer.data(), 80);
  if (!is) {
      std::cerr << "Error: Could not read the file (PatientId) " << path << std::endl;
      throw std::runtime_error("Failed to read file");
  }
  PatientId = std::string(buffer.data(), 80);

  // Read RecordingID
  is.read(buffer.data(), 80);
  if (!is) {
      std::cerr << "Error: Could not read the file (RecordingID) " << path << std::endl;
      throw std::runtime_error("Failed to read file");
  }
  RecordingId = std::string(buffer.data(), 80);

  // Read StartDate
  is.read(buffer.data(), 8);
  if (!is) {
      std::cerr << "Error: Could not read the file (StartDate) " << path << std::endl;
      throw std::runtime_error("Failed to read file");
  }
  StartDate = std::string(buffer.data(), 8);


  // Read StartDate
  is.read(buffer.data(), 8);
  if (!is) {
      std::cerr << "Error: Could not read the file (StartTime) " << path << std::endl;
      throw std::runtime_error("Failed to read file");
  }
  StartTime = std::string(buffer.data(), 8);

  // Read SizeHeader
  is.read(buffer.data(), 8);
  if (!is) {
      std::cerr << "Error: Could not read the file (Number of bytes) " << path << std::endl;
      throw std::runtime_error("Failed to read file");
  }
  SizeHeader = std::string(buffer.data(), 8);

 // Read Reserved
  is.read(buffer.data(), 44);
  if (!is) {
      std::cerr << "Error: Could not read the file (Reserved) " << path << std::endl;
      throw std::runtime_error("Failed to read file");
  }
  Reserved = std::string(buffer.data(), 44);

 // Read NumberDataRecords
  is.read(buffer.data(), 8);
  if (!is) {
      std::cerr << "Error: Could not read the file (NumberDataRecords) " << path << std::endl;
      throw std::runtime_error("Failed to read file");
  }
  NumDataRecords = std::stoi(std::string(buffer.data(), 8));

 // Read DurationDataRecords
  is.read(buffer.data(), 8);
  if (!is) {
      std::cerr << "Error: Could not read the file (DurationDataRecords) " << path << std::endl;
      throw std::runtime_error("Failed to read file");
  }
  DurationDataRecords = std::string(buffer.data(), 8);

  // Read DurationDataRecords
  is.read(buffer.data(), 4);
  if (!is) {
      std::cerr << "Error: Could not read the file (NumberSignals) " << path << std::endl;
      throw std::runtime_error("Failed to read file");
  }
  NumberSignals = std::stoi(std::string(buffer.data(), 4));
  //-------------------------------------------------------------------------------------------------------//
  // parsing  data records
  SignalsInfo.resize(NumberSignals);
  std::cout<< "Number of signals: "<< NumberSignals<<std::endl;
  int ns = NumberSignals;
  //read label(s)
  for( int i = 0 ; i< ns ; i++ ){

    is.read(buffer.data(), 16);
    if (!is) {
        std::cerr << "Error: Could not read the file (label) in signal number: " << i << std::endl;
        throw std::runtime_error("Failed to read file");
    }
    SignalsInfo[i].Label =  std::string(buffer.data(), 16);

  }

  //read transducer
  for( int i = 0 ; i<ns ; i++ ){

    is.read(buffer.data(), 80);
    if (!is) {
        std::cerr << "Error: Could not read the file (transducer) in signal number: " << i << std::endl;
        throw std::runtime_error("Failed to read file");
    }
    SignalsInfo[i].Transducer =  std::string(buffer.data(), 80);

  }

  //read physical dimension
  for( int i = 0 ; i<ns ; i++ ){

    is.read(buffer.data(), 8);
    if (!is) {
        std::cerr << "Error: Could not read the file (physical dimension) in signal number: " << i << std::endl;
        throw std::runtime_error("Failed to read file");
    }
    SignalsInfo[i].Label =  std::string(buffer.data(), 8);

  }

  //read physical minimum
  for( int i = 0 ; i<ns ; i++ ){

    is.read(buffer.data(), 8);
    if (!is) {
        std::cerr << "Error: Could not read the file (physical minimum) in signal number: " << i << std::endl;
        throw std::runtime_error("Failed to read file");
    }
    SignalsInfo[i].PhysicalMinimum =  std::string(buffer.data(), 8);

  }

  //read physical maximum
  for( int i = 0 ; i<ns ; i++ ){

    is.read(buffer.data(), 8);
    if (!is) {
        std::cerr << "Error: Could not read the file (physical maximum) in signal number: " << i << std::endl;
        throw std::runtime_error("Failed to read file");
    }
    SignalsInfo[i].PhysicalMaximum =  std::string(buffer.data(), 8);

  }

  //read digital minimum
  for( int i = 0 ; i<ns ; i++ ){

    is.read(buffer.data(), 8);
    if (!is) {
        std::cerr << "Error: Could not read the file (digital minimum) in signal number: " << i << std::endl;
        throw std::runtime_error("Failed to read file");
    }
    SignalsInfo[i].DigitalMinimum =  std::string(buffer.data(), 8);

  }

  //read digital maximum
  for( int i = 0 ; i<ns ; i++ ){

    is.read(buffer.data(), 8);
    if (!is) {
        std::cerr << "Error: Could not read the file (digital maximum) in signal number: " << i << std::endl;
        throw std::runtime_error("Failed to read file");
    }
    SignalsInfo[i].DigitalMaximum =  std::string(buffer.data(), 8);

  }

  //read prefiltering
  for( int i = 0 ; i<ns ; i++ ){

    is.read(buffer.data(), 80);
    if (!is) {
        std::cerr << "Error: Could not read the file (prefiltering) in signal number: " << i << std::endl;
        throw std::runtime_error("Failed to read file");
    }
    SignalsInfo[i].Prefiltering =  std::string(buffer.data(), 80);

  }

  //read nr
  for( int i = 0 ; i<ns ; i++ ){

    is.read(buffer.data(), 8);
    if (!is) {
        std::cerr << "Error: Could not read the file (nr) in signal number: " << i << std::endl;
        throw std::runtime_error("Failed to read file");
    }
    SignalsInfo[i].Nr = std::stoi(std::string(buffer.data(), 8)) ;

  }

  //read reserved
  for( int i = 0 ; i<ns ; i++ ){

    is.read(buffer.data(), 32);
    if (!is) {
        std::cerr << "Error: Could not read the file (reserved) in signal number: " << i << std::endl;
        throw std::runtime_error("Failed to read file");
    }
    SignalsInfo[i].ReservedDos =  std::string(buffer.data(), 32);

  }

  //reading Signals


  /*for(int i= 0 ;i < ns; i++){

    std::vector<int16_t> signalData;
    for(int j = 0 ; j < SignalsInfo[i].Nr; j++){
      
      is.read(buffer.data(), 2);
      if (!is) {
          std::cerr << "Error: Could not read the file (Reading Signal) "<< std::endl;
          throw std::runtime_error("Failed to read file");
      }

      int16_t value = (static_cast<uint8_t>(buffer[1]) << 8) | static_cast<uint8_t>(buffer[0]);
      std::cout << "value casted: "<< value <<std::endl;
      signalData.push_back(value);


    }

    Signals.push_back(signalData);

  }*/
  std::vector< std::vector<std::complex<double> > > vectors(NumberSignals);

  for(int i=0 ;i<NumDataRecords;i++){
    for (int j = 0 ; j < NumberSignals;j++){

      for(int k = 0 ; k< SignalsInfo[j].Nr ;k++){
        std::int16_t value;
        is.read(reinterpret_cast<char*>(&value), sizeof(int16_t));
        if (!is) {
          std::cerr << "Error: Could not read the file (Reading Signal) "<< std::endl;
          throw std::runtime_error("Failed to read file");
        }
        
        //std::cout << "value casted: "<< value <<std::endl;
        double realPart = static_cast<double>(value); // Convert to double
    
        std::complex<double> complexNum(realPart, 0.0); // Create complex number with real part and imaginary part as 0

        vectors[j].push_back(value);
      }

    }

  }
  for(auto const newSignal : vectors){

    Eigen::VectorXcd eigenVec(newSignal.size());
    for (size_t i = 0; i < newSignal.size(); ++i) {
        eigenVec[i] = newSignal[i];
    }
    Signals.Signals.push_back(eigenVec);
  } 



  /*std::cout<< "APARTIR DE AQUI ES FALTANTE" <<std::endl;
  int16_t value;
  while (!is.eof()) {
    is.read(buffer.data(),1);
    if (!is) {
      std::cerr << "Error: Could not read the file (Reading Signal) "<< std::endl;
      throw std::runtime_error("Failed to read file");
    }
    std::cout << "falto: "<< std::string(buffer.data(), 1)<<std::endl;
  }*/

}



void edf::PrintHeaderRecords(){
  std::cout <<"Version: " << Version<<std::endl;
  std::cout<<"Local patient identification: " << PatientId<<std::endl;
  std::cout<<"Local recording identification: "<< RecordingId<<std::endl;
  std::cout<<"Startdate of recording: "<< StartDate<<std::endl;
  std::cout<<"Start time of recording: " << StartTime<<std::endl;
  std::cout<<"Number of bytes in header record: " << SizeHeader <<std::endl;
  std::cout<<"Reserved: " << Reserved <<std::endl;
  std::cout<<"Number of data records: " << NumDataRecords<<std::endl;
  std::cout<<"Duration of a data: " <<DurationDataRecords <<std::endl;
  std::cout<<"Number of signals: " << NumberSignals<<std::endl;

}

void edf::PrintDataRecords(){

  for(const auto Data : SignalsInfo){
    std::cout << "------------------------------------------------------------------"<<std::endl;
    std::cout <<"Label: "<< Data.Label<<std::endl;
    std::cout <<"Transducer: "<< Data.Transducer<<std::endl;
    std::cout <<"Physical dimension: " << Data.PhyisicalDimension<<std::endl;
    std::cout <<"Physical minimum: " << Data.PhysicalMinimum<<std::endl;
    std::cout <<"Physical maximum: " << Data.PhysicalMaximum<<std::endl;
    std::cout <<"Digital minimum: " << Data.DigitalMinimum<<std::endl;
    std::cout <<"Digital maximum: " << Data.DigitalMaximum<<std::endl;
    std::cout <<"Prefiltering:" << Data.Prefiltering<<std::endl;
    std::cout <<"Nr: " << Data.Nr<<std::endl;
    std::cout <<"Reserved: " << Data.ReservedDos<<std::endl;
    std::cout << "------------------------------------------------------------------"<<std::endl;

  }
  

}

void edf::PrintSizeSignals(){

  for(const auto signali : Signals.Signals){
    
    std::cout<<"Size of the vector is: "<< signali.size() <<std::endl;

  }


}

void edf::WriteRawCsv(const std::string filename){

  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open the file " << filename << std::endl;
    throw std::runtime_error("Failed to open file");
  }

  
    // Write the header
    for (int i = 0; i < NumberSignals; ++i) {
        file << "Signal " << i;
        if (i < NumberSignals - 1) {
            file << ",";
        }
    }
    file << "\n";

    // Find the maximum length of signals
    size_t max_length = 0;
    for (const auto& signal : Signals.Signals) {
        if (signal.size() > max_length) {
            max_length = signal.size();
        }
    }

    // Write the signals
    for (size_t i = 0; i < max_length; ++i) {
        for (int j = 0; j < NumberSignals; ++j) {
            if (i < Signals.Signals[j].size()) {
                file << Signals.Signals[j][i];
            }
            if (j < NumberSignals - 1) {
                file << ",";
            }
        }
        file << "\n";
    }




}


void edf::PrintTopValues(int i){


  for(const auto signali : Signals.Signals){
    
    for(int j=0;j<i;j++){
      std::cout<< signali[j] <<" ";
      
    }

    std::cout<< std::endl;

  }

}

