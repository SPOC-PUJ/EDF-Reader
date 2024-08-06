#include "edf.hpp"
#include <iostream>
#include <ostream>
#include <string>
int main (int argc, char *argv[]) {
  
    if(argc != 2){
        std::cerr << "Usage: " << argv[0] << " <file_path>" << std::endl;
        return 1;
    }

    std::cout << "Opening file ... " << argv[1] << std::endl;

    try {
        edf edfFile { std::string(argv[1]) };
        std::cout << "File opened and read successfully." << std::endl;
        edfFile.PrintHeaderRecords();
        edfFile.PrintDataRecords();
        edfFile.PrintSizeSignals();
        //edfFile.WriteRawCsv("signal.csv");
        edfFile.PrintTopValues(10);
        edfFile.Signals.CalculateDeviation();
        edfFile.Signals.PrintMeanAndDeviation();
        /*auto sum = edfFile.Signals.RuningSum(edfFile.Signals.Signals[0]);
        std::cout << "Sum size = " << sum.size() <<std::endl;
        for(const auto num : sum ){
          std::cout << num << " "; 
        }
        std::cout << std::endl;
        auto Diff = edfFile.Signals.FirstDifference(edfFile.Signals.Signals[0]);
        std::cout << "Diff size = " << Diff.size() <<std::endl;
        for(const auto di : Diff ){
          std::cout << di << " "; 
        }
        std::cout << std::endl;*/



    } catch (const std::runtime_error& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }
  return 0;
}
