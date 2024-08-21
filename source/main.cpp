#include "edf.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
const double PI = acos(-1);



void ExportToCSV(const Eigen::VectorXcd& approximation, const Eigen::VectorXcd& detail, const std::string& filename) {
    std::ofstream file(filename);
    
    file << "Index,Approximation,Detail\n";
    for (size_t i = 0; i < approximation.size(); ++i) {
        file << i << "," << approximation[i].real() << "," << detail[i].real() << "\n";
    }
    
    file.close();
}
void ExportToCSVDecomp(const std::vector<Eigen::VectorXcd>& approximations, 
                 const std::vector<Eigen::VectorXcd>& details, 
                 const std::string& filename) {
    std::ofstream file(filename);

    file << "Level,Index,Approximation,Detail\n";
    for (size_t level = 0; level < approximations.size(); ++level) {
        const Eigen::VectorXcd& approximation = approximations[level];
        const Eigen::VectorXcd& detail = details[level];
        for (size_t i = 0; i < approximation.size(); ++i) {
            file << level + 1 << "," << i << "," << approximation[i].real() << "," << detail[i].real() << "\n";
        }
    }

    file.close();
}

void write_to_csv(const std::string& filename, const Eigen::VectorXcd& data) {
    std::ofstream file(filename);

    if (file.is_open()) {
        file << "Index,Real,Imaginary\n";
        for (int i = 0; i < data.size(); ++i) {
            file << i << "," << data[i].real() << "," << data[i].imag() << "\n";
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
Eigen::VectorXcd GenerateChirpSignal(size_t numSamples, double startFreq, double endFreq, double duration, double samplingRate) {
    Eigen::VectorXcd chirpSignal(numSamples);
    double t;

    for (size_t i = 0; i < numSamples; ++i) {
        t = i / samplingRate;  // Current time
        double freq = startFreq + ((endFreq - startFreq) / duration) * t;  // Instantaneous frequency
        chirpSignal[i] = std::polar(1.0, 2 * PI * freq * t);  // Generate the chirp signal
    }

    return chirpSignal;
}
Eigen::VectorXcd CreateGaborFilter(int size, double center_freq, double sigma) {
    Eigen::VectorXcd filter(size);
    double T = 1.0 / size; // Sampling period
    double mid = size / 2.0;

    for (int i = 0; i < size; ++i) {
        double t = (i - mid) * T; // Time index centered at zero
        // Gabor filter formula: Gaussian envelope * cosine wave
        filter[i] = std::exp(-t * t / (2 * sigma * sigma)) * std::cos(2 * PI * center_freq * t);
    }

    return filter;
}
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
        std::cout << std::endl;
        
        int N = 32768; // Number of samples
        double T = 1.0 / N; // Sampling period
        double freq = 2.0; // Frequency of the sine wave

        Eigen::VectorXcd sine_wave(N);

        // Generate sine wave
        for (int i = 0; i < N; ++i) {
            double t = i * T; // Time index
            sine_wave[i] = std::sin(2 * PI * freq * t); // Sine wave formula
        }

        // Write original sine wave to CSV
        write_to_csv("sine_wave_before_fft.csv", sine_wave);

        auto after_fft = edfFile.Signals.FFT(sine_wave);

        // Write FFT result to CSV
        write_to_csv("sine_wave_after_fft.csv", after_fft);

        auto after_iFFT = edfFile.Signals.IFFT(after_fft);

        write_to_csv("sine_wave_after_Ifft.csv", after_iFFT);

        std::cout << "Sine wave data has been written to CSV files.\n";

        write_to_csv("test_wave_before_fft.csv", edfFile.Signals.Signals[0]);

        auto after_fft2 = edfFile.Signals.FFT(edfFile.Signals.Signals[0]);

        // Write FFT result to CSV
        write_to_csv("test_wave_after_fft.csv", after_fft2);

        auto after_iFFT2 = edfFile.Signals.IFFT(after_fft2);

        write_to_csv("test_wave_after_Ifft.csv", after_iFFT2);

        auto aftter_moving_avergae = edfFile.Signals.MovingAverage(edfFile.Signals.Signals[0], 11);

        write_to_csv("test_wave_after_Avergae.csv", aftter_moving_avergae);

        std::cout << "Sine wave data has been written to CSV files.\n";
      // Create Gabor filter
      int filter_size = 1024; // Choose a smaller size for the filter
      double center_freq = 5.0; // Center frequency of the Gabor filter
      double sigma = 0.005; // Width of the Gaussian envelope
      Eigen::VectorXcd gabor_filter = CreateGaborFilter(filter_size, center_freq, sigma);

      // Perform convolution using the filter
      SignalData signalData; // Assuming SignalData is your class
      auto start = std::chrono::high_resolution_clock::now();
      Eigen::VectorXcd convolved_signal = signalData.Convolve(sine_wave, gabor_filter);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> direct_conv_time = end - start;
      std::cout << "Direct Convolution Time: " << direct_conv_time.count() << " seconds" << std::endl;
      write_to_csv("test_wave_after_convolve.csv", convolved_signal);

      // Alternatively, you can use FFT-based convolution
      start = std::chrono::high_resolution_clock::now();
      Eigen::VectorXcd fft_convolved_signal = signalData.FFTConvolve(sine_wave, gabor_filter);
      end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> fft_conv_time = end - start;
      std::cout << "FFT Convolution Time: " << fft_conv_time.count() << " seconds" << std::endl;

      write_to_csv("test_wave_after_fftConvolve.csv", fft_convolved_signal);*/

      size_t numSamples = 1024;  // Number of samples
      double startFreq = 20.0;   // Start frequency in Hz
      double endFreq = 100.0;    // End frequency in Hz
      double duration = 1.0;     // Duration of the signal in seconds
      double samplingRate = 1024.0; // Sampling rate in Hz

      // Generate the chirp signal

      Eigen::VectorXcd chirpSignal = GenerateChirpSignal(numSamples, startFreq, endFreq, duration, samplingRate);
      write_to_csv("chirp_wave_before_.csv", chirpSignal);
      auto [approximation, detail] =edfFile.Signals.FastWaveletTransformHaar(chirpSignal);
      ExportToCSV(approximation, detail, "wavelet_output.csv");
      
      auto [approximations, details] = edfFile.Signals.FastWaveletTransform(chirpSignal, 3, "bior3.1");
      ExportToCSVDecomp(approximations, details, "wavelet_output_decomp.csv");


    } catch (const std::runtime_error& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }
  return 0;
}
