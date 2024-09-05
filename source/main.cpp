#include "edf.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
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
// Function to write CWT results to CSV
void ExportCWTToCSV(const std::vector<Eigen::VectorXcd>& cwtCoeffs, const std::vector<double>& scales, const std::string& filename) {
    std::ofstream file(filename);

    // Write header
    file << "Scale,Index,Real,Imaginary\n";

    // Write data
    for (size_t i = 0; i < scales.size(); ++i) {
        const auto& coeffs = cwtCoeffs[i];
        for (int j = 0; j < coeffs.size(); ++j) {
            file << scales[i] << "," << j << "," << coeffs[j].real() << "," << coeffs[j].imag() << "\n";
        }
    }

    file.close();
}

std::vector<double> GenerateLogScales(double start, double end, int numScales) {
    std::vector<double> scales(numScales);
    double logStart = std::log10(start);
    double logEnd = std::log10(end);
    double logStep = (logEnd - logStart) / (numScales - 1);

    for (int i = 0; i < numScales; ++i) {
        scales[i] = std::pow(10, logStart + i * logStep);
    }

    return scales;
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
        /*auto moving = edfFile.Signals.MovingAverage(edfFile.Signals.Signals[0], 10);
        std::cout<< moving(Eigen::seq(0,10))<<std::endl;
        auto sum = edfFile.Signals.RuningSum(edfFile.Signals.Signals[0]);
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
      edfFile.Signals.FastWaveletTransform(edfFile.Signals.Signals[0], 3, "bior3.1");
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

      auto morlet = edfFile.Signals.MorletWavelet(chirpSignal.size(), 50);
      write_to_csv("Morlet_Wavelet.csv", morlet);
      
      auto morlet_after_fft = edfFile.Signals.FFT(morlet);
      
      write_to_csv("Morlet_Wavelet_afterFFT.csv", morlet_after_fft);
      auto morlet_after_fft_Eigen = edfFile.Signals.FFTEigen(morlet);
      write_to_csv("Morlet_Wavelet_afterFFT_Eigem.csv", morlet_after_fft_Eigen);
      auto morlet_after_Ifft = edfFile.Signals.IFFT(morlet_after_fft);
        
      write_to_csv("Morlet_Wavelet_afterIFFT.csv", morlet_after_Ifft);

      double start = 1.0;
      double end = 1024.0;
      int numScales = 100;

      std::vector<double> scales = GenerateLogScales(start, end, numScales);
      auto start2 = std::chrono::high_resolution_clock::now();

      auto coeffs = edfFile.Signals.CWT(chirpSignal, scales);
      auto end2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> direct_conv_time2 = end2 - start2;
      std::cout << "CWT time: " << direct_conv_time2.count() << " seconds" << std::endl;

      ExportCWTToCSV(coeffs,scales, "CWT.csv");

      start2 = std::chrono::high_resolution_clock::now();

      coeffs = edfFile.Signals.CWTEigen(chirpSignal, scales);
      end2 = std::chrono::high_resolution_clock::now();
      direct_conv_time2 = end2 - start2;
      std::cout << "CWT time: " << direct_conv_time2.count() << " seconds" << std::endl;

      ExportCWTToCSV(coeffs,scales, "CWTEigen.csv");     

      numSamples = 2000;  // Number of samples
      startFreq = 0.2;   // Start frequency in Hz
      endFreq = 9;    // End frequency in Hz
      duration = 2.0;     // Duration of the signal in seconds
      samplingRate = 1024.0; // Sampling rate in Hz

      // Generate the chirp signal

      Eigen::VectorXcd chirp1 = GenerateChirpSignal(numSamples, startFreq, endFreq, duration, samplingRate);

       numSamples = 2000;  // Number of samples
      startFreq = 0.1;   // Start frequency in Hz
      endFreq = 5;    // End frequency in Hz
      duration = 2.0;     // Duration of the signal in seconds
      samplingRate = 1024.0; // Sampling rate in Hz

      // Generate the chirp signal

      Eigen::VectorXcd chirp2 = GenerateChirpSignal(numSamples, startFreq, endFreq, duration, samplingRate);

       Eigen::VectorXcd chirp = chirp1 + 0.6 * chirp2;
       write_to_csv("chirp_strange.csv", chirp);

      coeffs = edfFile.Signals.CWTEigen(chirp, scales);

      ExportCWTToCSV(coeffs,scales, "CWT2.csv");
      start2 = std::chrono::high_resolution_clock::now();

      coeffs = edfFile.Signals.CWTEigen(chirp, scales);
      end2 = std::chrono::high_resolution_clock::now();
      direct_conv_time2 = end2 - start2;
      std::cout << "CWT2Eigen time: " << direct_conv_time2.count() << " seconds" << std::endl;

      ExportCWTToCSV(coeffs,scales, "CWTEigen2.csv");     

    } catch (const std::runtime_error& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }
  return 0;
}
