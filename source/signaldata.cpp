#include "signaldata.hpp"
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <ostream>

const double PI = acos(-1);

SignalData::SignalData() : Signals{}, Means{}, StdDeviation{} {}

SignalData::SignalData(std::vector<Eigen::VectorXcd> InputSignal)
    : Signals{InputSignal} {}

void SignalData::CalculateMeans() {
    Means.clear();
    for (const auto& signal : Signals) {
        if (signal.size() > 0) {
            std::complex<double> sum = std::accumulate(signal.begin(), signal.end(), std::complex<double>(0.0, 0.0));
            Means.push_back(sum / static_cast<double>(signal.size()));
        } else {
            Means.push_back(0.0);
        }
    }
}

void SignalData::CalculateDeviation() {
    StdDeviation.clear();
    CalculateMeans(); // Asegurarse de que los medios estén actualizados
    for (size_t i = 0; i < Signals.size(); ++i) {
        const auto& signal = Signals[i];
        auto mean = Means[i];
        if (signal.size() > 0) {
            double variance = std::accumulate(signal.begin(), signal.end(), 0.0, 
                                              [mean](double acc, const std::complex<double>& value) {
                                                  auto diff = value - mean;
                                                  return acc + std::norm(diff);
                                              }) / static_cast<double>(signal.size());
            StdDeviation.push_back(std::sqrt(variance));
        } else {
            StdDeviation.push_back(0.0);
        }
    }
}

void SignalData::PrintMeanAndDeviation(){

  if(!(Signals.size() == Means.size() == StdDeviation.size())){
    std::cout<< "Mean or standar deviation not propertly initialized" <<std::endl;
  }
  int Size = Signals.size();
  for( int i = 0 ; i < Size ; i++ ){
    std::cout<< "signal number "<< i << " mean = "<< Means[i] << " standart deviation = " << StdDeviation[i] << std::endl;
  }


}


void SignalData::GenerateRandomSignals(size_t numSignals, size_t numSamples, std::complex<double> mean, double stddev) {
    Signals.clear();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist_real(mean.real(), stddev);
    std::normal_distribution<> dist_imag(mean.imag(), stddev);

    for (size_t i = 0; i < numSignals; ++i) {
        Eigen::VectorXcd signal(numSamples);
        for (size_t j = 0; j < numSamples; ++j) {
            signal[j] = std::complex<double>(dist_real(gen), dist_imag(gen));
        }
        Signals.push_back(signal);
    }
}

Eigen::VectorXcd SignalData::FirstDifference(const Eigen::VectorXcd& Input){
  size_t Size = Input.size();

  Eigen::VectorXcd Output(Size);

  if(Input.size()<0) return Output;

  Output[0] = Input[0];

  for(int i = 1;i < Size ; i++ ){
    Output[i] = Input[i] - Input[i-1];
  }

  return Output;
}

Eigen::VectorXcd SignalData::RuningSum(const Eigen::VectorXcd& Input){
  int Size = Input.size();

  Eigen::VectorXcd Output(Size);

  if(Input.size()< 0) return Output;
  Output[0] = Input[0];
  for (size_t i = 1; i < Size; ++i) {
      Output[i] = Output[i - 1] + Input[i];
  }
  return Output;
}


Eigen::VectorXcd SignalData::FFTAux(Eigen::VectorXcd a, bool invert){



  if ((a.size() & (a.size() - 1)) != 0) {
      a = ZeroPadPowerTwo(a);
  }


  int n = a.size();
  if (n<=1) return a;
    // Crear los vectores even y odd
    // Ensure even and odd are correctly sized
  Eigen::VectorXcd even = a(Eigen::seq(0, n - 1, 2));
  Eigen::VectorXcd odd = a(Eigen::seq(1, n - 1, 2));
    // Recursive FFT calls
  even = FFTAux(even, invert);
  odd = FFTAux(odd, invert);

  double angle = 2 * PI / n * (invert ? 1 : -1);
  std::complex<double> w(1) ,wn(cos(angle),sin(angle));
  Eigen::VectorXcd result(n);

    for (int k = 0; k < n / 2; ++k) {
        std::complex<double> t = w * odd[k];
        result[k] = even[k] + t;
        result[k + n / 2] = even[k] - t;

        if (invert) {
            result[k] /= 2;
            result[k + n / 2] /= 2;
        }

        w *= wn;
    }

    return result;




}

Eigen::VectorXcd SignalData::FFT(Eigen::VectorXcd &a){
    int OrgSize = a.size();
    return FFTAux(a, false).head(OrgSize);
}


Eigen::VectorXcd SignalData::IFFT(Eigen::VectorXcd &a){
    int OrgSize = a.size();
    return FFTAux(a, true).head(OrgSize);
}


Eigen::VectorXcd SignalData::ZeroPadPowerTwo(const Eigen::VectorXcd &a) {
    int n = a.size();
    int m = 1;

    // Find the next power of 2 greater than or equal to n
    while (m < n) m <<= 1;
    // If n is already a power of two, return the original array
    if (m == n) {
        return a;
    }
    // Directly initializing the new VectorXcd with zero-padding
    Eigen::VectorXcd padded_a = Eigen::VectorXcd::Zero(m);

    // Copy the input vector into the beginning of the padded vector
    padded_a.head(n) = a;

    return padded_a;
}


Eigen::VectorXcd SignalData::MovingAverage(const Eigen::VectorXcd &a, int window_size) {
    int n = a.size();
    Eigen::VectorXcd result(n);

    // Asegurarse de que el tamaño de la ventana sea válido
    if (window_size <= 0 || window_size > n) {
        throw std::invalid_argument("El tamaño de la ventana debe ser mayor que 0 y menor o igual al tamaño de la señal.");
    }

    // Inicializar la suma para la primera ventana
    std::complex<double> sum = 0.0;
    for (int i = 0; i < window_size; ++i) {
        sum += a[i];
    }

    // Calcular el promedio para el primer elemento del resultado
    result[0] = sum / static_cast<double>(window_size);

    // Aplicar el filtro
    for (int i = 1; i <= n - window_size; ++i) {
        sum += a[i + window_size - 1];  // Añadir el nuevo elemento en la ventana
        sum -= a[i - 1];  // Restar el elemento que sale de la ventana
        result[i] = sum / static_cast<double>(window_size);
    }

    // Rellenar los últimos elementos (si la ventana no cubre toda la señal)
    for (int i = n - window_size + 1; i < n; ++i) {
        result[i] = result[i - 1];
    }

    return result;
}
std::pair<Eigen::VectorXcd, Eigen::VectorXcd> SignalData::FastWaveletTransformHaar(const Eigen::VectorXcd& input) {
    size_t n = input.size();
    Eigen::VectorXcd approximation(n / 2);
    Eigen::VectorXcd detail(n / 2);

    // Haar Wavelet Filters
    std::vector<double> low_pass_filter = {0.5, 0.5};
    std::vector<double> high_pass_filter = {0.5, -0.5};

    for (size_t i = 0; i < n / 2; ++i) {
        approximation[i] = low_pass_filter[0] * input[2 * i] + low_pass_filter[1] * input[2 * i + 1];
        detail[i] = high_pass_filter[0] * input[2 * i] + high_pass_filter[1] * input[2 * i + 1];
    }

    return std::make_pair(approximation, detail);
}
