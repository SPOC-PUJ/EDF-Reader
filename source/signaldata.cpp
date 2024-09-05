#include "signaldata.hpp"
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include "../unsupported/Eigen/FFT"

const double PI = acos(-1);
    // db1 :
    const Eigen::VectorXcd Db1DecLow {{std::complex<double>(0.7071067812, 0.0), std::complex<double>(0.7071067812, 0.0)}};
    
    const Eigen::VectorXcd Db1DecHigh   {{std::complex<double>(-0.7071067812, 0.0), std::complex<double>(0.7071067812, 0.0)}};

    // db2 :
    const Eigen::VectorXcd Db2DecLow  {{std::complex<double>(-0.1294095226, 0.0), std::complex<double>(0.2241438680, 0.0),std::complex<double>(0.8365163037, 0.0), std::complex<double>(0.4829629131, 0.0)}};

    const Eigen::VectorXcd Db2DecHigh {{std::complex<double>(-0.4829629131, 0.0), std::complex<double>(0.8365163037, 0.0),std::complex<double>(-0.2241438680, 0.0), std::complex<double>(-0.1294095226, 0.0)}};


    // db3 :
    const Eigen::VectorXcd Db3DecLow  {{ std::complex<double>(0.0352262919, 0.0), std::complex<double>(-0.0854412739, 0.0),std::complex<double>(-0.1350110200, 0.0), std::complex<double>(0.4598775021, 0.0),std::complex<double>(0.8068915093, 0.0), std::complex<double>(0.3326705530, 0.0)}};

    const Eigen::VectorXcd Db3DecHigh  {{ std::complex<double>(-0.3326705530, 0.0), std::complex<double>(0.8068915093, 0.0),std::complex<double>(-0.4598775021, 0.0), std::complex<double>(-0.1350110200, 0.0),std::complex<double>(0.0854412739, 0.0), std::complex<double>(0.0352262919, 0.0)}};

    // db4 :
    const Eigen::VectorXcd Db4DecLow  {{std::complex<double>(-0.0105974018, 0.0), std::complex<double>(0.0328830117, 0.0),std::complex<double>(0.0308413818, 0.0), std::complex<double>(-0.1870348117, 0.0),std::complex<double>(-0.0279837694, 0.0), std::complex<double>(0.6308807679, 0.0),std::complex<double>(0.7148465706, 0.0), std::complex<double>(0.2303778133, 0.0)}};

    const Eigen::VectorXcd Db4DecHigh  {{std::complex<double>(-0.2303778133, 0.0), std::complex<double>(0.7148465706, 0.0),std::complex<double>(-0.6308807679, 0.0), std::complex<double>(-0.0279837694, 0.0),std::complex<double>(0.1870348117, 0.0), std::complex<double>(0.0308413818, 0.0),std::complex<double>(-0.0328830117, 0.0), std::complex<double>(-0.0105974018, 0.0)}};

    // bior 3.1
    const Eigen::VectorXcd Bio31DecLow  {{std::complex<double>(-0.3535533905932738, 0.0), std::complex<double>(1.0606601717798212, 0.0),std::complex<double>(1.0606601717798212, 0.0), std::complex<double>(-0.3535533905932738, 0.0)}};

    const Eigen::VectorXcd Bio31DecHigh  {{std::complex<double>(-0.1767766952966369, 0.0), std::complex<double>(0.5303300858899106, 0.0),std::complex<double>(-0.5303300858899106, 0.0), std::complex<double>(0.1767766952966369, 0.0)}};


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

Eigen::VectorXcd SignalData::ZeroPadGivenSize(const Eigen::VectorXcd &a,int m) {
    int n = a.size();

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
// Assuming this method is part of SignalData class
Eigen::VectorXcd SignalData::Convolve(const Eigen::VectorXcd& x, const Eigen::VectorXcd& h) {
    int x_size = x.size();
    int h_size = h.size();
    int y_size = x_size + h_size - 1;

    Eigen::VectorXcd y = Eigen::VectorXcd::Zero(y_size);

    for (int n = 0; n < y_size; ++n) {
        for (int m = std::max(0, n + 1 - h_size); m <= std::min(n, x_size - 1); ++m) {
            y[n] += x[m] * h[n - m];
        }
    }

    return y;
}


Eigen::VectorXcd SignalData::FFTConvolve(const Eigen::VectorXcd& x, const Eigen::VectorXcd& h) {
    // Calculate the size for zero padding (next power of two greater than x_size + h_size - 1)
    int n = x.size() + h.size() - 1;
    int m = 1;
    while (m < n) m <<= 1;

    // Zero pad both signals
    Eigen::VectorXcd x_padded = ZeroPadGivenSize(x,m);
    Eigen::VectorXcd h_padded = ZeroPadGivenSize(h,m);

    // Perform FFT on both signals
    Eigen::VectorXcd X = FFT(x_padded);
    Eigen::VectorXcd H = FFT(h_padded);

    // Multiply the FFT results element-wise
    Eigen::VectorXcd Y = X.cwiseProduct(H);

    // Perform IFFT to obtain the convolved signal
    Eigen::VectorXcd y_padded = IFFT(Y);

    // Trim the padded result to the expected convolution length
    return y_padded.head(n);
}

Eigen::VectorXcd SignalData::FFTEigen(Eigen::VectorXcd &a){
  
  Eigen::FFT<double> ft; 
  Eigen::VectorXcd signalFFT;
  ft.fwd(signalFFT,a);
  return signalFFT;
}

Eigen::VectorXcd SignalData::IFFTEigen(Eigen::VectorXcd &a){
  Eigen::FFT<double> ft; 
  Eigen::VectorXcd signalFFT;
  ft.inv(signalFFT,a);
  return signalFFT;
}


Eigen::VectorXcd SignalData::FFTconvolveEigen(const Eigen::VectorXcd& x, const Eigen::VectorXcd& h, bool shift){
  int N = x.size();

  auto h_padded = ZeroPadGivenSize(h, N);
  Eigen::VectorXcd signalFFT;
  Eigen::VectorXcd Kernel;

  Eigen::FFT<double> ft;

  ft.fwd(signalFFT, x);
  ft.fwd(Kernel, h_padded);
  if(shift){
    // Frequency domain shift to center the wavelet
    Eigen::VectorXcd phaseShift(N);
    for (int k = 0; k < N; ++k) {
        double freq = 2.0 * M_PI * k / N;
        phaseShift[k] = std::exp(std::complex<double>(0, -freq * (static_cast<double>(N) / 2)));
    }
    Kernel = Kernel.cwiseProduct(phaseShift);

  }
  Eigen::VectorXcd convFFT = signalFFT.cwiseProduct(Kernel);
  Eigen::VectorXcd convolutionResult;
  ft.inv(convolutionResult, convFFT);
  return convolutionResult.head(N);
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

std::pair< std::vector< Eigen::VectorXcd >, std::vector< Eigen::VectorXcd > > SignalData::FastWaveletTransform(const Eigen::VectorXcd& input,int DecLevel,std::string WaveName){

    std::vector<Eigen::VectorXcd>Aproximations;
    std::vector<Eigen::VectorXcd>Details;
    const auto getWaveletFilters = [](const std::string& waveName) -> std::pair<Eigen::VectorXcd, Eigen::VectorXcd> {
        if (waveName == "db1") {
            return std::make_pair(Db1DecLow, Db1DecHigh);
        } else if (waveName == "db2") {
            return std::make_pair(Db2DecLow, Db2DecHigh);
        } else if (waveName == "db3") {
            return std::make_pair(Db3DecLow, Db3DecHigh);
        } else if (waveName == "db4") {
            return std::make_pair(Db4DecLow, Db4DecHigh);
        } else if (waveName == "bior3.1") {
            return std::make_pair(Bio31DecLow, Bio31DecHigh);
        } else {
            throw std::invalid_argument("Unknown wavelet name");
        }
    };

    auto [Lo_d, Hi_d] = getWaveletFilters(WaveName);


    auto [aprox ,detail ] = FastWaveletTransformAux(input,Lo_d,Hi_d );
    Aproximations.push_back(aprox);
    Details.push_back(detail);
    if(DecLevel > 1){

      for( int i = 1 ; i < DecLevel ; i++ ){

        auto descompocision = FastWaveletTransformAux(aprox,Lo_d,Hi_d );
        aprox = descompocision.first;
        detail = descompocision.second;
        Aproximations.push_back(aprox);
        Details.push_back(detail);
      }


    }
  return std::make_pair(Aproximations, Details);
}

std::pair<Eigen::VectorXcd,Eigen::VectorXcd> SignalData::FastWaveletTransformAux(const Eigen::VectorXcd& input,Eigen::VectorXcd Lo_d , Eigen::VectorXcd Hi_d){

  // conseguir filtros segun el nombre
  // Tipos de filtros

  // conseguir aproximation FFTConvolve con low pass 
    auto aprox = FFTconvolveEigen(input, Lo_d );
  // conseguir detail FFTConvolve con High pass
    auto detail = FFTconvolveEigen(input, Hi_d );
  // hacer down sample de aproximation y de detail

 
    int n = aprox.size();
    int m = detail.size();

    auto aproxi = aprox(Eigen::seq(0,n-1,2));
    auto detaili = detail(Eigen::seq(0,m-1,2));

  // llamar la descompocision denuevo hata que se acben los niveles

  return std::make_pair(aproxi, detaili);
}



Eigen::VectorXcd SignalData::MorletWavelet(int N, double scale, double f0, double fb) {
  Eigen::VectorXcd wavelet(N);
    double t;
    double normalization = 1 / sqrt(M_PI*fb);
    for (int n = 0; n < N; ++n) {
        t = (n - N / 2.0) / scale;
        wavelet[n] = normalization * std::exp(-(1/fb) * t * t) * std::exp(std::complex<double>(0, 2 * M_PI * f0 * t));
    }
    return wavelet;
}


std::vector< Eigen::VectorXcd> SignalData::CWT(const Eigen::VectorXcd& signal,const std::vector<double>& scales){
    int n = signal.size();
    std::vector<Eigen::VectorXcd> coeffs;
    for(const auto scale : scales){
      const auto morlet = MorletWavelet(n, scale);
      auto coef = FFTConvolve(signal, morlet);
      coeffs.push_back(coef);
    }

  return coeffs;
}

std::vector< Eigen::VectorXcd> SignalData::CWTEigen(const Eigen::VectorXcd& signal,const std::vector<double>& scales){
    int n = signal.size();
    std::vector<Eigen::VectorXcd> coeffs;
    for(const auto scale : scales){
      const auto morlet = MorletWavelet(n, scale);
      auto coef = FFTconvolveEigen(signal, morlet,true);
      coeffs.push_back(coef);
    }

  return coeffs;
}
