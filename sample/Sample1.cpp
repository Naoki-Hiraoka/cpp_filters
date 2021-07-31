#include <cpp_filters/FirstOrderLowpassFilter.h>
#include <cpp_filters/IIRFilter.h>
#include <cpp_filters/TwoPointInterpolator.h>

int main(void){
  filters::FirstOrderLowPassFilter<double> firstOrderLowPassFilter(1.0,0.1,0.0);
  std::cout << "FirstOrderLowPassFilter"<<std::endl;
  std::cout << 0.0 << "\t" << firstOrderLowPassFilter.passFilter(0.0)<<std::endl;
  for(double t=0.1;t<=2.0;t+=0.1){
    std::cout << t << "\t" << firstOrderLowPassFilter.passFilter(1.0)<<std::endl;
  }

  filters::IIRFilter<double> iirFilter;
  iirFilter.setParameterAsBiquad(1.0,0.5,10,0.0);
  std::cout << "iirFilter"<<std::endl;
  std::cout << 0.0 << "\t" << iirFilter.passFilter(0.0)<<std::endl;
  for(double t=0.1;t<=2.0;t+=0.1){
    std::cout << t << "\t" << iirFilter.passFilter(1.0)<<std::endl;
  }

  filters::TwoPointInterpolator<double> twoPointInterpolator(0.0,0.0,0.0,filters::TwoPointInterpolator<double>::LINEAR);
  twoPointInterpolator.setGoal(1.0,1.0);
  std::cout << "twoPointInterpolator LINEAR"<<std::endl;
  for(double t=0.1;t<2.0;t+=0.1){
    double x;
    twoPointInterpolator.get(x,0.1);
    std::cout << t << "\t" << x<<std::endl;
  }

  twoPointInterpolator.reset(0.0);
  twoPointInterpolator.setInterpolationMode(filters::TwoPointInterpolator<double>::HOFFARBIB);
  twoPointInterpolator.setGoal(1.0,1.0);
  std::cout << "twoPointInterpolator HOFFARBIB"<<std::endl;
  for(double t=0.1;t<2.0;t+=0.1){
    double x;
    twoPointInterpolator.get(x,0.1);
    std::cout << t << "\t" << x<<std::endl;
  }

}
