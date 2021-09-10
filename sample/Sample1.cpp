#include <cpp_filters/FirstOrderLowpassFilter.h>
#include <cpp_filters/IIRFilter.h>
#include <cpp_filters/TwoPointInterpolator.h>

int main(void){
   cpp_filters::FirstOrderLowPassFilter<double> firstOrderLowPassFilter(1.0,0.1,0.0);
  std::cout << "FirstOrderLowPassFilter"<<std::endl;
  std::cout << 0.0 << "\t" << firstOrderLowPassFilter.passFilter(0.0)<<std::endl;
  for(double t=0.1;t<=2.0;t+=0.1){
    std::cout << t << "\t" << firstOrderLowPassFilter.passFilter(1.0)<<std::endl;
  }

   cpp_filters::IIRFilter<double> iirFilter;
  iirFilter.setParameterAsBiquad(1.0,0.5,10,0.0);
  std::cout << "iirFilter"<<std::endl;
  std::cout << 0.0 << "\t" << iirFilter.passFilter(0.0)<<std::endl;
  for(double t=0.1;t<=2.0;t+=0.1){
    std::cout << t << "\t" << iirFilter.passFilter(1.0)<<std::endl;
  }

   cpp_filters::TwoPointInterpolator<double> twoPointInterpolator(0.0,0.0,0.0, cpp_filters::LINEAR);
  twoPointInterpolator.setGoal(1.0,1.0);
  std::cout << "twoPointInterpolator LINEAR"<<std::endl;
  for(double t=0.1;t<2.0;t+=0.1){
    double x;
    twoPointInterpolator.get(x,0.1);
    std::cout << t << "\t" << x<<std::endl;
  }

  twoPointInterpolator.reset(0.0);
  twoPointInterpolator.setInterpolationMode( cpp_filters::HOFFARBIB);
  twoPointInterpolator.setGoal(1.0,1.0);
  std::cout << "twoPointInterpolator HOFFARBIB"<<std::endl;
  for(double t=0.1;t<2.0;t+=0.1){
    double x;
    twoPointInterpolator.get(x,0.1);
    std::cout << t << "\t" << x<<std::endl;
  }

   cpp_filters::TwoPointInterpolatorSO3 twoPointInterpolatorSO3(Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),Eigen::Vector3d::Zero(), cpp_filters::HOFFARBIB);
  twoPointInterpolatorSO3.setGoal(Eigen::Matrix3d(Eigen::AngleAxisd(1.7,Eigen::Vector3d::UnitY())),1.0);
  std::cout << "twoPointInterpolatorSO3 HOFFARBIB"<<std::endl;
  for(double t=0.1;t<=1.0;t+=0.1){
    Eigen::Matrix3d x;
    Eigen::Vector3d v;
    twoPointInterpolatorSO3.get(x,v,0.1);
    std::cout << t << std::endl;
    std::cout << x << std::endl;
    std::cout << v.transpose() <<std::endl;
  }

  twoPointInterpolatorSO3.reset(Eigen::Matrix3d::Identity());
  twoPointInterpolatorSO3.setGoal(Eigen::Matrix3d(Eigen::AngleAxisd(1.5708,Eigen::Vector3d::UnitZ())*Eigen::AngleAxisd(1.7,Eigen::Vector3d::UnitY())),Eigen::Vector3d::UnitY(),Eigen::Vector3d::UnitY(),1.0);
  std::cout << "twoPointInterpolatorSO3 HOFFARBIB with vel acc"<<std::endl;
  for(double t=0.1;t<=1.0;t+=0.1){
    Eigen::Matrix3d x;
    Eigen::Vector3d v, a;
    twoPointInterpolatorSO3.get(x,v,a,0.1);
    std::cout << t << std::endl;
    std::cout << x << std::endl;
    std::cout << v.transpose() <<std::endl;
    std::cout << a.transpose() <<std::endl;
  }

}
