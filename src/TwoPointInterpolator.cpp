#include <cpp_filters/TwoPointInterpolator.h>

namespace cpp_filters {
  inline Eigen::Matrix3d hat(const Eigen::Vector3d& x) {
    Eigen::Matrix3d M;
    M <<  0.0, -x(2),   x(1),
      x(2),   0.0,  -x(0),
      -x(1),  x(0),   0.0;
    return M;
  }

  // Smooth Attitude Interpolation
  // https://github.com/scipy/scipy/files/2932755/attitude_interpolation.pdf
  template<>
  void TwoPointInterpolatorBase<Eigen::Matrix3d,Eigen::Vector3d>::getImpl(Eigen::Matrix3d& x, Eigen::Vector3d& v, Eigen::Vector3d& a, double t) {
    Eigen::Vector3d theta, dtheta, ddtheta;
    this->calcPolynomial(theta,dtheta,ddtheta,t);
    double th = theta.norm();
    Eigen::Matrix3d thetaX = hat(theta);

    Eigen::Matrix3d Ainv = Eigen::Matrix3d::Identity();
    if(th>0.0){
      Ainv =
        Eigen::Matrix3d::Identity()
        - (1-std::cos(th))/std::pow(th,2) * thetaX
        + (th - std::sin(th))/std::pow(th,3) * thetaX * thetaX;
    }
    Eigen::Vector3d dAinv_dtheta = Eigen::Vector3d::Zero();
    if(th>0.0){
      dAinv_dtheta =
        - (th*std::sin(th)+2*(std::cos(th)-1))/std::pow(th,4) * theta.dot(dtheta) * theta.cross(dtheta)
        - (2*th+th*std::cos(th)-3*std::sin(th))/std::pow(th,5) * theta.dot(dtheta) * theta.cross(theta.cross(dtheta))
        + (th-std::sin(th))/std::pow(th,3) * dtheta.cross(theta.cross(dtheta));
    }

    // 単純に3x3行列の空間でRを積算していると、だんだん数値誤差によってユニタリ行列でなくなってしまう
    if(th>0.0) x = Eigen::Matrix3d(Eigen::AngleAxisd(this->startx_) * Eigen::AngleAxisd(th,theta.normalized()));
    else x = this->startx_;
    v = Ainv * dtheta;
    a = Ainv * ddtheta + dAinv_dtheta;
  }
  template<>
  void TwoPointInterpolatorBase<Eigen::Matrix3d,Eigen::Vector3d>::resetImpl(const Eigen::Matrix3d& x, const Eigen::Vector3d& v, const Eigen::Vector3d& a) {
    this->startx_ = x;
    this->a0_ = Eigen::Vector3d::Zero();
    this->a1_ = v;
    this->a2_ = a/2;
    this->a3_ = Eigen::Vector3d::Zero();
    this->a4_ = Eigen::Vector3d::Zero();
    this->a5_ = Eigen::Vector3d::Zero();
  }
  template<>
  void TwoPointInterpolatorBase<Eigen::Matrix3d,Eigen::Vector3d>::setGoalImpl(const Eigen::Matrix3d& startx, const Eigen::Vector3d& startv, const Eigen::Vector3d& starta, const Eigen::Matrix3d& goalx, const Eigen::Vector3d& goalv, const Eigen::Vector3d& goala, double t) {
    this->startx_ = startx;

    Eigen::Vector3d starttheta = Eigen::Vector3d::Zero();
    Eigen::AngleAxisd angleaxis(startx.transpose()*goalx);
    Eigen::Vector3d goaltheta = angleaxis.angle() * angleaxis.axis();
    double th = goaltheta.norm();
    Eigen::Matrix3d thetaX = hat(goaltheta);
    Eigen::Vector3d startdtheta = startv;
    Eigen::Matrix3d A = Eigen::Matrix3d::Identity();
    if(th > 0.0){
      A =
        Eigen::Matrix3d::Identity()
        + 0.5*thetaX
        + (1-th/2*std::cos(th/2)/std::sin(th/2))/std::pow(th,2) * thetaX * thetaX;
    }
    Eigen::Vector3d goaldtheta = A*goalv;
    double dth = goaldtheta.norm();
    Eigen::Matrix3d dthetaX = hat(goaldtheta);
    Eigen::Vector3d startddtheta = starta;
    Eigen::Matrix3d dA = 0.5*thetaX;
    if(th > 0.0){
      dA =
        0.5*dthetaX
        + (dth*std::cos(th/2)/std::sin(th/2)/(2*std::pow(th,2)) + dth/(4*th*std::pow(std::sin(th/2),2)) - 2*dth/std::pow(th,3)) * thetaX * thetaX
        + (1-th/2*std::cos(th/2)/std::sin(th/2))/std::pow(th,2) * (dthetaX*thetaX + thetaX*dthetaX);
    }
    Eigen::Vector3d goalddtheta = dA * goalv + A * goala;

    this->calcCoeff(starttheta, startdtheta, startddtheta, goaltheta, goaldtheta, goalddtheta, t);
  }

}
