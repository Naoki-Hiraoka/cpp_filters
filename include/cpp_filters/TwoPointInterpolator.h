#ifndef __CPP_FILTERS_TWOPOINTINTERPOLATOR_H__
#define __CPP_FILTERS_TWOPOINTINTERPOLATOR_H__

#include <string>
#include <Eigen/Eigen>

namespace cpp_filters {
  typedef enum {LINEAR, HOFFARBIB,QUINTICSPLINE,CUBICSPLINE} interpolation_mode;

  template<typename T1, typename T2> class TwoPointInterpolatorBase
  {
    // interpolator class is to interpolate from current value to goal value considering position, velocities, and accelerations.
    //   Two status : empty or not
    //                Interpolator interpolates based on remaining time (remain_t) and pushes value to queue (q, dq, ddq).
    //                Users can get interpolated results from queue (q, dq, ddq).
    //                If remain_t <= 0 and queue is empty, interpolator is "empty", otherwise "not empty".
    //                This is related with isEmpty() function.
    //   Setting goal value : setGoal(), go(), and load()
    //   Getting current value : get()
    //   Resetting current value : set()
    //   Interpolate : interpolate()
  public:
    TwoPointInterpolatorBase(const T1& init_x, const T2& init_v, const T2& init_a, interpolation_mode imode=HOFFARBIB) :
      imode_(imode)
    {
      this->reset(init_x,init_v,init_a);
    }
    // Getter function.
    void get(T1& x, double dt=0.0) {
      T2 v, a;
      get(x, v, a, dt);
    }
    void get(T1& x, T2& v, double dt=0.0) {
      T2 a;
      get(x, v, a, dt);
    }
    void get(T1& x, T2& v, T2& a, double dt=0.0) {
      current_time_ += dt;
      if (current_time_ < 0.0) current_time_ = 0.0;
      if (current_time_ > goal_time_) current_time_ = goal_time_;
      double t = current_time_;
      this->getImpl(x,v,a,t);
    }
    // Reset current value.
    void reset(const T1& x) {
      T2 tmp;
      this->reset(x,tmp*0,tmp*0);
    }
    void reset(const T1& x, const T2& v) {
      T2 tmp;
      this->reset(x,v,tmp*0);
    }
    void reset(const T1& x, const T2& v, const T2& a)
    {
      goal_time_ = 0.0;
      current_time_ = 0.0;
      this->resetImpl(x,v,a);
    }
    // Stop to current value
    void clear() {
      T1 x;
      T2 v, a;
      this->get(x,v,a,0.0);
      this->reset(x,v*0,a*0);
    }
    bool isEmpty() {
      return current_time_ == goal_time_;
    }
    double remain_time() {
      return goal_time_ - current_time_;
    }
    bool setInterpolationMode (interpolation_mode i_mode){
      if (i_mode != LINEAR && i_mode != HOFFARBIB &&
          i_mode != QUINTICSPLINE && i_mode != CUBICSPLINE) return false;
      imode_ = i_mode;
      return true;
    };
    // Set goal
    void setGoal(const T1& goalx, double t) {
      T2 tmp;
      this->setGoal(goalx,tmp*0,tmp*0,t);
    }
    void setGoal(const T1& goalx, const T2& goalv, double t) {
      T2 tmp;
      this->setGoal(goalx,goalv,tmp*0,t);
    }
    void setGoal(const T1& goalx, const T2& goalv, const T2& goala, double t) {
      if(t == 0.0) {
        this->reset(goalx,goalv,goala);
        return;
      }

      T1 x;
      T2 v, a;
      this->get(x,v,a,0.0);
      this->current_time_ = 0.0;
      this->goal_time_ = t;

      this->setGoalImpl(x,v,a,goalx,goalv,goala,t);
    }
    std::string& name() { return name_; };
    const std::string& name() const { return name_; };

  protected:
    void getImpl(T1& x, T2& v, T2& a, double t);
    void resetImpl(const T1& x, const T2& v, const T2& a);
    void setGoalImpl(const T1& startx, const T2& startv, const T2& starta, const T1& goalx, const T2& goalv, const T2& goala, double t);
    void calcCoeff(const T2& startx, const T2& startv, const T2& starta, const T2& goalx, const T2& goalv, const T2& goala, double t){
      T2 A,B,C;
      switch(imode_){
      case LINEAR:
        a0_=startx;
        a1_=(goalx-startx)/t;
        a2_*=0;
        a3_*=0;
        a4_*=0;
        a5_*=0;
        break;
      case HOFFARBIB:
        A=(goalx-(startx+startv*t+(starta/2.0)*t*t))/(t*t*t);
        B=(goalv-(startv+starta*t))/(t*t);
        C=(goala-starta)/t;

        a0_=startx;
        a1_=startv;
        a2_=starta/2.0;
        a3_=10*A-4*B+0.5*C;
        a4_=(-15*A+7*B-C)/t;
        a5_=(6*A-3*B+0.5*C)/(t*t);
        break;
      case QUINTICSPLINE:
        a0_=startx;
        a1_=startv;
        a2_=0.5*starta;
        a3_=(-20*startx + 20*goalx - 3*starta*t*t + goala*t*t -
             12*startv*t - 8*goalv*t) / (2*t*t*t);
        a4_=(30*startx - 30*goalx + 3*starta*t*t - 2*goala*t*t +
             16*startv*t + 14*goalv*t) / (2*t*t*t*t);
        a5_=(-12*startx + 12*goalx - starta*t*t + goala*t*t -
             6*startv*t - 6*goalv*t) / (2*t*t*t*t*t);
        break;
      case CUBICSPLINE:
        a0_=startx;
        a1_=startv;
        a2_=(-3*startx + 3*goalx - 2*startv*t - goalv*t) / (t*t);
        a3_=( 2*startx - 2*goalx + startv*t + goalv*t) / (t*t*t);
        a4_*=0;
        a5_*=0;
        break;
      }
    }
    void calcPolynomial(T2& x, T2& v, T2& a, double t){
      x=a0_+a1_*t+a2_*t*t+a3_*t*t*t+a4_*t*t*t*t+a5_*t*t*t*t*t;
      v=a1_+2*a2_*t+3*a3_*t*t+4*a4_*t*t*t+5*a5_*t*t*t*t;
      a=2*a2_+6*a3_*t+12*a4_*t*t+20*a5_*t*t*t;
    }
    // Current interpolation mode
    interpolation_mode imode_;
    double goal_time_;
    double current_time_;
    // Coefficients for interpolation polynomials.
    T2 a0_, a1_, a2_, a3_, a4_, a5_;
    T1 startx_;
    // Interpolator name
    std::string name_;
  };


  // for Euclid
  template<typename T1, typename T2>
  void TwoPointInterpolatorBase<T1,T2>::getImpl(T1& x, T2& v, T2& a, double t) {
    this->calcPolynomial(x,v,a,t);
  }
  template<typename T1, typename T2>
  void TwoPointInterpolatorBase<T1,T2>::resetImpl(const T1& x, const T2& v, const T2& a) {
    this->startx_ = x;
    this->a0_ = x;
    this->a1_ = v;
    this->a2_ = a/2;
    this->a3_ = x*0;
    this->a4_ = x*0;
    this->a5_ = x*0;
  }
  template<typename T1, typename T2>
  void TwoPointInterpolatorBase<T1,T2>::setGoalImpl(const T1& startx, const T2& startv, const T2& starta, const T1& goalx, const T2& goalv, const T2& goala, double t) {
    this->calcCoeff(startx, startv, starta, goalx, goalv, goala, t);
  }
  template<typename T> using TwoPointInterpolator = TwoPointInterpolatorBase<T,T>;


  // for SO3. v and a are represented in local frame.
  template<>
  void TwoPointInterpolatorBase<Eigen::Matrix3d,Eigen::Vector3d>::getImpl(Eigen::Matrix3d& x, Eigen::Vector3d& v, Eigen::Vector3d& a, double t);
  template<>
  void TwoPointInterpolatorBase<Eigen::Matrix3d,Eigen::Vector3d>::resetImpl(const Eigen::Matrix3d& x, const Eigen::Vector3d& v, const Eigen::Vector3d& a);
  template<>
  void TwoPointInterpolatorBase<Eigen::Matrix3d,Eigen::Vector3d>::setGoalImpl(const Eigen::Matrix3d& startx, const Eigen::Vector3d& startv, const Eigen::Vector3d& starta, const Eigen::Matrix3d& goalx, const Eigen::Vector3d& goalv, const Eigen::Vector3d& goala, double t);
  using TwoPointInterpolatorSO3 = TwoPointInterpolatorBase<Eigen::Matrix3d,Eigen::Vector3d>;

}

#endif
