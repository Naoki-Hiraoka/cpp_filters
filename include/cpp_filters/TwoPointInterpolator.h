#ifndef __CPP_FILTERS_TWOPOINTINTERPOLATOR_H__
#define __CPP_FILTERS_TWOPOINTINTERPOLATOR_H__

#include <string>

namespace filters {
  template<typename T> class TwoPointInterpolator
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
    typedef enum {LINEAR, HOFFARBIB,QUINTICSPLINE,CUBICSPLINE} interpolation_mode;
    TwoPointInterpolator(const T& init_x, const T& init_v, const T& init_a, interpolation_mode imode=HOFFARBIB) :
      imode_(imode)
    {
      this->reset(init_x,init_v,init_a);
    }
    // Getter function.
    void get(T& x, double dt=0.0) {
      T v, a;
      get(x, v, a, dt);
    }
    void get(T& x, T& v, double dt=0.0) {
      T a;
      get(x, v, a, dt);
    }
    void get(T& x, T& v, T& a, double dt=0.0) {
      current_time_ += dt;
      if (current_time_ < 0.0) current_time_ = 0.0;
      if (current_time_ > goal_time_) current_time_ = goal_time_;
      double t = current_time_;
      x=a0_+a1_*t+a2_*t*t+a3_*t*t*t+a4_*t*t*t*t+a5_*t*t*t*t*t;
      v=a1_+2*a2_*t+3*a3_*t*t+4*a4_*t*t*t+5*a5_*t*t*t*t;
      a=2*a2_+6*a3_*t+12*a4_*t*t+20*a5_*t*t*t;
    }
    // Reset current value.
    void reset(const T& x) {
      this->reset(x,x*0,x*0);
    }
    void reset(const T& x, const T& v) {
      this->reset(x,v,x*0);
    }
    void reset(const T& x, const T& v, const T& a)
    {
      goal_time_ = 0.0;
      current_time_ = 0.0;
      a0_ = x;
      a1_ = v;
      a2_ = a/2;
      a3_ = x*0;
      a4_ = x*0;
      a5_ = x*0;
    }
    // Stop to current value
    void clear() {
      T x, v, a;
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
    void setGoal(const T& goalx, double t) {
      this->setGoal(goalx,goalx*0,goalx*0,t);
    }
    void setGoal(const T& goalx, const T& goalv, double t) {
      this->setGoal(goalx,goalv,goalx*0,t);
    }
    void setGoal(const T& goalx, const T& goalv, const T& goala, double t) {
      if(t == 0.0) {
        this->reset(goalx,goalv,goala);
        return;
      }

      T x, v, a;
      this->get(x,v,a,0.0);
      this->current_time_ = 0.0;
      this->goal_time_ = t;

      T A,B,C;
      switch(imode_){
      case LINEAR:
        a0_=x;
        a1_=(goalx-x)/t;
        a2_*=0;
        a3_*=0;
        a4_*=0;
        a5_*=0;
        break;
      case HOFFARBIB:
        A=(goalx-(x+v*t+(a/2.0)*t*t))/(t*t*t);
        B=(goalv-(v+a*t))/(t*t);
        C=(goala-a)/t;

        a0_=x;
        a1_=v;
        a2_=a/2.0;
        a3_=10*A-4*B+0.5*C;
        a4_=(-15*A+7*B-C)/t;
        a5_=(6*A-3*B+0.5*C)/(t*t);
        break;
      case QUINTICSPLINE:
        a0_=x;
        a1_=v;
        a2_=0.5*a;
        a3_=(-20*x + 20*goalx - 3*a*t*t + goala*t*t -
             12*v*t - 8*goalv*t) / (2*t*t*t);
        a4_=(30*x - 30*goalx + 3*a*t*t - 2*goala*t*t +
             16*v*t + 14*goalv*t) / (2*t*t*t*t);
        a5_=(-12*x + 12*goalx - a*t*t + goala*t*t -
             6*v*t - 6*goalv*t) / (2*t*t*t*t*t);
        break;
      case CUBICSPLINE:
        a0_=x;
        a1_=v;
        a2_=(-3*x + 3*goalx - 2*v*t - goalv*t) / (t*t);
        a3_=( 2*x - 2*goalx + v*t + goalv*t) / (t*t*t);
        a4_*=0;
        a5_*=0;
        break;
      }
    }
    std::string& name() { return name_; };
    const std::string& name() const { return name_; };

  private:
    // Current interpolation mode
    interpolation_mode imode_;
    double goal_time_;
    double current_time_;
    // Coefficients for interpolation polynomials.
    T a0_, a1_, a2_, a3_, a4_, a5_;
    // Interpolator name
    std::string name_;
  };
}

#endif
