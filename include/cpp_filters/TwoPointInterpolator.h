#ifndef __CPP_FILTERS_INTERPOLATOR_H__
#define __CPP_FILTERS_INTERPOLATOR_H__

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
    interpolator(const T& init_x, const T& init_v, const T& init_a, interpolation_mode imode=HOFFARBIB) :
      imode_(imode),
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
    void get(double *x, double *v, double *a, double dt=0.0) {
      current_time_ += dt;
      if (current_time_ < 0.0) current_time_ = 0.0;
      if (current_time_ > goal_time_) current_time_ = goal_time_;
      double t = current_time_;
      x=a0+a1*t+a2*t*t+a3*t*t*t+a4*t*t*t*t+a5*t*t*t*t*t;
      v=a1+2*a2*t+3*a3*t*t+4*a4*t*t*t+5*a5*t*t*t*t;
      a=2*a2+6*a3*t+12*a4*t*t+20*a5*t*t*t;
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
      a0_ = init_x;
      a1_ = init_v;
      a2_ = init_a/2;
      a3_ = init_x*0;
      a4_ = init_x*0;
      a5_ = init_x*0;
    }
    // Stop to current value
    void clear() {
      T x, v, a;
      this->get(x,v,a,0.0);
      this->reset(x,v,a);
    }
    bool isEmpty() {
      return current_time_ == goal_time_;
    }
    double remain_time() {
      return goal_time_ - current_time_;
    }
    bool setInterpolationMode (interpolation_mode i_mode_){
      if (i_mode_ != LINEAR && i_mode_ != HOFFARBIB &&
          i_mode_ != QUINTICSPLINE && i_mode_ != CUBICSPLINE) return false;
      imode = i_mode_;
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
      switch(imode){
      case LINEAR:
        a0=x;
        a1=(goalx-x)/t;
        a2*=0;
        a3*=0;
        a4*=0;
        a5*=0;
        break;
      case HOFFARBIB:
        A=(goalx-(x+v*t+(a/2.0)*t*t))/(t*t*t);
        B=(goalv-(v+a*t))/(t*t);
        C=(goala-a)/t;

        a0=x;
        a1=v;
        a2=a/2.0;
        a3=10*A-4*B+0.5*C;
        a4=(-15*A+7*B-C)/t;
        a5=(6*A-3*B+0.5*C)/(t*t);
        break;
      case QUINTICSPLINE:
        a0=x;
        a1=v;
        a2=0.5*a;
        a3=(-20*x + 20*goalx - 3*a*t*t + goala*t*t -
               12*v*t - 8*goalv*t) / (2*t*t*t);
        a4=(30*x - 30*goalx + 3*a*t*t - 2*goala*t*t +
               16*v*t + 14*goalv*t) / (2*t*t*t*t);
        a5=(-12*x + 12*goalx - a*t*t + goala*t*t -
               6*v*t - 6*goalv*t) / (2*t*t*t*t*t);
        break;
      case CUBICSPLINE:
        a0=x;
        a1=v;
        a2=(-3*x + 3*goalx - 2*v*t - goalv*t) / (t*t);
        a3=( 2*x - 2*goalx + v*t + goalv*t) / (t*t*t);
        a4*=0;
        a5*=0;
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
