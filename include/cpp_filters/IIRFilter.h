#ifndef __CPP_FILTERS_IIR_FILTER_H__
#define __CPP_FILTERS_IIR_FILTER_H__

#include <cmath>

namespace filters{
  /**
     Infinite Impulse Filter
     y[n] = sum(0, dim, ff_coeffs[i] * x[n - i]) + sum(1, dim, fb_coeffs[i] * y[n - i])
  */
  template <class T> class IIRFilter
  {
  public:
    /**
       \brief Constructor
    */
    IIRFilter(const std::string& error_prefix = "") :
      m_initialized(false),
      m_error_prefix(error_prefix)
    {
    }

    /**
       \brief Set parameters
       Y[n] = B[0] * X[n] + B[1] * X[n-1] + ... + B[dim] * X[n-dim] - A[1] * Y[n-1] ... - A[dim] * Y[n-dim]
       A[0] would be 1.0

       How to generete parameter by octave
       butterworth filter (dimension = 2, cutoff_freq = 8Hz)
       [B, A] = butter(2, 2 * 0.004 * 8) ;;; dimension=2, 2 * dt * cutoff_freq
    */
    bool setParameter(int dim, std::vector<double>& A, std::vector<double>& B, const T& initial_input) {
      m_dimension = dim;

      // init coefficients
      if((A.size() != dim && A.size() != dim + 1) || B.size() != dim + 1) {
        std::cout << "[" <<  m_error_prefix << "]" << "IIRFilter coefficients size error" << std::endl;
        return false;
      }

      // clear previous coefficients
      m_fb_coefficients.clear();
      m_ff_coefficients.clear();

      if (A.size() == dim) {
        m_fb_coefficients.push_back(1.0);
      }
      for(std::vector<double>::iterator it = A.begin(); it != A.end(); it++){
        if (it == A.begin()) {
          if( *it != 1.0 ) {
            std::cout << "[" <<  m_error_prefix << "]" << "IIRFilter : parameter A[0] is not 1.0 !!!" << std::endl;
          }
          m_fb_coefficients.push_back(*it);
        } else {
          m_fb_coefficients.push_back(- *it);
        }
      }
      for(std::vector<double>::iterator it = B.begin(); it != B.end(); it++){
        m_ff_coefficients.push_back(*it);
      }

      // init previous values
      this->reset(init_value);
      m_initialized = true;
      return true;
    }

    /**
       \brief Simple user interface of setParameter
       \param f_cutoff cut off frequency
       \param Q quality factor: 1/2 = no overshoot, 1/sqrt(2) = Butterworth
       \param hz sampling rate
    */
    bool setParameterAsBiquad(const double f_cutoff, const double Q, const double hz) {
      std::vector<double> fb_coeffs(3), ff_coeffs(3);
      const double omega = 2 * 3.14159265 * f_cutoff / hz;
      const double alpha = std::sin(omega) / (2 * Q);
      const double denom = 1 + alpha;
      fb_coeffs[0] = 1;
      fb_coeffs[1] = -2 * std::cos(omega) / denom;
      fb_coeffs[2] = (1 - alpha) / denom;
      ff_coeffs[0] = (1 - std::cos(omega)) / 2 / denom;
      ff_coeffs[1] = (1 - std::cos(omega)) / denom;
      ff_coeffs[2] = (1 - std::cos(omega)) / 2 / denom;
      return this->setParameter(2, fb_coeffs, ff_coeffs);
    }

    /**
     */
    void getParameter(int &dim, std::vector<double>&A, std::vector<double>& B) {
      dim = m_dimension;
      B.resize(m_ff_coefficients.size());
      std::copy(m_ff_coefficients.begin(), m_ff_coefficients.end(), B.begin());
      A.resize(0);
      for(std::vector<double>::iterator it = m_fb_coefficients.begin();
          it != m_fb_coefficients.end(); it++) {
        if (it == m_fb_coefficients.begin()) {
          A.push_back(*it);
        } else {
          A.push_back(- *it);
        }
      }
    }

    /**
     */
    void reset(const T& initial_input) {
      // y[n] = b[0]*w[n] + b[1]*w[n-1] + ... + b[m]*w[n-m] in DirectForm-II.
      // When n->inf, y[n]->initial_input and w[n], w[n-1], ..., w[n-m] -> w,
      // m_previous_values should preserve w
      double sum_ff_coeffs = std::accumulate(m_ff_coefficients.begin(), m_ff_coefficients.end(), 0.0);
      T reset_val = initial_input / sum_ff_coeffs;
      m_previous_values.assign(m_dimension, reset_val);
    }

    /**
       \brief passFilter
    */
    double passFilter(const T& input) {
      // IIRFilter implementation based on DirectForm-II.
      // Cf. https://en.wikipedia.org/wiki/Digital_filter
      if (! m_initialized) {
        return 0.0;
      }
      T feedback, filtered;
      // calcurate retval
      feedback = m_fb_coefficients[0] * input;
      for (int i = 0; i < m_dimension; i++) {
        feedback += m_fb_coefficients[i + 1] * m_previous_values[i];
      }
      filtered = m_ff_coefficients[0] * feedback;
      for(int i = 0; i < m_dimension; i++) {
        filtered += m_ff_coefficients[i + 1] * m_previous_values[i];
      }
      // update previous values
      m_previous_values.push_front(feedback);
      m_previous_values.pop_back();

      return filtered;
    }

  private:
    int m_dimension;
    std::vector<double> m_fb_coefficients; // fb parameters (dim must be m_dimension + 1, m_fb_coefficients[0] would be 1.0)
    std::vector<double> m_ff_coefficients; // ff parameters (dim must be m_dimension + 1)
    std::deque<T> m_previous_values;

    bool m_initialized;
    std::string m_error_prefix;
  };

}
#endif
