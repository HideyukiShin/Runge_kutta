#ifndef rkm_
#define rkm_

#include "matplotlibcpp.h"
#include "Eigen/Core"
#include <vector>
#include <map>

namespace plt = matplotlibcpp;

template <size_t N>
class Runge_Kutta
{
   using X = Eigen::Matrix<double, N, 1>;
public:
 
  constexpr explicit Runge_Kutta(X init)
      : state{init}, t{0}, states{}, time{} {}

  void draw() const
  {
    for(size_t i = 0; i<N; ++i)
    {
      plt::figure(i+1);
      plt::plot(time, states[i]);
    }
  }

  void draw_state() const 
  {
    if(N!=2) 
    {
      std::cout << "this method is only applied when N=2";
      return;
    }
    plt::figure();
    plt::plot(states[0], states[1]);
  }

  size_t size() const 
  {
    return time.size();
  }

  void record()
  {
    for (size_t i = 0; i < N; ++i)
    {
      states[i].push_back(state(i));
    }
    time.push_back(t);
  }

  template <class F>
  void run(F &&f)
  {
    while (t < T_end)
    {
      auto x = state;
      record();
      auto k_1 = f(x);
      auto k_2 = f(x + 0.5 * T * k_1);
      auto k_3 = f(x + 0.5 * T * k_2);
      auto k_4 = f(x + T * k_3);
      state = x + T / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
      t += T;
    }
  }

private:
  
  X state;
  double t;
  std::array<std::vector<double>, N> states;
  std::vector<double> time;
  double const T = 1e-3;
  double const T_end = 100;
};

#endif