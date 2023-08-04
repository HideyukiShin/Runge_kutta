#include <rkm.hpp>
#include <matplotlibcpp.h>
#include <Eigen/Core>

Eigen::Matrix<double, 2, 1>
  Func(Eigen::Matrix<double, 2, 1> x)
  {
    Eigen::Matrix<double, 2,2> A;
    A << 0, 1, -0.1, -0.3;
    return A*x;
  }

  int main()
  {
    Eigen::Matrix<double, 2, 1> x{};
    x << 0.3, 0.7;
    auto rkm = Runge_Kutta<2>(x);
    rkm.run(Func);
    rkm.draw();
    rkm.draw_state();
    plt::show();
  }
