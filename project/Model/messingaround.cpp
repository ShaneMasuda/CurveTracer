#include <iostream>
#include <sstream>
#include "alglib/interpolation.h"

int main()
{
  double x[] = {1, 2, 3, 4, 5};
  double y[] = {1, 4, 9, 16, 25};
  alglib::real_1d_array w;
  alglib::real_1d_array z;
  w.setcontent(5, x);
  z.setcontent(5, y);
  alglib::ae_int_t info;
  alglib::barycentricinterpolant b;

  alglib::polynomialfitreport report;
  alglib::polynomialfit(w, z, 5, alglib::ae_int_t(3), info, b, report);
  alglib::real_1d_array curve;
  alglib::xparams degree;
  alglib::polynomialbar2pow(b, curve, degree);
  // std::cout << curve.tostring(3) << std::endl;
  std::ostringstream oss;
  oss << std::fixed;
  oss << 234.234234234234;
  std::cout << oss.str() << std::endl;
  oss.str("");
  oss << 444.444444444;
  std::cout << oss.str() << std::endl;
}
