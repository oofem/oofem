/* Traits that allow auto-promotion of int's to double's, double's to
   complex, etc.  Shamelessly stolen from

   http://extreme.indiana.edu/~tveldhui/papers/techniques/

   For now, only int, double, and complex<double> are covered.  If you
   want more, just insert a DECLARE_PROMOTE(A,B,C), where A and B are
   the two types, and C is what they should be coerced to.  */

#pragma once

namespace FTensor
{
  template <class T1, class T2> class promote
  {
  public:
    using V = T1;
  };

#define DECLARE_PROMOTE(A, B, C)                                              \
  template <> class promote<A, B>                                             \
  {                                                                           \
  public:                                                                     \
    using V = C;                                                              \
  }

  DECLARE_PROMOTE(int, double, double);
  DECLARE_PROMOTE(double, int, double);
  DECLARE_PROMOTE(int, std::complex<double>, std::complex<double>);
  DECLARE_PROMOTE(std::complex<double>, int, std::complex<double>);
  DECLARE_PROMOTE(double, std::complex<double>, std::complex<double>);
  DECLARE_PROMOTE(std::complex<double>, double, std::complex<double>);

#undef DECLARE_PROMOTE
}
