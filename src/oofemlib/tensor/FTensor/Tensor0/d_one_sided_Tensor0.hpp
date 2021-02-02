/* Takes a one-sided derivative of a Tensor0 in a particular
   direction, yielding a typename promote<T,double>::V. */

#pragma once

namespace FTensor
{
  template <class T>
  typename promote<T, double>::V
  d_one_sided(const Tensor0<T *> &a, const Number<0> n1, const int &di,
              const int &dj, const int &dk, const double &dx, const double &dy,
              const double &dz)
  {
    return (a - *(&a - di)) * dx;
  }

  template <class T>
  typename promote<T, double>::V
  d_one_sided(const Tensor0<T *> &a, const Number<1> n1, const int &di,
              const int &dj, const int &dk, const double &dx, const double &dy,
              const double &dz)
  {
    return (a - *(&a - dj)) * dy;
  }

  template <class T>
  typename promote<T, double>::V
  d_one_sided(const Tensor0<T *> &a, const Number<2> n1, const int &di,
              const int &dj, const int &dk, const double &dx, const double &dy,
              const double &dz)
  {
    return (a - *(&a - dk)) * dz;
  }
}
