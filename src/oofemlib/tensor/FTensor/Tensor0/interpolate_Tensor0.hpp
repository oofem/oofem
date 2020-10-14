/* Interpolates to (i0+distance[0], j0+distance[1], k0+distance[2]),
   yielding a double.  (i0,j0,k0) are 3D array coordinates,
   conjugate==1-distance, and (di,dj,dk) are the stride of the array
   coordinates.  It is assumed that the Tensor0_ptr is zero
   centered. */

#pragma once

namespace FTensor
{
  template <class T>
  typename promote<T, double>::V
  interpolate(const Tensor0<T *> &a, const int &di, const int &dj,
              const int &dk, const int &i0, const int &j0, const int &k0,
              const double distance[3], const double conjugate[3])
  {
    return conjugate[0] * conjugate[1] * conjugate[2]
             * (*(&a + di * i0 + dj * j0 + dk * k0))
           + distance[0] * conjugate[1] * conjugate[2]
               * (*(&a + di * (i0 + 1) + dj * (j0) + dk * (k0)))
           + conjugate[0] * distance[1] * conjugate[2]
               * (*(&a + di * (i0) + dj * (j0 + 1) + dk * (k0)))
           + distance[0] * distance[1] * conjugate[2]
               * (*(&a + di * (i0 + 1) + dj * (j0 + 1) + dk * (k0)))
           + conjugate[0] * conjugate[1] * distance[2]
               * (*(&a + di * (i0) + dj * (j0) + dk * (k0 + 1)))
           + distance[0] * conjugate[1] * distance[2]
               * (*(&a + di * (i0 + 1) + dj * (j0) + dk * (k0 + 1)))
           + conjugate[0] * distance[1] * distance[2]
               * (*(&a + di * (i0) + dj * (j0 + 1) + dk * (k0 + 1)))
           + distance[0] * distance[1] * distance[2]
               * (*(&a + di * (i0 + 1) + dj * (j0 + 1) + dk * (k0 + 1)));
  }
}
