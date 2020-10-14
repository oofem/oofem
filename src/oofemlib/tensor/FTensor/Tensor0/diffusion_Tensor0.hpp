/* Computes 2*del^2 of a Tensor0_ptr but uses diagonal derivatives for
   half of it. */

#pragma once

namespace FTensor
{
  template <class T>
  typename promote<T, double>::V
  diffusion(const Tensor0<T *> &a, const int &di, const int &dj, const int &dk,
            const double &dx)
  {
    return ((*(&a + di) - 2 * a + *(&a - di))
            + (*(&a + dj) - 2 * a + *(&a - dj))
            + (*(&a + dk) - 2 * a + *(&a - dk))
            + ((*(&a + di + dj) + *(&a + di - dj) + *(&a - di + dj)
                + *(&a - di - dj) - 4 * a)
               + (*(&a + di + dk) + *(&a + di - dk) + *(&a - di + dk)
                  + *(&a - di - dk) - 4 * a)
               + (*(&a + dj + dk) + *(&a + dj - dk) + *(&a - dj + dk)
                  + *(&a - dj - dk) - 4 * a))
                / (std::sqrt(2.0)))
           * dx * dx;
  }
}
