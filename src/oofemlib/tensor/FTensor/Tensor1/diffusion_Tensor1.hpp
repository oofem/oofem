/* Computes 2*del^2 of a Tensor1_ptr but uses diagonal derivatives for
   half of it.  Yields a Tensor1. */

#pragma once

namespace FTensor
{
  template <class T, int Dim, char i> class diffusion_Tensor1
  {
    const Tensor1<T *, Dim> &a;
    const int di, dj, dk;
    const double dx;

  public:
    typename promote<T, double>::V operator()(const int N1) const
    {
      return ((*(a.ptr(N1) + di) - 2 * a(N1) + *(a.ptr(N1) - di))
              + (*(a.ptr(N1) + dj) - 2 * a(N1) + *(a.ptr(N1) - dj))
              + (*(a.ptr(N1) + dk) - 2 * a(N1) + *(a.ptr(N1) - dk))
              + ((*(a.ptr(N1) + di + dj) + *(a.ptr(N1) + di - dj)
                  + *(a.ptr(N1) - di + dj) + *(a.ptr(N1) - di - dj)
                  - 4 * a(N1))
                 + (*(a.ptr(N1) + di + dk) + *(a.ptr(N1) + di - dk)
                    + *(a.ptr(N1) - di + dk) + *(a.ptr(N1) - di - dk)
                    - 4 * a(N1))
                 + (*(a.ptr(N1) + dj + dk) + *(a.ptr(N1) + dj - dk)
                    + *(a.ptr(N1) - dj + dk) + *(a.ptr(N1) - dj - dk)
                    - 4 * a(N1)))
                  / (std::sqrt(2.0)))
             * dx * dx;
    }
    diffusion_Tensor1(const Tensor1<T *, Dim> &A, const int Di, const int Dj,
                      const int Dk, const double Dx)
        : a(A), di(Di), dj(Dj), dk(Dk), dx(Dx)
    {}
  };

  template <class T, int Dim, char i>
  const Tensor1_Expr<const diffusion_Tensor1<T, Dim, i>,
                     typename promote<T, double>::V, Dim, i>
  diffusion(const Tensor1<T *, Dim> &a, const Index<i, Dim> index1,
            const int &di, const int &dj, const int &dk, const double &dx)
  {
    using Tensor_Expr = diffusion_Tensor1<T, Dim, i>;
    return Tensor1_Expr<Tensor_Expr, typename promote<T, double>::V, Dim, i>(
      Tensor_Expr(a, di, dj, dk, dx));
  }
}
