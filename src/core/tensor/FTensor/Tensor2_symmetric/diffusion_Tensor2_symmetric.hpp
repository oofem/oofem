/* Computes 2*del^2 of a Tensor2_symmetric<T*,Dim> but uses diagonal
   derivatives for half of it.  Yields a Tensor2_symmetric. */

#pragma once

namespace FTensor
{
  template <class T, int Dim, char i, char j> class diffusion_Tensor2_symmetric
  {
    const Tensor2_symmetric<T *, Dim> &a;
    const int di, dj, dk;
    const double dx;

  public:
    typename promote<T, double>::V operator()(const int N1, const int N2) const
    {
      return ((*(a.ptr(N1, N2) + di) - 2 * a(N1, N2) + *(a.ptr(N1, N2) - di))
              + (*(a.ptr(N1, N2) + dj) - 2 * a(N1, N2) + *(a.ptr(N1, N2) - dj))
              + (*(a.ptr(N1, N2) + dk) - 2 * a(N1, N2) + *(a.ptr(N1, N2) - dk))
              + ((*(a.ptr(N1, N2) + di + dj) + *(a.ptr(N1, N2) + di - dj)
                  + *(a.ptr(N1, N2) - di + dj) + *(a.ptr(N1, N2) - di - dj)
                  - 4 * a(N1, N2))
                 + (*(a.ptr(N1, N2) + di + dk) + *(a.ptr(N1, N2) + di - dk)
                    + *(a.ptr(N1, N2) - di + dk) + *(a.ptr(N1, N2) - di - dk)
                    - 4 * a(N1, N2))
                 + (*(a.ptr(N1, N2) + dj + dk) + *(a.ptr(N1, N2) + dj - dk)
                    + *(a.ptr(N1, N2) - dj + dk) + *(a.ptr(N1, N2) - dj - dk)
                    - 4 * a(N1, N2)))
                  / (std::sqrt(2.0)))
             * dx * dx;
    }
    diffusion_Tensor2_symmetric(const Tensor2_symmetric<T *, Dim> &A,
                                const int Di, const int Dj, const int Dk,
                                const double Dx)
        : a(A), di(Di), dj(Dj), dk(Dk), dx(Dx)
    {}
  };

  template <class T, int Dim, char i, char j>
  const Tensor2_symmetric_Expr<const diffusion_Tensor2_symmetric<T, Dim, i, j>,
                               typename promote<T, double>::V, Dim, i, j>
  diffusion(const Tensor2_symmetric<T *, Dim> &a, const Index<i, Dim> index1,
            const Index<j, Dim> index2, const int &di, const int &dj,
            const int &dk, const double &dx)
  {
    using Tensor_Expr = diffusion_Tensor2_symmetric<T, Dim, i, j>;
    return Tensor2_symmetric_Expr<Tensor_Expr, typename promote<T, double>::V,
                                  Dim, i, j>(Tensor_Expr(a, di, dj, dk, dx));
  }
}
