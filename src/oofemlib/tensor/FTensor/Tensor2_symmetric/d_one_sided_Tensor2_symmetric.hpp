/* Takes a one-sided derivative of a Tensor2_symmetric in a particular
   direction, yielding a Tensor2_symmetric. */

#pragma once

namespace FTensor
{
  template <class T, int Dim, char i, char j, int axis>
  class d_one_sided_Tensor2_symmetric
  {
    const Tensor2_symmetric<T *, Dim> &a;
    const int di, dj, dk;
    const double dx, dy, dz;

  public:
    typename promote<T, double>::V operator()(const int N1, const int N2) const
    {
      return axis == 0
               ? (a(N1, N2) - *(a.ptr(N1, N2) - di)) * dx
               : (axis == 1 ? (a(N1, N2) - *(a.ptr(N1, N2) - dj)) * dy
                            : (a(N1, N2) - *(a.ptr(N1, N2) - dk)) * dz);
    }
    d_one_sided_Tensor2_symmetric(const Tensor2_symmetric<T *, Dim> &A,
                                  const int Di, const int Dj, const int Dk,
                                  const double Dx, const double Dy,
                                  const double Dz)
        : a(A), di(Di), dj(Dj), dk(Dk), dx(Dx), dy(Dy), dz(Dz)
    {}
  };

  template <class T, int Dim, char i, char j, int axis>
  const Tensor2_symmetric_Expr<
    const d_one_sided_Tensor2_symmetric<T, Dim, i, j, axis>,
    typename promote<T, double>::V, Dim, i, j>
  d_one_sided(const Tensor2_symmetric<T *, Dim> &a, const Number<axis> n1,
              const Index<i, Dim> index1, const Index<j, Dim> index2,
              const int &di, const int &dj, const int &dk, const double &dx,
              const double &dy, const double &dz)
  {
    using TensorExpr = d_one_sided_Tensor2_symmetric<T, Dim, i, j, axis>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, double>::V,
                                  Dim, i, j>(
      TensorExpr(a, di, dj, dk, dx, dy, dz));
  }
}
