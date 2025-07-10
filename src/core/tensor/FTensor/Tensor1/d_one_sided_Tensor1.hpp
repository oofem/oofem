/* Takes a one-sided derivative of a Tensor1 in a particular
   direction, yieldind a Tensor1. */

#pragma once

namespace FTensor
{
  template <class T, int Dim, char i, int axis> class d_one_sided_Tensor1
  {
    const Tensor1<T *, Dim> &a;
    const int di, dj, dk;
    const double dx, dy, dz;

  public:
    typename promote<T, double>::V operator()(const int N1) const
    {
      return axis == 0 ? (a(N1) - *(a.ptr(N1) - di)) * dx
                       : (axis == 1 ? (a(N1) - *(a.ptr(N1) - dj)) * dy
                                    : (a(N1) - *(a.ptr(N1) - dk)) * dz);
    }
    d_one_sided_Tensor1(const Tensor1<T *, Dim> &A, const int Di, const int Dj,
                        const int Dk, const double Dx, const double Dy,
                        const double Dz)
        : a(A), di(Di), dj(Dj), dk(Dk), dx(Dx), dy(Dy), dz(Dz)
    {}
  };

  template <class T, int Dim, char i, int axis>
  const Tensor1_Expr<const d_one_sided_Tensor1<T, Dim, i, axis>,
                     typename promote<T, double>::V, Dim, i>
  d_one_sided(const Tensor1<T *, Dim> &a, const Number<axis> n1,
              const Index<i, Dim> index1, const int &di, const int &dj,
              const int &dk, const double &dx, const double &dy,
              const double &dz)
  {
    using TensorExpr = d_one_sided_Tensor1<T, Dim, i, axis>;
    return Tensor1_Expr<TensorExpr, typename promote<T, double>::V, Dim, i>(
      TensorExpr(a, di, dj, dk, dx, dy, dz));
  }
}
