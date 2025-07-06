/* Multiplies a Tensor2 with a generic, yielding a Tensor2. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim0, int Dim1, char i, char j>
  class Tensor2_times_generic
  {
    const Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;
    const U d;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) * d;
    }

    Tensor2_times_generic(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                          const U &d0)
        : iterA(a), d(d0)
    {}
  };

  template <class A, class T, class U, int Dim0, int Dim1, char i, char j>
  Tensor2_Expr<Tensor2_times_generic<A, T, U, Dim0, Dim1, i, j>,
               typename promote<T, U>::V, Dim0, Dim1, i, j>
  operator*(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a, const U &d0)
  {
    using TensorExpr = Tensor2_times_generic<A, T, U, Dim0, Dim1, i, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, d0));
  }

  template <class A, class T, class U, int Dim0, int Dim1, char i, char j>
  Tensor2_Expr<Tensor2_times_generic<A, T, U, Dim0, Dim1, i, j>,
               typename promote<T, U>::V, Dim0, Dim1, i, j>
  operator*(const U &d0, const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a)
  {
    using TensorExpr = Tensor2_times_generic<A, T, U, Dim0, Dim1, i, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, d0));
  }
}
