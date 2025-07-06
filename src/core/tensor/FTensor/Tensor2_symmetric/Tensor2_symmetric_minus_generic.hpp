/* Subtracts a generic from a Tensor2_symmetric, yielding a
   Tensor2_symmetric. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim, char i, char j>
  class Tensor2_symmetric_minus_generic
  {
    const Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    const U d;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) - d;
    }

    Tensor2_symmetric_minus_generic(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a, const U &d0)
        : iterA(a), d(d0)
    {}
  };

  template <class A, class T, class U, int Dim, char i, char j>
  Tensor2_symmetric_Expr<Tensor2_symmetric_minus_generic<A, T, U, Dim, i, j>,
                         typename promote<T, U>::V, Dim, i, j>
  operator-(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a, const U &d0)
  {
    using TensorExpr = Tensor2_symmetric_minus_generic<A, T, U, Dim, i, j>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim,
                                  i, j>(TensorExpr(a, d0));
  }
}
