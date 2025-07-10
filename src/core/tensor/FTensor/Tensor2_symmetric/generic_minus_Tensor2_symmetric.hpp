/* Subtracts a generic from a Tensor2_symmetric, yielding a
   Tensor2_symmetric. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim, char i, char j>
  class generic_minus_Tensor2_symmetric
  {
    Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    U d;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return d - iterA(N1, N2);
    }
    generic_minus_Tensor2_symmetric(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a, const U &d0)
        : iterA(a), d(d0)
    {}
  };

  template <class A, class T, class U, int Dim, char i, char j>
  Tensor2_symmetric_Expr<generic_minus_Tensor2_symmetric<A, T, U, Dim, i, j>,
                         typename promote<T, U>::V, Dim, i, j>
  operator-(const U &d0, const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a)
  {
    using TensorExpr = generic_minus_Tensor2_symmetric<A, T, U, Dim, i, j>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim,
                                  i, j>(TensorExpr(a, d0));
  }
}
