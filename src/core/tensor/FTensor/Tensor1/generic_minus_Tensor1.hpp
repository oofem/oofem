/* Subtracts a Tensor1 from a generic, yielding a Tensor1. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim, char i>
  class generic_minus_Tensor1
  {
    Tensor1_Expr<A, T, Dim, i> iterA;
    U d;

  public:
    typename promote<T, U>::V operator()(const int N) const
    {
      return d - iterA(N);
    }

    generic_minus_Tensor1(const Tensor1_Expr<A, T, Dim, i> &a, const U &d0)
        : iterA(a), d(d0)
    {}
  };

  template <class A, class T, class U, int Dim, char i>
  Tensor1_Expr<generic_minus_Tensor1<A, T, U, Dim, i>,
               typename promote<T, U>::V, Dim, i>
  operator-(const U &d0, const Tensor1_Expr<A, T, Dim, i> &a)
  {
    using TensorExpr = generic_minus_Tensor1<A, T, U, Dim, i>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, d0));
  }
}
