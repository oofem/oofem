/* Declares a wrapper class for the unary minus (-) operator. */

#pragma once

namespace FTensor
{
  template <class A, class T, int Dim, char i> class minus_Tensor1
  {
    Tensor1_Expr<A, T, Dim, i> iterA;

  public:
    T operator()(const int N) const { return -iterA(N); }

    minus_Tensor1(const Tensor1_Expr<A, T, Dim, i> &a) : iterA(a) {}
  };

  template <class A, class T, int Dim, char i>
  Tensor1_Expr<minus_Tensor1<A, T, Dim, i>, T, Dim, i>
  operator-(const Tensor1_Expr<A, T, Dim, i> &a)
  {
    using TensorExpr = minus_Tensor1<A, T, Dim, i>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(a));
  }
}
