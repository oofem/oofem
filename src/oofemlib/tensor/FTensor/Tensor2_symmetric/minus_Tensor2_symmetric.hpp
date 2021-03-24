/* Unary minus operator. */

#pragma once

namespace FTensor
{
  template <class A, class T, int Dim, char i, char j>
  class minus_Tensor2_symmetric
  {
  public:
    Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;

  public:
    T operator()(const int N1, const int N2) const { return -iterA(N1, N2); }
    minus_Tensor2_symmetric(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a)
        : iterA(a)
    {}
  };

  template <class A, class T, int Dim, char i, char j>
  Tensor2_symmetric_Expr<minus_Tensor2_symmetric<A, T, Dim, i, j>, T, Dim, i, j>
  operator-(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a)
  {
    using TensorExpr = minus_Tensor2_symmetric<A, T, Dim, i, j>;
    return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(TensorExpr(a));
  }
}
