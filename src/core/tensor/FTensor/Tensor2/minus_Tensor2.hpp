/* Unary minus operator. */

#pragma once

namespace FTensor
{
  template <class A, class T, int Dim0, int Dim1, char i, char j>
  class minus_Tensor2
  {
    const Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;

  public:
    T operator()(const int N1, const int N2) const { return -iterA(N1, N2); }

    minus_Tensor2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a) : iterA(a) {}
  };

  template <class A, class T, int Dim0, int Dim1, char i, char j>
  Tensor2_Expr<minus_Tensor2<A, T, Dim0, Dim1, i, j>, T, Dim0, Dim1, i, j>
  operator-(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a)
  {
    using TensorExpr = minus_Tensor2<A, T, Dim0, Dim1, i, j>;
    return Tensor2_Expr<TensorExpr, T, Dim0, Dim1, i, j>(TensorExpr(a));
  }
}
