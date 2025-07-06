/* complex conjugate operator. */

#pragma once

namespace FTensor
{
  template <class A, class T, int Dim0, int Dim1, char i, char j>
  class conj_Tensor2
  {
    const Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;

  public:
    T operator()(const int N1, const int N2) const
    {
      return conj(iterA(N1, N2));
    }

    conj_Tensor2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a) : iterA(a) {}
  };

  template <class A, class T, int Dim0, int Dim1, char i, char j>
  Tensor2_Expr<conj_Tensor2<A, T, Dim0, Dim1, i, j>, T, Dim0, Dim1, i, j>
  conj(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a)
  {
    using TensorExpr = conj_Tensor2<A, T, Dim0, Dim1, i, j>;
    return Tensor2_Expr<TensorExpr, T, Dim0, Dim1, i, j>(TensorExpr(a));
  }
}
