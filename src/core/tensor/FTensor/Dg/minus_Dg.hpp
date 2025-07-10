/* Unary minus operator. */

#pragma once

namespace FTensor
{
  template <class A, class T, int Dim01, int Dim2, char i, char j, char k>
  class minus_Dg
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;

  public:
    T operator()(const int N1, const int N2, const int N3) const
    {
      return -iterA(N1, N2, N3);
    }

    minus_Dg(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a) : iterA(a) {}
  };

  template <class A, class T, int Dim01, int Dim2, char i, char j, char k>
  Dg_Expr<minus_Dg<A, T, Dim01, Dim2, i, j, k>, T, Dim01, Dim2, i, j, k>
  operator-(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a)
  {
    using TensorExpr = minus_Dg<A, T, Dim01, Dim2, i, j, k>;
    return Dg_Expr<TensorExpr, T, Dim01, Dim2, i, j, k>(TensorExpr(a));
  }
}
