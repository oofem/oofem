/* Unary minus operator. */

#pragma once

namespace FTensor
{
  template <class A, class T, int Dim01, int Dim23, char i, char j, char k,
            char l>
  class minus_Ddg
  {
  public:
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;

  public:
    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return -iterA(N1, N2, N3, N4);
    }
    minus_Ddg(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a) : iterA(a) {}
  };

  template <class A, class T, int Dim01, int Dim23, char i, char j, char k,
            char l>
  Ddg_Expr<minus_Ddg<A, T, Dim01, Dim23, i, j, k, l>, T, Dim01, Dim23, i, j, k,
           l>
  operator-(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    using TensorExpr = minus_Ddg<A, T, Dim01, Dim23, i, j, k, l>;
    return Ddg_Expr<TensorExpr, T, Dim01, Dim23, i, j, k, l>(TensorExpr(a));
  }
}
