/* Subtracts a Riemann from a Riemann, yielding a
   Riemann. */

#pragma once

namespace FTensor
{
  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  class Riemann_minus_Riemann
  {
    Riemann_Expr<A, T, Dim, i, j, k, l> iterA;
    Riemann_Expr<B, U, Dim, i, j, k, l> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) - iterB(N1, N2, N3, N4);
    }

    Riemann_minus_Riemann(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
                          const Riemann_Expr<B, U, Dim, i, j, k, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Riemann_Expr<Riemann_minus_Riemann<A, B, T, U, Dim, i, j, k, l>,
               typename promote<T, U>::V, Dim, i, j, k, l>
  operator-(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
            const Riemann_Expr<B, U, Dim, i, j, k, l> &b)
  {
    using TensorExpr = Riemann_minus_Riemann<A, B, T, U, Dim, i, j, k, l>;
    return Riemann_Expr<TensorExpr, typename promote<T, U>::V, Dim, i, j, k, l>(
      TensorExpr(a, b));
  }
}
