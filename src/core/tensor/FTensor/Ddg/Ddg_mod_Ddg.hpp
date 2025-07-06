/* Subtracts Ddg-Ddg -> Riemann */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l) - B(i,l,k,j) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  class Ddg_mod_Ddg
  {
    Ddg_Expr<A, T, Dim, Dim, i, j, k, l> iterA;
    Ddg_Expr<B, U, Dim, Dim, i, l, k, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) - iterB(N1, N4, N3, N2);
    }

    Ddg_mod_Ddg(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                const Ddg_Expr<B, U, Dim, Dim, i, l, k, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Riemann_Expr<Ddg_mod_Ddg<A, B, T, U, Dim, i, j, k, l>,
               typename promote<T, U>::V, Dim, i, j, k, l>
  operator%(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
            const Ddg_Expr<B, U, Dim, Dim, i, l, k, j> &b)
  {
    using TensorExpr = Ddg_mod_Ddg<A, B, T, U, Dim, i, j, k, l>;
    return Riemann_Expr<TensorExpr, typename promote<T, U>::V, Dim, i, j, k, l>(
      TensorExpr(a, b));
  }
}
