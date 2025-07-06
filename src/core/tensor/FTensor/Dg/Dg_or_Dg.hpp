/* Adds two Dg's to make a Dg, but with different
   symmetries. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)+B(i,k,j) -> Dg(j,k,i) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Dg_or_Dg_12
  {
    Dg_Expr<A, T, Dim, Dim, i, j, k> iterA;
    Dg_Expr<B, U, Dim, Dim, i, k, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N3, N1, N2) + iterB(N3, N2, N1);
    }

    Dg_or_Dg_12(const Dg_Expr<A, T, Dim, Dim, i, j, k> &a,
                const Dg_Expr<B, U, Dim, Dim, i, k, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Dg_Expr<Dg_or_Dg_12<A, B, T, U, Dim, i, j, k>, typename promote<T, U>::V,
          Dim, Dim, j, k, i>
  operator||(const Dg_Expr<A, T, Dim, Dim, i, j, k> &a,
             const Dg_Expr<B, U, Dim, Dim, i, k, j> &b)
  {
    using TensorExpr = Dg_or_Dg_12<A, B, T, U, Dim, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim, j, k, i>(
      TensorExpr(a, b));
  }

  /* A(j,i,k)+B(k,i,j) -> Dg(j,k,i) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Dg_or_Dg_02
  {
    Dg_Expr<A, T, Dim, Dim, j, i, k> iterA;
    Dg_Expr<B, U, Dim, Dim, k, i, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N3, N2) + iterB(N2, N3, N1);
    }

    Dg_or_Dg_02(const Dg_Expr<A, T, Dim, Dim, j, i, k> &a,
                const Dg_Expr<B, U, Dim, Dim, k, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Dg_Expr<Dg_or_Dg_02<A, B, T, U, Dim, i, j, k>, typename promote<T, U>::V,
          Dim, Dim, j, k, i>
  operator||(const Dg_Expr<A, T, Dim, Dim, j, i, k> &a,
             const Dg_Expr<B, U, Dim, Dim, k, i, j> &b)
  {
    using TensorExpr = Dg_or_Dg_02<A, B, T, U, Dim, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim, j, k, i>(
      TensorExpr(a, b));
  }
}
