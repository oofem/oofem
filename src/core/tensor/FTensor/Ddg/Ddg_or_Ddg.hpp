/* Declares expressions of Ddg || Ddg.  This adds them in
   a different way, but still ending up with a Ddg. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)+B(i,l,k,j) -> Ddg */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  class Ddg_or_Ddg
  {
    Ddg_Expr<A, T, Dim, Dim, i, j, k, l> iterA;
    Ddg_Expr<B, U, Dim, Dim, i, l, k, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N3, N2, N4) + iterB(N1, N4, N2, N3);
    }
    Ddg_or_Ddg(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
               const Ddg_Expr<B, U, Dim, Dim, i, l, k, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Ddg_Expr<Ddg_or_Ddg<A, B, T, U, Dim, i, j, k, l>, typename promote<T, U>::V,
           Dim, Dim, i, k, j, l>
  operator||(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
             const Ddg_Expr<B, U, Dim, Dim, i, l, k, j> &b)
  {
    using TensorExpr = Ddg_or_Ddg<A, B, T, U, Dim, i, j, k, l>;
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim, i, k, j,
                    l>(TensorExpr(a, b));
  }
}
