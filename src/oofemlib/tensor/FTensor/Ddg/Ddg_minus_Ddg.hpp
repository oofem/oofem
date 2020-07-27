/* Subtracts Ddg-Ddg -> Ddg */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l) - B(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_minus_Ddg
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) - iterB(N1, N2, N3, N4);
    }

    Ddg_minus_Ddg(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                  const Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Ddg_Expr<Ddg_minus_Ddg<A, B, T, U, Dim01, Dim23, i, j, k, l>,
           typename promote<T, U>::V, Dim01, Dim23, i, j, k, l>
  operator-(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> &b)
  {
    using TensorExpr = Ddg_minus_Ddg<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j,
                    k, l>(TensorExpr(a, b));
  }
}
