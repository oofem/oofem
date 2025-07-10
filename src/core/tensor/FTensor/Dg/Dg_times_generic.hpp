/* Multiplies a Dg with a generic, yielding a Dg. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim01, int Dim2, char i, char j,
            char k>
  class Dg_times_generic
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    U d;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) * d;
    }

    Dg_times_generic(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a, const U &d0)
        : iterA(a), d(d0)
    {}
  };

  /* A(i,j,k)*d0->Dg */

  template <class A, class T, class U, int Dim01, int Dim2, char i, char j,
            char k>
  Dg_Expr<Dg_times_generic<A, T, U, Dim01, Dim2, i, j, k>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, k>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a, const U &d0)
  {
    using TensorExpr = Dg_times_generic<A, T, U, Dim01, Dim2, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, k>(
      TensorExpr(a, d0));
  }

  /* d0*A(i,j,k)->Dg */

  template <class A, class T, class U, int Dim01, int Dim2, char i, char j,
            char k>
  Dg_Expr<Dg_times_generic<A, T, U, Dim01, Dim2, i, j, k>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, k>
  operator*(const U &d0, const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a)
  {
    using TensorExpr = Dg_times_generic<A, T, U, Dim01, Dim2, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, k>(
      TensorExpr(a, d0));
  }
}
