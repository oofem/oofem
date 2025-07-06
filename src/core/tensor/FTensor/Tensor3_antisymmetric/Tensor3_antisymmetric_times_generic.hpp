/* Multiplies a Tensor3_antisymmetric with a generic, yielding a
   Tensor3_antisymmetric. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)*generic */

  template <class A, class T, class U, int Dim0, int Dim12, char i, char j,
            char k>
  class Tensor3_antisymmetric_times_generic
  {
    Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> iterA;
    U d;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) * d;
    }

    Tensor3_antisymmetric_times_generic(
      const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
      const U &d0)
        : iterA(a), d(d0)
    {}
  };

  template <class A, class T, class U, int Dim0, int Dim12, char i, char j,
            char k>
  Tensor3_antisymmetric_Expr<
    Tensor3_antisymmetric_times_generic<A, T, U, Dim0, Dim12, i, j, k>,
    typename promote<T, U>::V, Dim0, Dim12, i, j, k>
  operator*(const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
            const U &d0)
  {
    using TensorExpr
      = Tensor3_antisymmetric_times_generic<A, T, U, Dim0, Dim12, i, j, k>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim0, Dim12, i, j, k>(TensorExpr(a, d0));
  }

  /* generic*A(i,j,k) */

  template <class A, class T, class U, int Dim0, int Dim12, char i, char j,
            char k>
  Tensor3_antisymmetric_Expr<
    Tensor3_antisymmetric_times_generic<A, T, U, Dim0, Dim12, i, j, k>,
    typename promote<T, U>::V, Dim0, Dim12, i, j, k>
  operator*(const U &d0,
            const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a)
  {
    using TensorExpr
      = Tensor3_antisymmetric_times_generic<A, T, U, Dim0, Dim12, i, j, k>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim0, Dim12, i, j, k>(TensorExpr(a, d0));
  }
}
