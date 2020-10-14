/* Adds two Tensor3_antisymmetrics together to make a Dg. */

#pragma once

namespace FTensor
{
  /* A(i,j,k) + B(k,j,i) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Tensor3_antisymmetric_or_Tensor3_antisymmetric
  {
    Tensor3_antisymmetric_Expr<A, T, Dim, Dim, i, j, k> iterA;
    Tensor3_antisymmetric_Expr<B, U, Dim, Dim, k, j, i> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N3, N2) + iterB(N2, N3, N1);
    }

    Tensor3_antisymmetric_or_Tensor3_antisymmetric(
      const Tensor3_antisymmetric_Expr<A, T, Dim, Dim, i, j, k> &a,
      const Tensor3_antisymmetric_Expr<B, U, Dim, Dim, k, j, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Dg_Expr<
    Tensor3_antisymmetric_or_Tensor3_antisymmetric<A, B, T, U, Dim, i, j, k>,
    typename promote<T, U>::V, Dim, Dim, i, k, j>
  operator||(const Tensor3_antisymmetric_Expr<A, T, Dim, Dim, i, j, k> &a,
             const Tensor3_antisymmetric_Expr<B, U, Dim, Dim, k, j, i> &b)
  {
    using TensorExpr
      = Tensor3_antisymmetric_or_Tensor3_antisymmetric<A, B, T, U, Dim, i, j,
                                                       k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim, i, k, j>(
      TensorExpr(a, b));
  }
}
