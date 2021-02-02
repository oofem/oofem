/* Multiplies a Tensor3_antisymmetric and a Tensor1, yielding a
   Riemann. */

#pragma once

namespace FTensor
{
  /* A(j,k,l)*B(i) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  class Tensor3_antisymmetric_mod_Tensor1
  {
    Tensor3_antisymmetric_Expr<A, T, Dim, Dim, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim, i> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N2, N3, N4) * iterB(N1);
    }

    Tensor3_antisymmetric_mod_Tensor1(
      const Tensor3_antisymmetric_Expr<A, T, Dim, Dim, j, k, l> &a,
      const Tensor1_Expr<B, U, Dim, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Riemann_Expr<Tensor3_antisymmetric_mod_Tensor1<A, B, T, U, Dim, i, j, k, l>,
               typename promote<T, U>::V, Dim, i, j, k, l>
  operator%(const Tensor3_antisymmetric_Expr<A, T, Dim, Dim, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim, i> &b)
  {
    using TensorExpr
      = Tensor3_antisymmetric_mod_Tensor1<A, B, T, U, Dim, i, j, k, l>;
    return Riemann_Expr<TensorExpr, typename promote<T, U>::V, Dim, i, j, k, l>(
      TensorExpr(a, b));
  }

  /* B(i)*A(j,k,l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Riemann_Expr<Tensor3_antisymmetric_mod_Tensor1<A, B, T, U, Dim, i, j, k, l>,
               typename promote<T, U>::V, Dim, i, j, k, l>
  operator%(const Tensor1_Expr<B, U, Dim, i> &b,
            const Tensor3_antisymmetric_Expr<A, T, Dim, Dim, j, k, l> &a)
  {
    using TensorExpr
      = Tensor3_antisymmetric_mod_Tensor1<A, B, T, U, Dim, i, j, k, l>;
    return Riemann_Expr<TensorExpr, typename promote<T, U>::V, Dim, i, j, k, l>(
      TensorExpr(a, b));
  }
}
