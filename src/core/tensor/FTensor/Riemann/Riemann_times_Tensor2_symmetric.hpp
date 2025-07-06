/* This file has all of the declarations for expressions like
   Riemann*Tensor2_symmetric and
   Tensor2_symmetric*Riemann, yielding a Tensor2_symmetric. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(i,k) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  class Riemann_times_Tensor2_symmetric_0
  {
    Riemann_Expr<A, T, Dim, i, j, k, l> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, i, k> iterB;

  public:
    Riemann_times_Tensor2_symmetric_0(
      const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, i, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(0, N1, 0, N2) * iterB(0, 0)
             + iterA(0, N1, 1, N2) * iterB(0, 1)
             + iterA(0, N1, 2, N2) * iterB(0, 2)
             + iterA(1, N1, 0, N2) * iterB(1, 0)
             + iterA(1, N1, 1, N2) * iterB(1, 1)
             + iterA(1, N1, 2, N2) * iterB(1, 2)
             + iterA(2, N1, 0, N2) * iterB(2, 0)
             + iterA(2, N1, 1, N2) * iterB(2, 1)
             + iterA(2, N1, 2, N2) * iterB(2, 2);
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Tensor2_symmetric_Expr<
    Riemann_times_Tensor2_symmetric_0<A, B, T, U, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim, j, l>
  operator*(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, i, k> &b)
  {
    using TensorExpr
      = Riemann_times_Tensor2_symmetric_0<A, B, T, U, Dim, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim,
                                  j, l>(TensorExpr(a, b));
  }

  /* B(i,k)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Tensor2_symmetric_Expr<
    Riemann_times_Tensor2_symmetric_0<A, B, T, U, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim, j, l>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, i, k> &b,
            const Riemann_Expr<A, T, Dim, i, j, k, l> &a)
  {
    using TensorExpr
      = Riemann_times_Tensor2_symmetric_0<A, B, T, U, Dim, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim,
                                  j, l>(TensorExpr(a, b));
  }
}
