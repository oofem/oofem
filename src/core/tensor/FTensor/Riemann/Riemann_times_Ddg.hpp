/* This file has the declaration for Riemann*Ddg
   yielding a typename promote<T,U>::V.  I simplify the expression by
   removing the identically zero components of Riemann. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  typename promote<T, U>::V
  operator*(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
            const Ddg_Expr<B, U, Dim, Dim, i, j, k, l> &b)
  {
    return a(1, 0, 0, 1) * b(1, 0, 0, 1) + a(2, 0, 0, 1) * b(2, 0, 0, 1)
           + a(0, 1, 0, 1) * b(0, 1, 0, 1) + a(2, 1, 0, 1) * b(2, 1, 0, 1)
           + a(0, 2, 0, 1) * b(0, 2, 0, 1) + a(1, 2, 0, 1) * b(1, 2, 0, 1)
           + a(1, 0, 0, 2) * b(1, 0, 0, 2) + a(2, 0, 0, 2) * b(2, 0, 0, 2)
           + a(0, 1, 0, 2) * b(0, 1, 0, 2) + a(2, 1, 0, 2) * b(2, 1, 0, 2)
           + a(0, 2, 0, 2) * b(0, 2, 0, 2) + a(1, 2, 0, 2) * b(1, 2, 0, 2)
           + a(1, 0, 1, 0) * b(1, 0, 1, 0) + a(2, 0, 1, 0) * b(2, 0, 1, 0)
           + a(0, 1, 1, 0) * b(0, 1, 1, 0) + a(2, 1, 1, 0) * b(2, 1, 1, 0)
           + a(0, 2, 1, 0) * b(0, 2, 1, 0) + a(1, 2, 1, 0) * b(1, 2, 1, 0)
           + a(1, 0, 1, 2) * b(1, 0, 1, 2) + a(2, 0, 1, 2) * b(2, 0, 1, 2)
           + a(0, 1, 1, 2) * b(0, 1, 1, 2) + a(2, 1, 1, 2) * b(2, 1, 1, 2)
           + a(0, 2, 1, 2) * b(0, 2, 1, 2) + a(1, 2, 1, 2) * b(1, 2, 1, 2)
           + a(1, 0, 2, 0) * b(1, 0, 2, 0) + a(2, 0, 2, 0) * b(2, 0, 2, 0)
           + a(0, 1, 2, 0) * b(0, 1, 2, 0) + a(2, 1, 2, 0) * b(2, 1, 2, 0)
           + a(0, 2, 2, 0) * b(0, 2, 2, 0) + a(1, 2, 2, 0) * b(1, 2, 2, 0)
           + a(1, 0, 2, 1) * b(1, 0, 2, 1) + a(2, 0, 2, 1) * b(2, 0, 2, 1)
           + a(0, 1, 2, 1) * b(0, 1, 2, 1) + a(2, 1, 2, 1) * b(2, 1, 2, 1)
           + a(0, 2, 2, 1) * b(0, 2, 2, 1) + a(1, 2, 2, 1) * b(1, 2, 2, 1);
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  typename promote<T, U>::V
  operator*(const Ddg_Expr<B, U, Dim, Dim, i, j, k, l> &b,
            const Riemann_Expr<A, T, Dim, i, j, k, l> &a)

  {
    return operator*(a, b);
  }

  /* A(i,j,k,l)*B(i,k,j,l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  typename promote<T, U>::V
  operator*(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
            const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b)
  {
    return a(1, 0, 0, 1) * b(1, 0, 0, 1) + a(2, 0, 0, 1) * b(2, 0, 0, 1)
           + a(0, 1, 0, 1) * b(0, 0, 1, 1) + a(2, 1, 0, 1) * b(2, 0, 1, 1)
           + a(0, 2, 0, 1) * b(0, 0, 2, 1) + a(1, 2, 0, 1) * b(1, 0, 2, 1)
           + a(1, 0, 0, 2) * b(1, 0, 0, 2) + a(2, 0, 0, 2) * b(2, 0, 0, 2)
           + a(0, 1, 0, 2) * b(0, 0, 1, 2) + a(2, 1, 0, 2) * b(2, 0, 1, 2)
           + a(0, 2, 0, 2) * b(0, 0, 2, 2) + a(1, 2, 0, 2) * b(1, 0, 2, 2)
           + a(1, 0, 1, 0) * b(1, 1, 0, 0) + a(2, 0, 1, 0) * b(2, 1, 0, 0)
           + a(0, 1, 1, 0) * b(0, 1, 1, 0) + a(2, 1, 1, 0) * b(2, 1, 1, 0)
           + a(0, 2, 1, 0) * b(0, 1, 2, 0) + a(1, 2, 1, 0) * b(1, 1, 2, 0)
           + a(1, 0, 1, 2) * b(1, 1, 0, 2) + a(2, 0, 1, 2) * b(2, 1, 0, 2)
           + a(0, 1, 1, 2) * b(0, 1, 1, 2) + a(2, 1, 1, 2) * b(2, 1, 1, 2)
           + a(0, 2, 1, 2) * b(0, 1, 2, 2) + a(1, 2, 1, 2) * b(1, 1, 2, 2)
           + a(1, 0, 2, 0) * b(1, 2, 0, 0) + a(2, 0, 2, 0) * b(2, 2, 0, 0)
           + a(0, 1, 2, 0) * b(0, 2, 1, 0) + a(2, 1, 2, 0) * b(2, 2, 1, 0)
           + a(0, 2, 2, 0) * b(0, 2, 2, 0) + a(1, 2, 2, 0) * b(1, 2, 2, 0)
           + a(1, 0, 2, 1) * b(1, 2, 0, 1) + a(2, 0, 2, 1) * b(2, 2, 0, 1)
           + a(0, 1, 2, 1) * b(0, 2, 1, 1) + a(2, 1, 2, 1) * b(2, 2, 1, 1)
           + a(0, 2, 2, 1) * b(0, 2, 2, 1) + a(1, 2, 2, 1) * b(1, 2, 2, 1);
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  typename promote<T, U>::V
  operator*(const Ddg_Expr<B, U, Dim, Dim, i, j, k, l> &b,
            const Riemann_Expr<A, T, Dim, i, k, j, l> &a)

  {
    return operator*(a, b);
  }
}
