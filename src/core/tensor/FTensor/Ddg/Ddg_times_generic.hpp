/* Multiplies a Ddg with a generic, yielding a
   Ddg. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim01, int Dim23, char i, char j,
            char k, char l>
  auto
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a, const U &d0)
  {
    auto TensorExpr
      = [&a, &d0](const int N1, const int N2, const int N3, const int N4) {
          return a(N1, N2, N3, N4) * d0;
        };
    return Ddg_Expr<decltype(TensorExpr), typename promote<T, U>::V, Dim01,
                    Dim23, i, j, k, l>(TensorExpr);
  }

  template <class A, class T, class U, int Dim01, int Dim23, char i, char j,
            char k, char l>
  auto
  operator*(const U &d0, const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    auto TensorExpr
      = [&a, &d0](const int N1, const int N2, const int N3, const int N4) {
          return d0 * a(N1, N2, N3, N4);
        };
    return Ddg_Expr<decltype(TensorExpr), typename promote<T, U>::V, Dim01,
                    Dim23, i, j, k, l>(TensorExpr);
  }
}
