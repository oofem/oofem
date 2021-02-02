/* Multiplies a Tensor3 with a generic, yielding a Tensor3. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)*generic */

  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, char i,
            char j, char k>
  auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                 const U &d0)
  {
    auto TensorExpr = [&a, &d0](const int N1, const int N2, const int N3) {
      return a.operator()(N1, N2, N3) * d0;
    };
    return Tensor3_Expr<decltype(TensorExpr), typename promote<T, U>::V, Dim0,
                        Dim1, Dim2, i, j, k>(TensorExpr);
  }

  /* generic*A(i,j,k) */

  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, char i,
            char j, char k>
  auto operator*(const U &d0,
                 const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    auto TensorExpr = [&a, &d0](const int N1, const int N2, const int N3) {
      return d0 * a.operator()(N1, N2, N3);
    };
    return Tensor3_Expr<decltype(TensorExpr), typename promote<T, U>::V, Dim0,
                        Dim1, Dim2, i, j, k>(TensorExpr);
  }
}
