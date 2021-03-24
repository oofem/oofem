/* Adds a Tensor2_symmetric to a generic, yielding a
   Tensor2_symmetric. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim, char i, char j>
  auto operator+(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a, const U &d0)
  {
    auto TensorExpr = [&a, d0](const int N1, const int N2) {
      return a.operator()(N1, N2) + d0;
    };
    return Tensor2_symmetric_Expr<decltype(TensorExpr),
                                  typename promote<T, U>::V, Dim, i, j>(
      TensorExpr);
  }

  template <class A, class T, class U, int Dim, char i, char j>
  auto operator+(const U &d0, const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a)
  {
    auto TensorExpr = [&a, d0](const int N1, const int N2) {
      return a.operator()(N1, N2) + d0;
    };
    return Tensor2_symmetric_Expr<decltype(TensorExpr),
                                  typename promote<T, U>::V, Dim, i, j>(
      TensorExpr);
  }
}
