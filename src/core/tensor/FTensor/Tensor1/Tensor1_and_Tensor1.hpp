/* Multiply two Tensor1's together but don't contract, yielding a
   Tensor1. */

#pragma once

namespace FTensor
{
  template <class A, class B, class T, class U, int Dim, char i>
  class Tensor1_and_Tensor1
  {
    Tensor1_Expr<A, T, Dim, i> iterA;
    Tensor1_Expr<B, U, Dim, i> iterB;

  public:
    typename promote<T, U>::V operator()(const int N) const
    {
      return iterA(N) * iterB(N);
    }

    Tensor1_and_Tensor1(const Tensor1_Expr<A, T, Dim, i> &a,
                        const Tensor1_Expr<B, U, Dim, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, char i>
  Tensor1_Expr<Tensor1_and_Tensor1<A, B, T, U, Dim, i>,
               typename promote<T, U>::V, Dim, i>
  operator&(const Tensor1_Expr<A, T, Dim, i> &a,
            const Tensor1_Expr<B, U, Dim, i> &b)
  {
    using TensorExpr = Tensor1_and_Tensor1<A, B, T, U, Dim, i>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, b));
  }
}
