/* Multiply a Tensor1 and a Dg together but don't contract, yielding a
   Dg. */

#pragma once

namespace FTensor
{
  /* A(i,j,k) & B(k) -> Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  class Dg_and_Tensor1
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Tensor1_Expr<B, U, Dim2, k> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) * iterB(N3);
    }

    Dg_and_Tensor1(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
                   const Tensor1_Expr<B, U, Dim2, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Dg_Expr<Dg_and_Tensor1<A, B, T, U, Dim01, Dim2, i, j, k>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, k>
  operator&(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Tensor1_Expr<B, U, Dim2, k> &b)
  {
    using TensorExpr = Dg_and_Tensor1<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, k>(
      TensorExpr(a, b));
  }

  /* B(k) & A(i,j,k) -> Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Dg_Expr<Dg_and_Tensor1<A, B, T, U, Dim01, Dim2, i, j, k>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, k>
  operator&(const Tensor1_Expr<B, U, Dim2, k> &b,
            const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a)
  {
    using TensorExpr = Dg_and_Tensor1<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, k>(
      TensorExpr(a, b));
  }
}
