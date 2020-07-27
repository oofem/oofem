/* Multiply a Tensor2_symmetric and a Dg together but don't
   contract, yielding a Dg. */

#pragma once

namespace FTensor
{
  /* A(i,j,k) & B(i,j) -> Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  class Dg_and_Tensor2_symmetric
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Tensor2_symmetric_Expr<B, U, Dim01, i, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) * iterB(N1, N2);
    }

    Dg_and_Tensor2_symmetric(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
                             const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Dg_Expr<Dg_and_Tensor2_symmetric<A, B, T, U, Dim01, Dim2, i, j, k>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, k>
  operator&(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &b)
  {
    using TensorExpr
      = Dg_and_Tensor2_symmetric<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, k>(
      TensorExpr(a, b));
  }

  /* B(i,j) & A(i,j,k) -> Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Dg_Expr<Dg_and_Tensor2_symmetric<A, B, T, U, Dim01, Dim2, i, j, k>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, k>
  operator&(const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &b,
            const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a)
  {
    using TensorExpr
      = Dg_and_Tensor2_symmetric<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, k>(
      TensorExpr(a, b));
  }
}
