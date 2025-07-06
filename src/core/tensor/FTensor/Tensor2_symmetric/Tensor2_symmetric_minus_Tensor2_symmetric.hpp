/* Subtracts a Tensor2_symmetric from a Tensor2_symmetric, yielding a
   Tensor2_symmetric. */

#pragma once

namespace FTensor
{
  /* Base Template */
  template <class A, class B, class T, class U, int Dim_0, int Dim_1, char i0,
            char j0, char i1, char j1>
  class Tensor2_symmetric_minus_Tensor2_symmetric
  {};

  /* A(i,j) - B(i,j) */

  template <class A, class B, class T, class U, int Dim, char i, char j>
  class Tensor2_symmetric_minus_Tensor2_symmetric<A, B, T, U, Dim, Dim, i, j,
                                                  i, j>
  {
    Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, i, j> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) - iterB(N1, N2);
    }

    Tensor2_symmetric_minus_Tensor2_symmetric(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j) - B(j,i) */

  template <class A, class B, class T, class U, int Dim, char i, char j>
  class Tensor2_symmetric_minus_Tensor2_symmetric<A, B, T, U, Dim, Dim, i, j,
                                                  j, i>
  {
    Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, j, i> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) - iterB(N1, N2);
    }

    Tensor2_symmetric_minus_Tensor2_symmetric(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim_0, int Dim_1, char i0,
            char j0, char i1, char j1>
  Tensor2_symmetric_Expr<Tensor2_symmetric_minus_Tensor2_symmetric<
                           A, B, T, U, Dim_0, Dim_1, i0, j0, i1, j1>,
                         typename promote<T, U>::V, Dim_0, i0, j0>
  operator-(const Tensor2_symmetric_Expr<A, T, Dim_0, i0, j0> &a,
            const Tensor2_symmetric_Expr<B, U, Dim_1, i1, j1> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_minus_Tensor2_symmetric<A, B, T, U, Dim_0, Dim_1, i0,
                                                  j0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the - operator");
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim_0,
                                  i0, j0>(TensorExpr(a, b));
  }
}
