/* Subtracts a Tensor2 from a Tensor2_symmetric (or vice versa),
   yielding a Tensor2. */

#pragma once

namespace FTensor
{
  /* Base Template */
  template <class A, class B, class T, class U, int Dim_0, int Dim0_1,
            int Dim1_1, char i0, char j0, char i1, char j1>
  class Tensor2_symmetric_minus_Tensor2
  {};

  /* A(i,j)-B(i,j), A is symmetric, B is not. */

  template <class A, class B, class T, class U, int Dim, char i, char j>
  class Tensor2_symmetric_minus_Tensor2<A, B, T, U, Dim, Dim, Dim, i, j, i, j>
  {
    Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    Tensor2_Expr<B, U, Dim, Dim, i, j> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) - iterB(N1, N2);
    }

    Tensor2_symmetric_minus_Tensor2(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
      const Tensor2_Expr<B, U, Dim, Dim, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j)-B(j,i), A is symmetric, B is not. */

  template <class A, class B, class T, class U, int Dim, char i, char j>
  class Tensor2_symmetric_minus_Tensor2<A, B, T, U, Dim, Dim, Dim, i, j, j, i>
  {
    Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    Tensor2_Expr<B, U, Dim, Dim, j, i> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) - iterB(N2, N1);
    }

    Tensor2_symmetric_minus_Tensor2(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
      const Tensor2_Expr<B, U, Dim, Dim, j, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim_0, int Dim0_1,
            int Dim1_1, char i0, char j0, char i1, char j1>
  Tensor2_Expr<Tensor2_symmetric_minus_Tensor2<A, B, T, U, Dim_0, Dim0_1,
                                               Dim1_1, i0, j0, i1, j1>,
               typename promote<T, U>::V, Dim_0, Dim_0, i0, j0>
  operator-(const Tensor2_symmetric_Expr<A, T, Dim_0, i0, j0> &a,
            const Tensor2_Expr<B, U, Dim0_1, Dim1_1, i1, j1> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_minus_Tensor2<A, B, T, U, Dim_0, Dim0_1, Dim1_1, i0,
                                        j0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the - operator");
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim_0, Dim_0,
                        i0, j0>(TensorExpr(a, b));
  }

  // TODO: Maybe we should consider moving the following functions to the space
  // of Tensor2 for organization sake
  /* B(i,j)-A(i,j), A is symmetric, B is not. */

  /* Base Template */
  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim_1, char i0, char j0, char i1, char j1>
  class Tensor2_minus_Tensor2_symmetric
  {};

  template <class A, class B, class T, class U, int Dim, char i, char j>
  class Tensor2_minus_Tensor2_symmetric<A, B, T, U, Dim, Dim, Dim, i, j, i, j>
  {
    Tensor2_Expr<A, T, Dim, Dim, i, j> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, i, j> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) - iterB(N1, N2);
    }

    Tensor2_minus_Tensor2_symmetric(
      const Tensor2_Expr<A, T, Dim, Dim, i, j> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* B(i,j)-A(j,i), A is symmetric, B is not. */

  template <class A, class B, class T, class U, int Dim, char i, char j>
  class Tensor2_minus_Tensor2_symmetric<A, B, T, U, Dim, Dim, Dim, i, j, j, i>
  {
    Tensor2_Expr<A, T, Dim, Dim, i, j> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, j, i> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) - iterB(N2, N1);
    }

    Tensor2_minus_Tensor2_symmetric(
      const Tensor2_Expr<A, T, Dim, Dim, i, j> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim_1, char i0, char j0, char i1, char j1>
  Tensor2_Expr<Tensor2_minus_Tensor2_symmetric<A, B, T, U, Dim0_0, Dim1_0,
                                               Dim_1, i0, j0, i1, j1>,
               typename promote<T, U>::V, Dim0_0, Dim1_0, i0, j0>
  operator-(const Tensor2_Expr<A, T, Dim0_0, Dim1_0, i0, j0> &a,
            const Tensor2_symmetric_Expr<B, U, Dim_1, i1, j1> &b)
  {
    using TensorExpr
      = Tensor2_minus_Tensor2_symmetric<A, B, T, U, Dim0_0, Dim1_0, Dim_1, i0,
                                        j0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the - operator");
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0_0, Dim1_0,
                        i0, j0>(TensorExpr(a, b));
  }
}
