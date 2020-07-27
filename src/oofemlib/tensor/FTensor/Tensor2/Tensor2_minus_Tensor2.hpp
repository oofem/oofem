/* Subtracts a Tensor2 from a Tensor2, yielding a Tensor2. */

#pragma once

namespace FTensor
{
  /* Base template */
  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim0_1, int Dim1_1, char i0, char j0, char i1, char j1>
  class Tensor2_minus_Tensor2
  {};

  /* A(i,j)-B(i,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  class Tensor2_minus_Tensor2<A, B, T, U, Dim0, Dim1, Dim0, Dim1, i, j, i, j>
  {
    const Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;
    const Tensor2_Expr<B, U, Dim0, Dim1, i, j> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) - iterB(N1, N2);
    }

    Tensor2_minus_Tensor2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                          const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j)-B(j,i) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  class Tensor2_minus_Tensor2<A, B, T, U, Dim0, Dim1, Dim1, Dim0, i, j, j, i>
  {
    const Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;
    const Tensor2_Expr<B, U, Dim1, Dim0, j, i> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) - iterB(N2, N1);
    }

    Tensor2_minus_Tensor2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                          const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim0_1, int Dim1_1, char i0, char j0, char i1, char j1>
  Tensor2_Expr<Tensor2_minus_Tensor2<A, B, T, U, Dim0_0, Dim1_0, Dim0_1,
                                     Dim1_1, i0, j0, i1, j1>,
               typename promote<T, U>::V, Dim0_0, Dim1_0, i0, j0>
  operator-(const Tensor2_Expr<A, T, Dim0_0, Dim1_0, i0, j0> &a,
            const Tensor2_Expr<B, U, Dim0_1, Dim1_1, i1, j1> &b)
  {
    using TensorExpr = Tensor2_minus_Tensor2<A, B, T, U, Dim0_0, Dim1_0,
                                             Dim0_1, Dim1_1, i0, j0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the - operator");
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0_0, Dim1_0,
                        i0, j0>(TensorExpr(a, b));
  }
}
