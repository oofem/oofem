/* Multiply a Tensor1 and a Tensor2 together but don't contract, yielding a
   Tensor2. */

/* A(i,j) & B(i) -> Tensor2 */

#pragma once

namespace FTensor
{
  // Base Template
  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim0_1, char i0, char j0, char i1>
  class Tensor2_and_Tensor1
  {};

  /* A(i,j) & B(i) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  class Tensor2_and_Tensor1<A, B, T, U, Dim0, Dim1, Dim0, i, j, i>
  {
    const Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;
    const Tensor1_Expr<B, U, Dim0, i> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) * iterB(N1);
    }

    Tensor2_and_Tensor1(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                        const Tensor1_Expr<B, U, Dim0, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j) & B(j) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  class Tensor2_and_Tensor1<A, B, T, U, Dim0, Dim1, Dim1, i, j, j>
  {
    const Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;
    const Tensor1_Expr<B, U, Dim1, j> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return iterA(N1, N2) * iterB(N2);
    }

    Tensor2_and_Tensor1(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                        const Tensor1_Expr<B, U, Dim1, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j) & B(i/j) -> Tensor2 */

  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim0_1, char i0, char j0, char i1>
  auto operator&(const Tensor2_Expr<A, T, Dim0_0, Dim1_0, i0, j0> &a,
                 const Tensor1_Expr<B, U, Dim0_1, i1> &b)
  {
    using TensorExpr
      = Tensor2_and_Tensor1<A, B, T, U, Dim0_0, Dim1_0, Dim0_1, i0, j0, i1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the & operator");
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0_0, Dim1_0,
                        i0, j0>(TensorExpr(a, b));
  }

  /* B(i/j) & A(i,j) -> Tensor2 */

  // TODO=> We are not respecting operation order, in the really odd case
  // someone uses this operation with non
  // TODO=>  commutable T or U. Any ideas on how should we handel this.

  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim0_1, char i0, char j0, char i1>
  auto operator&(const Tensor1_Expr<B, U, Dim0_1, i1> &b,
                 const Tensor2_Expr<A, T, Dim0_0, Dim1_0, i0, j0> &a)
  {
    using TensorExpr
      = Tensor2_and_Tensor1<A, B, T, U, Dim0_0, Dim1_0, Dim0_1, i0, j0, i1>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0_0, Dim1_0,
                        i0, j0>(TensorExpr(a, b));
  }
}
