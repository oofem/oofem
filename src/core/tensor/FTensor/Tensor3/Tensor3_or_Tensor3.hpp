/* Adds a Tensor3 to a Tensor3, yielding a Dg. */

#pragma once

namespace FTensor
{
  /* Base Template */
  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim2_0, int Dim0_1, int Dim1_1, int Dim2_1, char i0, char j0,
            char k0, char i1, char j1, char k1>
  class Tensor3_or_Tensor3
  {};

  /* A(i,j,k)+B(i,k,j)->Dg */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  class Tensor3_or_Tensor3<A, B, T, U, Dim0, Dim12, Dim12, Dim0, Dim12, Dim12,
                           i, j, k, i, k, j>
  {
    Tensor3_Expr<A, T, Dim0, Dim12, Dim12, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, k, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N3, N1, N2) + iterB(N3, N2, N1);
    }

    Tensor3_or_Tensor3(
      const Tensor3_Expr<A, T, Dim0, Dim12, Dim12, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, k, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j,k)+B(k,j,i)->Dg */

  template <class A, class B, class T, class U, int Dim02, int Dim1, char i,
            char j, char k>
  class Tensor3_or_Tensor3<A, B, T, U, Dim02, Dim1, Dim02, Dim02, Dim1, Dim02,
                           i, j, k, k, j, i>
  {
    Tensor3_Expr<A, T, Dim02, Dim1, Dim02, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim02, Dim1, Dim02, k, j, i> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N3, N2) + iterB(N2, N3, N1);
    }

    Tensor3_or_Tensor3(
      const Tensor3_Expr<A, T, Dim02, Dim1, Dim02, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim02, Dim1, Dim02, k, j, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim2_0, int Dim0_1, int Dim1_1, int Dim2_1, char i0, char j0,
            char k0, char i1, char j1, char k1>
  auto
  operator||(const Tensor3_Expr<A, T, Dim0_0, Dim1_0, Dim2_0, i0, j0, k0> &a,
             const Tensor3_Expr<B, U, Dim0_1, Dim1_1, Dim2_1, i1, j1, k1> &b)
  {
    using TensorExpr
      = Tensor3_or_Tensor3<A, B, T, U, Dim0_0, Dim1_0, Dim2_0, Dim0_1, Dim1_1,
                           Dim2_1, i0, j0, k0, i1, j1, k1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the || operator");

    // Definition of Helper constexpr variables
    constexpr char i = (i0 == i1) ? j0 : i0, j = (i0 == i1) ? i0 : j0;

    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim0_0, Dim1_0, i,
                   k0, j>(TensorExpr(a, b));
  }
}
