/* Divide a Ddg by a Tensor2_symmetric without contracting, yielding a
   Ddg. */

#pragma once

namespace FTensor
{
  /* Base Template */
  template <class A, class B, class T, class U, int Dim01_0, int Dim23_0,
            int Dim_1, char i0, char j0, char k0, char l0, char i1, char j1>
  class Ddg_mod_Tensor2_symmetric
  {};

  /* A(i,j,k,l) % B(i,j) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_mod_Tensor2_symmetric<A, B, T, U, Dim01, Dim23, Dim01, i, j, k, l,
                                  i, j>
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor2_symmetric_Expr<B, U, Dim01, i, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) / iterB(N1, N2);
    }

    Ddg_mod_Tensor2_symmetric(
      const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* B(i,j) % A(i,j,k,l) -> Ddg */

  /* A(i,j,k,l) % B(k,l) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_mod_Tensor2_symmetric<A, B, T, U, Dim01, Dim23, Dim23, i, j, k, l,
                                  k, l>
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor2_symmetric_Expr<B, U, Dim23, k, l> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) / iterB(N3, N4);
    }

    Ddg_mod_Tensor2_symmetric(
      const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim23, k, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim01_0, int Dim23_0,
            int Dim_1, char i0, char j0, char k0, char l0, char i1, char j1>
  Ddg_Expr<Ddg_mod_Tensor2_symmetric<A, B, T, U, Dim01_0, Dim23_0, Dim_1, i0,
                                     j0, k0, l0, i1, j1>,
           typename promote<T, U>::V, Dim01_0, Dim23_0, i0, j0, k0, l0>
  operator%(const Ddg_Expr<A, T, Dim01_0, Dim23_0, i0, j0, k0, l0> &a,
            const Tensor2_symmetric_Expr<B, U, Dim_1, i1, j1> &b)
  {
    using TensorExpr
      = Ddg_mod_Tensor2_symmetric<A, B, T, U, Dim01_0, Dim23_0, Dim_1, i0, j0,
                                  k0, l0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the % operator");
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01_0, Dim23_0,
                    i0, j0, k0, l0>(TensorExpr(a, b));
  }

  /* B(k,l) % A(i,j,k,l) -> Ddg */

  template <class A, class B, class T, class U, int Dim01_0, int Dim23_0,
            int Dim_1, char i0, char j0, char k0, char l0, char i1, char j1>
  Ddg_Expr<Ddg_mod_Tensor2_symmetric<A, B, T, U, Dim01_0, Dim23_0, Dim_1, i0,
                                     j0, k0, l0, i1, j1>,
           typename promote<T, U>::V, Dim01_0, Dim23_0, i0, j0, k0, l0>
  operator%(const Tensor2_symmetric_Expr<B, U, Dim_1, i1, j1> &b,
            const Ddg_Expr<A, T, Dim01_0, Dim23_0, i0, j0, k0, l0> &a)
  {
    using TensorExpr
      = Ddg_mod_Tensor2_symmetric<A, B, T, U, Dim01_0, Dim23_0, Dim_1, i0, j0,
                                  k0, l0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the % operator");
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01_0, Dim23_0,
                    i0, j0, k0, l0>(TensorExpr(a, b));
  }
}
