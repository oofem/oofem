/* Adds Ddg+Ddg -> Ddg */

#pragma once

namespace FTensor
{
  /* Base Template */
  template <class A, class B, class T, class U, int Dim01_0, int Dim23_0,
            int Dim01_1, int Dim23_1, char i0, char j0, char k0, char l0,
            char i1, char j1, char k1, char l1>
  class Ddg_plus_Ddg
  {};

  /* A(i,j,k,l)+B(i,j,k,l) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_plus_Ddg<A, B, T, U, Dim01, Dim23, Dim01, Dim23, i, j, k, l, i, j,
                     k, l>
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) + iterB(N1, N2, N3, N4);
    }

    Ddg_plus_Ddg(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                 const Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j,k,l)+B(k,l,i,j) -> Ddg */
  // TODO: Add rest of combinations (k,l,j,i), (l,k,i,j), (l,k,j,i)
  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_plus_Ddg<A, B, T, U, Dim01, Dim23, Dim23, Dim01, i, j, k, l, k, l,
                     i, j>
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Ddg_Expr<B, U, Dim23, Dim01, k, l, i, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) + iterB(N3, N4, N1, N2);
    }

    Ddg_plus_Ddg(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                 const Ddg_Expr<B, U, Dim23, Dim01, k, l, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim01_0, int Dim23_0,
            int Dim01_1, int Dim23_1, char i0, char j0, char k0, char l0,
            char i1, char j1, char k1, char l1>
  Ddg_Expr<Ddg_plus_Ddg<A, B, T, U, Dim01_0, Dim23_0, Dim01_1, Dim23_1, i0, j0,
                        k0, l0, i1, j1, k1, l1>,
           typename promote<T, U>::V, Dim01_0, Dim23_0, i0, j0, k0, l0>
  operator+(const Ddg_Expr<A, T, Dim01_0, Dim23_0, i0, j0, k0, l0> &a,
            const Ddg_Expr<B, U, Dim01_1, Dim23_1, i1, j1, k1, l1> &b)
  {
    using TensorExpr = Ddg_plus_Ddg<A, B, T, U, Dim01_0, Dim23_0, Dim01_1,
                                    Dim23_1, i0, j0, k0, l0, i1, j1, k1, l1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the + operator");
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01_0, Dim23_0,
                    i0, j0, k0, l0>(TensorExpr(a, b));
  }
}
