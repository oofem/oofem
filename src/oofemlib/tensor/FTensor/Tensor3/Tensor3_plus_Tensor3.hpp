/* Adds a Tensor3 to a Tensor3, yielding a Tensor3. */

#pragma once

namespace FTensor
{
  /* Base Template */
  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim2_0, int Dim0_1, int Dim1_1, int Dim2_1, char i0, char j0,
            char k0, char i1, char j1, char k1>
  class Tensor3_plus_Tensor3
  {};

  /* A(i,j,k)+B(i,j,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor3_plus_Tensor3<A, B, T, U, Dim0, Dim1, Dim2, Dim0, Dim1, Dim2, i,
                             j, k, i, j, k>
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) + iterB(N1, N2, N3);
    }

    Tensor3_plus_Tensor3(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j,k)+B(i,k,j)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor3_plus_Tensor3<A, B, T, U, Dim0, Dim1, Dim2, Dim0, Dim2, Dim1, i,
                             j, k, i, k, j>
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) + iterB(N1, N3, N2);
    }

    Tensor3_plus_Tensor3(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j,k)+B(j,i,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor3_plus_Tensor3<A, B, T, U, Dim0, Dim1, Dim2, Dim1, Dim0, Dim2, i,
                             j, k, j, i, k>
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) + iterB(N2, N1, N3);
    }

    Tensor3_plus_Tensor3(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j,k)+B(j,k,i)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor3_plus_Tensor3<A, B, T, U, Dim0, Dim1, Dim2, Dim1, Dim2, Dim0, i,
                             j, k, j, k, i>
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) + iterB(N2, N3, N1);
    }

    Tensor3_plus_Tensor3(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j,k)+B(k,i,j)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor3_plus_Tensor3<A, B, T, U, Dim0, Dim1, Dim2, Dim2, Dim0, Dim1, i,
                             j, k, k, i, j>
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) + iterB(N3, N1, N2);
    }

    Tensor3_plus_Tensor3(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* A(i,j,k)+B(k,j,i)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor3_plus_Tensor3<A, B, T, U, Dim0, Dim1, Dim2, Dim2, Dim1, Dim0, i,
                             j, k, k, j, i>
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) + iterB(N3, N2, N1);
    }

    Tensor3_plus_Tensor3(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim2_0, int Dim0_1, int Dim1_1, int Dim2_1, char i0, char j0,
            char k0, char i1, char j1, char k1>
  Tensor3_Expr<Tensor3_plus_Tensor3<A, B, T, U, Dim0_0, Dim1_0, Dim2_0, Dim0_1,
                                    Dim1_1, Dim2_1, i0, j0, k0, i1, j1, k1>,
               typename promote<T, U>::V, Dim0_0, Dim1_0, Dim2_0, i0, j0, k0>
  operator+(const Tensor3_Expr<A, T, Dim0_0, Dim1_0, Dim2_0, i0, j0, k0> &a,
            const Tensor3_Expr<B, U, Dim0_1, Dim1_1, Dim2_1, i1, j1, k1> &b)
  {
    using TensorExpr
      = Tensor3_plus_Tensor3<A, B, T, U, Dim0_0, Dim1_0, Dim2_0, Dim0_1,
                             Dim1_1, Dim2_1, i0, j0, k0, i1, j1, k1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the + operator");
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0_0, Dim1_0,
                        Dim2_0, i0, j0, k0>(TensorExpr(a, b));
  }
}
