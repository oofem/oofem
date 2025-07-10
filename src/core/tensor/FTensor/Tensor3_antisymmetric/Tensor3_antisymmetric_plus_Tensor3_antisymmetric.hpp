/* Adds two Tensor3_antisymmetric's together, yielding a
   Tensor3_antisymmetric. */

#pragma once

namespace FTensor
{
  /* A(i,j,k) + B(i,j,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  class Tensor3_antisymmetric_plus_Tensor3_antisymmetric
  {
    Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> iterA;
    Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, j, k> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) + iterB(N1, N2, N3);
    }

    Tensor3_antisymmetric_plus_Tensor3_antisymmetric(
      const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
      const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, j, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  Tensor3_antisymmetric_Expr<Tensor3_antisymmetric_plus_Tensor3_antisymmetric<
                               A, B, T, U, Dim0, Dim12, i, j, k>,
                             typename promote<T, U>::V, Dim0, Dim12, i, j, k>
  operator+(const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
            const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, j, k> &b)
  {
    using TensorExpr
      = Tensor3_antisymmetric_plus_Tensor3_antisymmetric<A, B, T, U, Dim0,
                                                         Dim12, i, j, k>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim0, Dim12, i, j, k>(TensorExpr(a, b));
  }

  /* A(i,j,k) + B(i,k,j) (which simplifies to A(i,j,k) - B(i,j,k)) */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  class Tensor3_antisymmetric_plus_Tensor3_antisymmetric_12
  {
    Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> iterA;
    Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, k, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) - iterB(N1, N2, N3);
    }

    Tensor3_antisymmetric_plus_Tensor3_antisymmetric_12(
      const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
      const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, k, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  Tensor3_antisymmetric_Expr<
    Tensor3_antisymmetric_plus_Tensor3_antisymmetric_12<A, B, T, U, Dim0,
                                                        Dim12, i, j, k>,
    typename promote<T, U>::V, Dim0, Dim12, i, j, k>
  operator+(const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
            const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, k, j> &b)
  {
    using TensorExpr
      = Tensor3_antisymmetric_plus_Tensor3_antisymmetric_12<A, B, T, U, Dim0,
                                                            Dim12, i, j, k>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim0, Dim12, i, j, k>(TensorExpr(a, b));
  }
}
