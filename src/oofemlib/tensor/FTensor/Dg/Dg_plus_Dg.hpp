/* Adds two Dg's together, yielding a Dg or Tensor3. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)+B(i,j,k)->Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  class Dg_plus_Dg
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Dg_Expr<B, U, Dim01, Dim2, i, j, k> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) + iterB(N1, N2, N3);
    }

    Dg_plus_Dg(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
               const Dg_Expr<B, U, Dim01, Dim2, i, j, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Dg_Expr<Dg_plus_Dg<A, B, T, U, Dim01, Dim2, i, j, k>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, k>
  operator+(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Dg_Expr<B, U, Dim01, Dim2, i, j, k> &b)
  {
    using TensorExpr = Dg_plus_Dg<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, k>(
      TensorExpr(a, b));
  }

  /* A(i,j,k)+B(j,i,k)->Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  class Dg_plus_Dg_10
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Dg_Expr<B, U, Dim01, Dim2, j, i, k> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) + iterB(N1, N2, N3);
    }

    Dg_plus_Dg_10(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
                  const Dg_Expr<B, U, Dim01, Dim2, j, i, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Dg_Expr<Dg_plus_Dg_10<A, B, T, U, Dim01, Dim2, i, j, k>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, k>
  operator+(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Dg_Expr<B, U, Dim01, Dim2, j, i, k> &b)
  {
    using TensorExpr = Dg_plus_Dg_10<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, k>(
      TensorExpr(a, b));
  }

  /* A(i,j,k)+B(k,j,i)->Tensor3 */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Dg_plus_Dg_210
  {
    Dg_Expr<A, T, Dim, Dim, i, j, k> iterA;
    Dg_Expr<B, U, Dim, Dim, k, j, i> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3) + iterB(N3, N2, N1);
    }

    Dg_plus_Dg_210(const Dg_Expr<A, T, Dim, Dim, i, j, k> &a,
                   const Dg_Expr<B, U, Dim, Dim, k, j, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor3_Expr<Dg_plus_Dg_210<A, B, T, U, Dim, i, j, k>,
               typename promote<T, U>::V, Dim, Dim, Dim, i, j, k>
  operator+(const Dg_Expr<A, T, Dim, Dim, i, j, k> &a,
            const Dg_Expr<B, U, Dim, Dim, k, j, i> &b)
  {
    using TensorExpr = Dg_plus_Dg_210<A, B, T, U, Dim, i, j, k>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim, Dim,
                        i, j, k>(TensorExpr(a, b));
  }
}
