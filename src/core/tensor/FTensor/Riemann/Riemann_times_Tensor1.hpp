/* This file has all of the declarations for expressions like
   Riemann*Tensor1, yielding a Tensor3_antisymmetric. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(i) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  class Riemann_times_Tensor1_0
  {
    Riemann_Expr<A, T, Dim, i, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim, i> iterB;

  public:
    Riemann_times_Tensor1_0(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
                            const Tensor1_Expr<B, U, Dim, i> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(0, N1, N2, N3) * iterB(0) + iterA(1, N1, N2, N3) * iterB(1)
             + iterA(2, N1, N2, N3) * iterB(2);
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Tensor3_antisymmetric_Expr<
    Riemann_times_Tensor1_0<A, B, T, U, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim, Dim, j, k, l>
  operator*(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim, i> &b)
  {
    using TensorExpr = Riemann_times_Tensor1_0<A, B, T, U, Dim, i, j, k, l>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim, Dim, j, k, l>(TensorExpr(a, b));
  }

  /* B(i)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Tensor3_antisymmetric_Expr<
    Riemann_times_Tensor1_0<A, B, T, U, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim, Dim, j, k, l>
  operator*(const Tensor1_Expr<B, U, Dim, i> &b,
            const Riemann_Expr<A, T, Dim, i, j, k, l> &a)
  {
    using TensorExpr = Riemann_times_Tensor1_0<A, B, T, U, Dim, i, j, k, l>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim, Dim, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(j) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  class Riemann_times_Tensor1_1
  {
    Riemann_Expr<A, T, Dim, i, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim, j> iterB;

  public:
    Riemann_times_Tensor1_1(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
                            const Tensor1_Expr<B, U, Dim, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, 0, N2, N3) * iterB(0) + iterA(N1, 1, N2, N3) * iterB(1)
             + iterA(N1, 2, N2, N3) * iterB(2);
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Tensor3_antisymmetric_Expr<
    Riemann_times_Tensor1_1<A, B, T, U, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim, Dim, i, k, l>
  operator*(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim, j> &b)
  {
    using TensorExpr = Riemann_times_Tensor1_1<A, B, T, U, Dim, i, j, k, l>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim, Dim, i, k, l>(TensorExpr(a, b));
  }

  /* B(j)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Tensor3_antisymmetric_Expr<
    Riemann_times_Tensor1_1<A, B, T, U, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim, Dim, i, k, l>
  operator*(const Tensor1_Expr<B, U, Dim, j> &b,
            const Riemann_Expr<A, T, Dim, i, j, k, l> &a)
  {
    using TensorExpr = Riemann_times_Tensor1_1<A, B, T, U, Dim, i, j, k, l>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim, Dim, i, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(k) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  class Riemann_times_Tensor1_2
  {
    Riemann_Expr<A, T, Dim, i, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim, k> iterB;

  public:
    Riemann_times_Tensor1_2(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
                            const Tensor1_Expr<B, U, Dim, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, 0, N3) * iterB(0) + iterA(N1, N2, 1, N3) * iterB(1)
             + iterA(N1, N2, 2, N3) * iterB(2);
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Tensor3_antisymmetric_Expr<
    Riemann_times_Tensor1_2<A, B, T, U, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim, Dim, i, j, l>
  operator*(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim, k> &b)
  {
    using TensorExpr = Riemann_times_Tensor1_2<A, B, T, U, Dim, i, j, k, l>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim, Dim, i, j, l>(TensorExpr(a, b));
  }

  /* B(k)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Tensor3_antisymmetric_Expr<
    Riemann_times_Tensor1_2<A, B, T, U, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim, Dim, i, j, l>
  operator*(const Tensor1_Expr<B, U, Dim, k> &b,
            const Riemann_Expr<A, T, Dim, i, j, k, l> &a)
  {
    using TensorExpr = Riemann_times_Tensor1_2<A, B, T, U, Dim, i, j, k, l>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim, Dim, i, j, l>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  class Riemann_times_Tensor1_3
  {
    Riemann_Expr<A, T, Dim, i, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim, l> iterB;

  public:
    Riemann_times_Tensor1_3(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
                            const Tensor1_Expr<B, U, Dim, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2, N3, 0) * iterB(0) + iterA(N1, N2, N3, 1) * iterB(1)
             + iterA(N1, N2, N3, 2) * iterB(2);
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Tensor3_antisymmetric_Expr<
    Riemann_times_Tensor1_3<A, B, T, U, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim, Dim, i, j, k>
  operator*(const Riemann_Expr<A, T, Dim, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim, l> &b)
  {
    using TensorExpr = Riemann_times_Tensor1_3<A, B, T, U, Dim, i, j, k, l>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim, Dim, i, j, k>(TensorExpr(a, b));
  }

  /* B(l)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  Tensor3_antisymmetric_Expr<
    Riemann_times_Tensor1_3<A, B, T, U, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim, Dim, i, j, k>
  operator*(const Tensor1_Expr<B, U, Dim, l> &b,
            const Riemann_Expr<A, T, Dim, i, j, k, l> &a)
  {
    using TensorExpr = Riemann_times_Tensor1_3<A, B, T, U, Dim, i, j, k, l>;
    return Tensor3_antisymmetric_Expr<TensorExpr, typename promote<T, U>::V,
                                      Dim, Dim, i, j, k>(TensorExpr(a, b));
  }
}
