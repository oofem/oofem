/* This file has all of the declarations for expressions like
   Dg*Dg yielding a Tensor2 or Ddg. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)*B(j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  class Dg_times_Dg_12_01
  {
    Dg_Expr<A, T, Dim01, Dim01, i, j, k> iterA;
    Dg_Expr<B, U, Dim01, Dim2, j, k, l> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1, N2)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, 0, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1, N2)
             + eval(N1, N2, Number<Dim01>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, 0, 0) * iterB(0, 0, N2);
    }

  public:
    Dg_times_Dg_12_01(const Dg_Expr<A, T, Dim01, Dim01, i, j, k> &a,
                      const Dg_Expr<B, U, Dim01, Dim2, j, k, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim01>(), Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Tensor2_Expr<Dg_times_Dg_12_01<A, B, T, U, Dim01, Dim2, i, j, k, l>,
               typename promote<T, U>::V, Dim01, Dim2, i, l>
  operator*(const Dg_Expr<A, T, Dim01, Dim01, i, j, k> &a,
            const Dg_Expr<B, U, Dim01, Dim2, j, k, l> &b)
  {
    using TensorExpr = Dg_times_Dg_12_01<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i,
                        l>(TensorExpr(a, b));
  }

  /* B(j,k,l)*A(i,j,k) */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Tensor2_Expr<Dg_times_Dg_12_01<A, B, T, U, Dim01, Dim2, i, j, k, l>,
               typename promote<T, U>::V, Dim01, Dim2, i, l>
  operator*(const Dg_Expr<B, U, Dim01, Dim2, j, k, l> &b,
            const Dg_Expr<A, T, Dim01, Dim01, i, j, k> &a)
  {
    using TensorExpr = Dg_times_Dg_12_01<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i,
                        l>(TensorExpr(a, b));
  }

  /* A(i,j,k)*B(k,l,j) */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  class Dg_times_Dg_12_20
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Dg_Expr<B, U, Dim2, Dim01, k, l, j> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim1 - 1, N2, Current_Dim0 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, 0, Current_Dim1 - 1) * iterB(Current_Dim1 - 1, N2, 0)
             + eval(N1, N2, Number<Dim01>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, 0, 0) * iterB(0, N2, 0);
    }

  public:
    Dg_times_Dg_12_20(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
                      const Dg_Expr<B, U, Dim2, Dim01, k, l, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim01>(), Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Tensor2_Expr<Dg_times_Dg_12_20<A, B, T, U, Dim01, Dim2, i, j, k, l>,
               typename promote<T, U>::V, Dim01, Dim2, i, l>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Dg_Expr<B, U, Dim2, Dim01, k, l, j> &b)
  {
    using TensorExpr = Dg_times_Dg_12_20<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i,
                        l>(TensorExpr(a, b));
  }

  /* B(k,l,j)*A(i,j,k) */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Tensor2_Expr<Dg_times_Dg_12_20<A, B, T, U, Dim01, Dim2, i, j, k, l>,
               typename promote<T, U>::V, Dim01, Dim2, i, l>
  operator*(const Dg_Expr<B, U, Dim2, Dim01, k, l, j> &b,
            const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a)
  {
    using TensorExpr = Dg_times_Dg_12_20<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i,
                        l>(TensorExpr(a, b));
  }

  /* A(i,j,k)*B(l,m,k) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
            char i, char j, char k, char l, char m>
  class Dg_times_Dg_2
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Dg_Expr<B, U, Dim23, Dim2, l, m, k> iterB;

    template <int Current_Dim0>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<Current_Dim0> &) const
    {
      return iterA(N1, N2, Current_Dim0 - 1) * iterB(N3, N4, Current_Dim0 - 1)
             + eval(N1, N2, N3, N4, Number<Current_Dim0 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(N3, N4, 0);
    }

  public:
    Dg_times_Dg_2(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
                  const Dg_Expr<B, U, Dim23, Dim2, l, m, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, int Dim2,
            char i, char j, char k, char l, char m>
  Ddg_Expr<Dg_times_Dg_2<A, B, T, U, Dim01, Dim23, Dim2, i, j, k, l, m>,
           typename promote<T, U>::V, Dim01, Dim23, i, j, l, m>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Dg_Expr<B, U, Dim23, Dim2, l, m, k> &b)
  {
    using TensorExpr
      = Dg_times_Dg_2<A, B, T, U, Dim01, Dim23, Dim2, i, j, k, l, m>;
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j,
                    l, m>(TensorExpr(a, b));
  }
}
