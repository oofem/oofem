/* Contracts a Tensor3 with a Dg, yielding a typename promote<T,U>::V,
   Tensor2, or Tensor4. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)*B(i,j,k) */

  template <class A, class B, class T, class U, int Dim, int Dim2, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_012(const Tensor3_Expr<A, T, Dim, Dim, Dim2, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim2, i, j, k> &b,
                    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
           + T3_times_T3dg_012(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim2, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_012(const Tensor3_Expr<A, T, Dim, Dim, Dim2, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim2, i, j, k> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, 0, Current_Dim2 - 1)
           + T3_times_T3dg_012(a, b, Number<Current_Dim0>(), Number<Dim>(),
                               Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim2, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3dg_012(const Tensor3_Expr<A, T, Dim, Dim, Dim2, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim2, i, j, k> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(Current_Dim0 - 1, 0, 0)
           + T3_times_T3dg_012(a, b, Number<Current_Dim0 - 1>(), Number<Dim>(),
                               Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim2, char i,
            char j, char k>
  typename promote<T, U>::V
  T3_times_T3dg_012(const Tensor3_Expr<A, T, Dim, Dim, Dim2, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim2, i, j, k> &b,
                    const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim, int Dim2, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim, Dim, Dim2, i, j, k> &a,
            const Dg_Expr<B, U, Dim, Dim2, i, j, k> &b)
  {
    return T3_times_T3dg_012(a, b, Number<Dim>(), Number<Dim>(),
                             Number<Dim2>());
  }

  /* A(i,j,k)*B(k,i,j) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_201(const Tensor3_Expr<A, T, Dim, Dim1, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim1, k, i, j> &b,
                    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim0 - 1, Current_Dim1 - 1)
           + T3_times_T3dg_201(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_201(const Tensor3_Expr<A, T, Dim, Dim1, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim1, k, i, j> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim0 - 1, 0)
           + T3_times_T3dg_201(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                               Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3dg_201(const Tensor3_Expr<A, T, Dim, Dim1, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim1, k, i, j> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, Current_Dim0 - 1, 0)
           + T3_times_T3dg_201(a, b, Number<Current_Dim0 - 1>(),
                               Number<Dim1>(), Number<Dim>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  typename promote<T, U>::V
  T3_times_T3dg_201(const Tensor3_Expr<A, T, Dim, Dim1, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim1, k, i, j> &b,
                    const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim, Dim1, Dim, i, j, k> &a,
            const Dg_Expr<B, U, Dim, Dim1, k, i, j> &b)
  {
    return T3_times_T3dg_201(a, b, Number<Dim>(), Number<Dim1>(),
                             Number<Dim>());
  }

  /* A(i,j,k)*B(j,k,i) */

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_120(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim0, j, k, i> &b,
                    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim1 - 1, Current_Dim2 - 1, Current_Dim0 - 1)
           + T3_times_T3dg_120(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_120(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim0, j, k, i> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(0, Current_Dim2 - 1, Current_Dim0 - 1)
           + T3_times_T3dg_120(a, b, Number<Current_Dim0>(), Number<Dim>(),
                               Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3dg_120(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim0, j, k, i> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, 0, Current_Dim0 - 1)
           + T3_times_T3dg_120(a, b, Number<Current_Dim0 - 1>(), Number<Dim>(),
                               Number<Dim>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k>
  typename promote<T, U>::V
  T3_times_T3dg_120(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim0, j, k, i> &b,
                    const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
            const Dg_Expr<B, U, Dim, Dim0, j, k, i> &b)
  {
    return T3_times_T3dg_120(a, b, Number<Dim0>(), Number<Dim>(),
                             Number<Dim>());
  }

  /* A(i,j,k)*B(j,i,k) */

  template <class A, class B, class T, class U, int Dim, int Dim2, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_102(const Tensor3_Expr<A, T, Dim, Dim, Dim2, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim2, j, i, k> &b,
                    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim1 - 1, Current_Dim0 - 1, Current_Dim2 - 1)
           + T3_times_T3dg_102(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim2, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_102(const Tensor3_Expr<A, T, Dim, Dim, Dim2, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim2, j, i, k> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(0, Current_Dim0 - 1, Current_Dim2 - 1)
           + T3_times_T3dg_102(a, b, Number<Current_Dim0>(), Number<Dim>(),
                               Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim2, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3dg_102(const Tensor3_Expr<A, T, Dim, Dim, Dim2, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim2, j, i, k> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, Current_Dim0 - 1, 0)
           + T3_times_T3dg_102(a, b, Number<Current_Dim0 - 1>(), Number<Dim>(),
                               Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim2, char i,
            char j, char k>
  typename promote<T, U>::V
  T3_times_T3dg_102(const Tensor3_Expr<A, T, Dim, Dim, Dim2, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim2, j, i, k> &b,
                    const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim, int Dim2, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim, Dim, Dim2, i, j, k> &a,
            const Dg_Expr<B, U, Dim, Dim2, j, i, k> &b)
  {
    return T3_times_T3dg_102(a, b, Number<Dim>(), Number<Dim>(),
                             Number<Dim2>());
  }

  /* A(i,j,k)*B(k,j,i) */

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_210(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim0, k, j, i> &b,
                    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim1 - 1, Current_Dim0 - 1)
           + T3_times_T3dg_210(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_210(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim0, k, j, i> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, 0, Current_Dim0 - 1)
           + T3_times_T3dg_210(a, b, Number<Current_Dim0>(), Number<Dim>(),
                               Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3dg_210(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim0, k, j, i> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, 0, Current_Dim0 - 1)
           + T3_times_T3dg_210(a, b, Number<Current_Dim0 - 1>(), Number<Dim>(),
                               Number<Dim>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k>
  typename promote<T, U>::V
  T3_times_T3dg_210(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim0, k, j, i> &b,
                    const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
            const Dg_Expr<B, U, Dim, Dim0, k, j, i> &b)
  {
    return T3_times_T3dg_210(a, b, Number<Dim0>(), Number<Dim>(),
                             Number<Dim>());
  }

  /* A(i,j,k)*B(i,k,j) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_021(const Tensor3_Expr<A, T, Dim, Dim1, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim1, i, k, j> &b,
                    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim2 - 1, Current_Dim1 - 1)
           + T3_times_T3dg_021(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3dg_021(const Tensor3_Expr<A, T, Dim, Dim1, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim1, i, k, j> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim2 - 1, 0)
           + T3_times_T3dg_021(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                               Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3dg_021(const Tensor3_Expr<A, T, Dim, Dim1, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim1, i, k, j> &b,
                    const Number<Current_Dim0> &, const Number<1> &,
                    const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(Current_Dim0 - 1, 0, 0)
           + T3_times_T3dg_021(a, b, Number<Current_Dim0 - 1>(),
                               Number<Dim1>(), Number<Dim>());
  }

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  typename promote<T, U>::V
  T3_times_T3dg_021(const Tensor3_Expr<A, T, Dim, Dim1, Dim, i, j, k> &a,
                    const Dg_Expr<B, U, Dim, Dim1, i, k, j> &b,
                    const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim, Dim1, Dim, i, j, k> &a,
            const Dg_Expr<B, U, Dim, Dim1, i, k, j> &b)
  {
    return T3_times_T3dg_021(a, b, Number<Dim>(), Number<Dim1>(),
                             Number<Dim>());
  }

  /* A(j,i,k)*B(k,l,j) */

  template <class A, class B, class T, class U, int Dim2, int Dim1, int Dim01,
            char i, char j, char k, char l>
  class Tensor3_times_Dg_02_20
  {
    Tensor3_Expr<A, T, Dim2, Dim1, Dim01, j, i, k> iterA;
    Dg_Expr<B, U, Dim01, Dim2, k, l, j> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, N1, Current_Dim1 - 1)
               * iterB(Current_Dim1 - 1, N2, Current_Dim0 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(0, N1, Current_Dim1 - 1) * iterB(Current_Dim1 - 1, N2, 0)
             + eval(N1, N2, Number<Dim2>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(0, N1, 0) * iterB(0, N2, 0);
    }

  public:
    Tensor3_times_Dg_02_20(
      const Tensor3_Expr<A, T, Dim2, Dim1, Dim01, j, i, k> &a,
      const Dg_Expr<B, U, Dim01, Dim2, k, l, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim2>(), Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim2, int Dim1, int Dim01,
            char i, char j, char k, char l>
  Tensor2_Expr<Tensor3_times_Dg_02_20<A, B, T, U, Dim2, Dim1, Dim01, i, j, k, l>,
               typename promote<T, U>::V, Dim1, Dim01, i, l>
  operator*(const Tensor3_Expr<A, T, Dim2, Dim1, Dim01, j, i, k> &a,
            const Dg_Expr<B, U, Dim01, Dim2, k, l, j> &b)
  {
    using TensorExpr
      = Tensor3_times_Dg_02_20<A, B, T, U, Dim2, Dim1, Dim01, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim01, i,
                        l>(TensorExpr(a, b));
  }

  /* B(k,l,j)*A(j,i,k) */

  /* A(i,j,k)*B(l,k,m) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim, char i, char j, char k, char l, char m>
  class Tensor3_times_Dg
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim, i, j, k> iterA;
    Dg_Expr<B, U, Dim, Dim2, l, k, m> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, Current_Dim - 1) * iterB(N3, Current_Dim - 1, N4)
             + eval(N1, N2, N3, N4, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(N3, 0, N4);
    }

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim>());
    }
    Tensor3_times_Dg(const Tensor3_Expr<A, T, Dim0, Dim1, Dim, i, j, k> &a,
                     const Dg_Expr<B, U, Dim, Dim2, l, k, m> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim, char i, char j, char k, char l, char m>
  Tensor4_Expr<
    Tensor3_times_Dg<A, B, T, U, Dim0, Dim1, Dim, Dim2, i, j, k, l, m>,
    typename promote<T, U>::V, Dim0, Dim1, Dim, Dim2, i, j, l, m>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim, i, j, k> &a,
            const Dg_Expr<B, U, Dim, Dim2, l, k, m> &b)
  {
    using TensorExpr
      = Tensor3_times_Dg<A, B, T, U, Dim0, Dim1, Dim2, Dim, i, j, k, l, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim,
                        Dim2, i, j, l, m>(TensorExpr(a, b));
  }
}
