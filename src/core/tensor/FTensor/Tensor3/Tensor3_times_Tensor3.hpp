/* Fully contracts a Tensor3 with a Tensor3, yielding a typename
   promote<T,U>::V. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)*B(i,j,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_012(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
           + T3_times_T3_012(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_012(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, 0, Current_Dim2 - 1)
           + T3_times_T3_012(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_012(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(Current_Dim0 - 1, 0, 0)
           + T3_times_T3_012(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_012(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b)
  {
    return T3_times_T3_012(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(k,i,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_201(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim0 - 1, Current_Dim1 - 1)
           + T3_times_T3_201(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_201(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim0 - 1, 0)
           + T3_times_T3_201(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_201(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, Current_Dim0 - 1, 0)
           + T3_times_T3_201(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_201(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b)
  {
    return T3_times_T3_201(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(j,k,i) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_120(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim1 - 1, Current_Dim2 - 1, Current_Dim0 - 1)
           + T3_times_T3_120(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_120(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(0, Current_Dim2 - 1, Current_Dim0 - 1)
           + T3_times_T3_120(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_120(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, 0, Current_Dim0 - 1)
           + T3_times_T3_120(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_120(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b)
  {
    return T3_times_T3_120(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(j,i,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_102(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim1 - 1, Current_Dim0 - 1, Current_Dim2 - 1)
           + T3_times_T3_102(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_102(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(0, Current_Dim0 - 1, Current_Dim2 - 1)
           + T3_times_T3_102(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_102(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, Current_Dim0 - 1, 0)
           + T3_times_T3_102(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_102(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b)
  {
    return T3_times_T3_102(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(k,j,i) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_210(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim1 - 1, Current_Dim0 - 1)
           + T3_times_T3_210(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_210(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, 0, Current_Dim0 - 1)
           + T3_times_T3_210(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_210(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(0, 0, Current_Dim0 - 1)
           + T3_times_T3_210(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_210(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b)
  {
    return T3_times_T3_210(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(i,k,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_021(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim2 - 1, Current_Dim1 - 1)
           + T3_times_T3_021(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim1 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V
  T3_times_T3_021(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim2 - 1, 0)
           + T3_times_T3_021(a, b, Number<Current_Dim0>(), Number<Dim1>(),
                             Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int Current_Dim0>
  typename promote<T, U>::V
  T3_times_T3_021(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
                  const Number<Current_Dim0> &, const Number<1> &,
                  const Number<1> &)
  {
    return a(Current_Dim0 - 1, 0, 0) * b(Current_Dim0 - 1, 0, 0)
           + T3_times_T3_021(a, b, Number<Current_Dim0 - 1>(), Number<Dim1>(),
                             Number<Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  T3_times_T3_021(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                  const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
                  const Number<1> &, const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0) * b(0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b)
  {
    return T3_times_T3_021(a, b, Number<Dim0>(), Number<Dim1>(),
                           Number<Dim2>());
  }

  /* A(i,j,k)*B(k,l,m) -> Tensor 4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim23,
            int Dim4, int Dim5, char i, char j, char k, char l, char m>
  class Tensor3_times_Tensor3_21
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim23, i, j, k> iterA;
    Tensor3_Expr<B, U, Dim23, Dim4, Dim5, k, l, m> iterB;

    template <int CurrentDim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const int N4,
         const Number<CurrentDim> &) const
    {
      return iterA(N1, N2, CurrentDim - 1) * iterB(CurrentDim - 1, N3, N4)
             + eval(N1, N2, N3, N4, Number<CurrentDim - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const int N4, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(0, N3, N4);
    }

  public:
    Tensor3_times_Tensor3_21(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim23, i, j, k> &a,
      const Tensor3_Expr<B, U, Dim23, Dim4, Dim5, k, l, m> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int &N1, const int &N2,
                                         const int &N3, const int &N4) const
    {
      return eval(N1, N2, N3, N4, Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim23,
            int Dim4, int Dim5, char i, char j, char k, char l, char m>
  Tensor4_Expr<Tensor3_times_Tensor3_21<A, B, T, U, Dim0, Dim1, Dim23, Dim4,
                                        Dim5, i, j, k, l, m>,
               typename promote<T, U>::V, Dim0, Dim1, Dim4, Dim5, i, j, l, m>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim23, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim23, Dim4, Dim5, k, l, m> &b)
  {
    using TensorExpr = Tensor3_times_Tensor3_21<A, B, T, U, Dim0, Dim1, Dim23,
                                                Dim4, Dim5, i, j, k, l, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim4, Dim5, i, j, l, m>(TensorExpr(a, b));
  };
}
