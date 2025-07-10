/* Fully contracts a Tensor3_antisymmetric with a Tensor3, yielding a
   typename promote<T,U>::V.  I removed the identically zero components of
   Tensor3_antisymmetric.*/

#pragma once

namespace FTensor
{
  /* A(i,j,k)*B(i,j,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_012(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, j, k> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
           + T3as_times_T3_012(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_012(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, j, k> &b,
    const Number<Current_Dim0> &, const Number<1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, 0, Current_Dim2 - 1)
           + T3as_times_T3_012(a, b, Number<Current_Dim0>(), Number<Dim12>(),
                               Number<Current_Dim2 - 1>());
  }

  /* A special case for when the last two indices are equal. */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_012(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, j, k> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim2> &,
    const Number<Current_Dim2> &)
  {
    return T3as_times_T3_012(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim2 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V T3as_times_T3_012(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, j, k> &b,
    const Number<Current_Dim0> &, const Number<2> &, const Number<1> &)
  {
    return a(Current_Dim0 - 1, 1, 0) * b(Current_Dim0 - 1, 1, 0)
           + T3as_times_T3_012(a, b, Number<Current_Dim0 - 1>(),
                               Number<Dim12 - 1>(), Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V T3as_times_T3_012(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, j, k> &b,
    const Number<1> &, const Number<2> &, const Number<1> &)
  {
    return a(0, 1, 0) * b(0, 1, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, j, k> &b)
  {
    return T3as_times_T3_012(a, b, Number<Dim0>(), Number<Dim12 - 1>(),
                             Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, j, k> &b,
            const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a)
  {
    return a * b;
  }

  /* A(i,j,k)*B(k,i,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_201(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, i, j> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim0 - 1, Current_Dim1 - 1)
           + T3as_times_T3_201(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_201(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, i, j> &b,
    const Number<Current_Dim0> &, const Number<1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim0 - 1, 0)
           + T3as_times_T3_201(a, b, Number<Current_Dim0>(), Number<Dim12>(),
                               Number<Current_Dim2 - 1>());
  }

  /* A special case for when the last two indices are equal. */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_201(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, i, j> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim2> &,
    const Number<Current_Dim2> &)
  {
    return T3as_times_T3_201(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim2 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V T3as_times_T3_201(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, i, j> &b,
    const Number<Current_Dim0> &, const Number<2> &, const Number<1> &)
  {
    return a(Current_Dim0 - 1, 1, 0) * b(0, Current_Dim0 - 1, 1)
           + T3as_times_T3_201(a, b, Number<Current_Dim0 - 1>(),
                               Number<Dim12 - 1>(), Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V T3as_times_T3_201(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, i, j> &b,
    const Number<1> &, const Number<2> &, const Number<1> &)
  {
    return a(0, 1, 0) * b(0, 0, 1);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, i, j> &b)
  {
    return T3as_times_T3_201(a, b, Number<Dim0>(), Number<Dim12 - 1>(),
                             Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<B, U, Dim12, Dim0, Dim12, k, i, j> &b,
            const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a)
  {
    return a * b;
  }

  /* A(i,j,k)*B(j,k,i) */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_120(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, k, i> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim1 - 1, Current_Dim2 - 1, Current_Dim0 - 1)
           + T3as_times_T3_120(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_120(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, k, i> &b,
    const Number<Current_Dim0> &, const Number<1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(0, Current_Dim2 - 1, Current_Dim0 - 1)
           + T3as_times_T3_120(a, b, Number<Current_Dim0>(), Number<Dim12>(),
                               Number<Current_Dim2 - 1>());
  }

  /* A special case for when the last two indices are equal. */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_120(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, k, i> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim2> &,
    const Number<Current_Dim2> &)
  {
    return T3as_times_T3_120(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim2 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V T3as_times_T3_120(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, k, i> &b,
    const Number<Current_Dim0> &, const Number<2> &, const Number<1> &)
  {
    return a(Current_Dim0 - 1, 1, 0) * b(1, 0, Current_Dim0 - 1)
           + T3as_times_T3_120(a, b, Number<Current_Dim0 - 1>(),
                               Number<Dim12 - 1>(), Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V T3as_times_T3_120(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, k, i> &b,
    const Number<1> &, const Number<2> &, const Number<1> &)
  {
    return a(0, 1, 0) * b(1, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, k, i> &b)
  {
    return T3as_times_T3_120(a, b, Number<Dim0>(), Number<Dim12 - 1>(),
                             Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<B, U, Dim12, Dim12, Dim0, j, k, i> &b,
            const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a)
  {
    return a * b;
  }

  /* A(i,j,k)*B(j,i,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_102(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, i, k> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim1 - 1, Current_Dim0 - 1, Current_Dim2 - 1)
           + T3as_times_T3_102(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_102(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, i, k> &b,
    const Number<Current_Dim0> &, const Number<1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(0, Current_Dim0 - 1, Current_Dim2 - 1)
           + T3as_times_T3_102(a, b, Number<Current_Dim0>(), Number<Dim12>(),
                               Number<Current_Dim2 - 1>());
  }

  /* A special case for when the last two indices are equal. */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_102(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, i, k> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim2> &,
    const Number<Current_Dim2> &)
  {
    return T3as_times_T3_102(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim2 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V T3as_times_T3_102(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, i, k> &b,
    const Number<Current_Dim0> &, const Number<2> &, const Number<1> &)
  {
    return a(Current_Dim0 - 1, 1, 0) * b(1, Current_Dim0 - 1, 0)
           + T3as_times_T3_102(a, b, Number<Current_Dim0 - 1>(),
                               Number<Dim12 - 1>(), Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V T3as_times_T3_102(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, i, k> &b,
    const Number<1> &, const Number<2> &, const Number<1> &)
  {
    return a(0, 1, 0) * b(1, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, j, i, k> &b)
  {
    return T3as_times_T3_102(a, b, Number<Dim0>(), Number<Dim12 - 1>(),
                             Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<B, U, Dim12, Dim0, Dim12, j, i, k> &b,
            const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a)
  {
    return a * b;
  }

  /* A(i,j,k)*B(k,j,i) */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_210(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, j, i> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, Current_Dim1 - 1, Current_Dim0 - 1)
           + T3as_times_T3_210(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_210(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, j, i> &b,
    const Number<Current_Dim0> &, const Number<1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim2 - 1, 0, Current_Dim0 - 1)
           + T3as_times_T3_210(a, b, Number<Current_Dim0>(), Number<Dim12>(),
                               Number<Current_Dim2 - 1>());
  }

  /* A special case for when the last two indices are equal. */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_210(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, j, i> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim2> &,
    const Number<Current_Dim2> &)
  {
    return T3as_times_T3_210(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim2 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V T3as_times_T3_210(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, j, i> &b,
    const Number<Current_Dim0> &, const Number<2> &, const Number<1> &)
  {
    return a(Current_Dim0 - 1, 1, 0) * b(0, 1, Current_Dim0 - 1)
           + T3as_times_T3_210(a, b, Number<Current_Dim0 - 1>(),
                               Number<Dim12 - 1>(), Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V T3as_times_T3_210(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, j, i> &b,
    const Number<1> &, const Number<2> &, const Number<1> &)
  {
    return a(0, 1, 0) * b(0, 1, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, k, j, i> &b)
  {
    return T3as_times_T3_210(a, b, Number<Dim0>(), Number<Dim12 - 1>(),
                             Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<B, U, Dim12, Dim12, Dim0, k, j, i> &b,
            const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a)
  {
    return a * b;
  }

  /* A(i,j,k)*B(i,k,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim1,
            int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_012(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, k, j> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim2 - 1, Current_Dim1 - 1)
           + T3as_times_T3_012(a, b, Number<Current_Dim0>(),
                               Number<Current_Dim1 - 1>(),
                               Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_012(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, k, j> &b,
    const Number<Current_Dim0> &, const Number<1> &,
    const Number<Current_Dim2> &)
  {
    return a(Current_Dim0 - 1, 0, Current_Dim2 - 1)
             * b(Current_Dim0 - 1, Current_Dim2 - 1, 0)
           + T3as_times_T3_012(a, b, Number<Current_Dim0>(), Number<Dim12>(),
                               Number<Current_Dim2 - 1>());
  }

  /* A special case for when the last two indices are equal. */

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0, int Current_Dim2>
  typename promote<T, U>::V T3as_times_T3_012(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, k, j> &b,
    const Number<Current_Dim0> &, const Number<Current_Dim2> &,
    const Number<Current_Dim2> &)
  {
    return T3as_times_T3_012(a, b, Number<Current_Dim0>(),
                             Number<Current_Dim2 - 1>(),
                             Number<Current_Dim2>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k, int Current_Dim0>
  typename promote<T, U>::V T3as_times_T3_012(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, k, j> &b,
    const Number<Current_Dim0> &, const Number<2> &, const Number<1> &)
  {
    return a(Current_Dim0 - 1, 1, 0) * b(Current_Dim0 - 1, 0, 1)
           + T3as_times_T3_012(a, b, Number<Current_Dim0 - 1>(),
                               Number<Dim12 - 1>(), Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V T3as_times_T3_012(
    const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
    const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, k, j> &b,
    const Number<1> &, const Number<2> &, const Number<1> &)
  {
    return a(0, 1, 0) * b(0, 0, 1);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a,
            const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, k, j> &b)
  {
    return T3as_times_T3_012(a, b, Number<Dim0>(), Number<Dim12 - 1>(),
                             Number<Dim12>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim12, char i,
            char j, char k>
  typename promote<T, U>::V
  operator*(const Tensor3_Expr<B, U, Dim0, Dim12, Dim12, i, k, j> &b,
            const Tensor3_antisymmetric_Expr<A, T, Dim0, Dim12, i, j, k> &a)
  {
    return a * b;
  }
}
