/* This file has all of the declarations for
   Tensor2_symmetric*Tensor2.  This includes the double
   contraction A(i,j)*B(i,j) (yielding a double) as well as the more
   complicated single contraction A(i,j)*B(j,k) (yielding a Tensor2
   expression). */

#pragma once

namespace FTensor
{
  /* Double contraction. */

  /* A(i,j)*B(i,j) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            int Current_Dim0, int Current_Dim1>
  typename promote<T, U>::V
  T2s_times_T2_01(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                  const Tensor2_Expr<B, U, Dim, Dim, i, j> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1)
             * b(Current_Dim0 - 1, Current_Dim1 - 1)
           + T2s_times_T2_01(a, b, Number<Current_Dim0 - 1>(),
                             Number<Current_Dim1>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            int Current_Dim1>
  typename promote<T, U>::V
  T2s_times_T2_01(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                  const Tensor2_Expr<B, U, Dim, Dim, i, j> &b,
                  const Number<1> &, const Number<Current_Dim1> &)
  {
    return a(0, Current_Dim1 - 1) * b(0, Current_Dim1 - 1)
           + T2s_times_T2_01(a, b, Number<Dim>(), Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j>
  typename promote<T, U>::V
  T2s_times_T2_01(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                  const Tensor2_Expr<B, U, Dim, Dim, i, j> &b,
                  const Number<1> &, const Number<1> &)
  {
    return a(0, 0) * b(0, 0);
  }

  template <class A, class B, class T, class U, int Dim, char i, char j>
  typename promote<T, U>::V
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
            const Tensor2_Expr<B, U, Dim, Dim, i, j> &b)
  {
    return T2s_times_T2_01(a, b, Number<Dim>(), Number<Dim>());
  }

  /* B(i,j)*A(i,j) */

  template <class A, class B, class T, class U, int Dim, char i, char j>
  typename promote<T, U>::V
  operator*(const Tensor2_Expr<B, U, Dim, Dim, i, j> &b,
            const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a)

  {
    return T2s_times_T2_01(a, b, Number<Dim>(), Number<Dim>());
  }

  /* Double contraction with switched indices.  The perspicacious reader
     will note that it just replicated code for the aligned indices case
     above, but with B(i,j) -> B(j,i).  This is ok, since A is
     symmetric. */

  /* A(i,j)*B(j,i) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            int Current_Dim0, int Current_Dim1>
  typename promote<T, U>::V
  T2s_times_T2_01(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                  const Tensor2_Expr<B, U, Dim, Dim, j, i> &b,
                  const Number<Current_Dim0> &, const Number<Current_Dim1> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1)
             * b(Current_Dim0 - 1, Current_Dim1 - 1)
           + T2s_times_T2_01(a, b, Number<Current_Dim0 - 1>(),
                             Number<Current_Dim1>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            int Current_Dim1>
  typename promote<T, U>::V
  T2s_times_T2_01(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                  const Tensor2_Expr<B, U, Dim, Dim, j, i> &b,
                  const Number<1> &, const Number<Current_Dim1> &)
  {
    return a(0, Current_Dim1 - 1) * b(0, Current_Dim1 - 1)
           + T2s_times_T2_01(a, b, Number<Dim>(), Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j>
  typename promote<T, U>::V
  T2s_times_T2_01(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                  const Tensor2_Expr<B, U, Dim, Dim, j, i> &b,
                  const Number<1> &, const Number<1> &)
  {
    return a(0, 0) * b(0, 0);
  }

  template <class A, class B, class T, class U, int Dim, char i, char j>
  typename promote<T, U>::V
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
            const Tensor2_Expr<B, U, Dim, Dim, j, i> &b)
  {
    return T2s_times_T2_01(a, b, Number<Dim>(), Number<Dim>());
  }

  /* B(j,i)*A(i,j) */

  template <class A, class B, class T, class U, int Dim, char i, char j>
  typename promote<T, U>::V
  operator*(const Tensor2_Expr<B, U, Dim, Dim, j, i> &b,
            const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a)

  {
    return T2s_times_T2_01(a, b, Number<Dim>(), Number<Dim>());
  }

  /* Single contraction.  The wrapper class has a different name for
     each possible placing of the indices (e.g. A(i,j)*B(j,k) has the
     number 10 because the contraction indices are on the second and
     first slots (counting from 0). */

  /* A(i,j)*B(j,k) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  class Tensor2_symmetric_times_Tensor2_10
  {
    Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    Tensor2_Expr<B, U, Dim, Dim1, j, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1) * iterB(Current_Dim - 1, N2)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(N1, 0) * iterB(0, N2);
    }

  public:
    Tensor2_symmetric_times_Tensor2_10(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
      const Tensor2_Expr<B, U, Dim, Dim1, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_10<A, B, T, U, Dim, Dim1, i, j, k>,
    typename promote<T, U>::V, Dim, Dim1, i, k>
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
            const Tensor2_Expr<B, U, Dim, Dim1, j, k> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_10<A, B, T, U, Dim, Dim1, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim1, i, k>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(i,j) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_10<A, B, T, U, Dim, Dim1, i, j, k>,
    typename promote<T, U>::V, Dim, Dim1, i, k>
  operator*(const Tensor2_Expr<B, U, Dim, Dim1, j, k> &b,
            const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_10<A, B, T, U, Dim, Dim1, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim1, i, k>(
      TensorExpr(a, b));
  }

  /* A(i,j)*B(k,j) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  class Tensor2_symmetric_times_Tensor2_11
  {
    Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    Tensor2_Expr<B, U, Dim1, Dim, k, j> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1) * iterB(N2, Current_Dim - 1)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(N1, 0) * iterB(N2, 0);
    }

  public:
    Tensor2_symmetric_times_Tensor2_11(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
      const Tensor2_Expr<B, U, Dim1, Dim, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_11<A, B, T, U, Dim, Dim1, i, j, k>,
    typename promote<T, U>::V, Dim, Dim1, i, k>
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
            const Tensor2_Expr<B, U, Dim1, Dim, k, j> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_11<A, B, T, U, Dim, Dim1, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim1, i, k>(
      TensorExpr(a, b));
  }

  /* B(k,j)*A(i,j) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_11<A, B, T, U, Dim, Dim1, i, j, k>,
    typename promote<T, U>::V, Dim, Dim1, i, k>
  operator*(const Tensor2_Expr<B, U, Dim1, Dim, k, j> &b,
            const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_11<A, B, T, U, Dim, Dim1, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim1, i, k>(
      TensorExpr(a, b));
  }

  /* A(j,i)*B(j,k) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  class Tensor2_symmetric_times_Tensor2_00
  {
    Tensor2_symmetric_Expr<A, T, Dim, j, i> iterA;
    Tensor2_Expr<B, U, Dim, Dim1, j, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(Current_Dim - 1, N1) * iterB(Current_Dim - 1, N2)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(0, N1) * iterB(0, N2);
    }

  public:
    Tensor2_symmetric_times_Tensor2_00(
      const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a,
      const Tensor2_Expr<B, U, Dim, Dim1, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_00<A, B, T, U, Dim, Dim1, i, j, k>,
    typename promote<T, U>::V, Dim, Dim1, i, k>
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a,
            const Tensor2_Expr<B, U, Dim, Dim1, j, k> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_00<A, B, T, U, Dim, Dim1, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim1, i, k>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(j,i) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_00<A, B, T, U, Dim, Dim1, i, j, k>,
    typename promote<T, U>::V, Dim, Dim1, i, k>
  operator*(const Tensor2_Expr<B, U, Dim, Dim1, j, k> &b,
            const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_00<A, B, T, U, Dim, Dim1, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim1, i, k>(
      TensorExpr(a, b));
  }

  /* A(j,i)*B(k,j) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  class Tensor2_symmetric_times_Tensor2_01
  {
    Tensor2_symmetric_Expr<A, T, Dim, j, i> iterA;
    Tensor2_Expr<B, U, Dim1, Dim, k, j> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(Current_Dim - 1, N1) * iterB(N2, Current_Dim - 1)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(0, N1) * iterB(N2, 0);
    }

  public:
    Tensor2_symmetric_times_Tensor2_01(
      const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a,
      const Tensor2_Expr<B, U, Dim1, Dim, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_01<A, B, T, U, Dim, Dim1, i, j, k>,
    typename promote<T, U>::V, Dim, Dim1, i, k>
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a,
            const Tensor2_Expr<B, U, Dim1, Dim, k, j> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_01<A, B, T, U, Dim, Dim1, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim1, i, k>(
      TensorExpr(a, b));
  }

  /* B(k,j)*A(j,i) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_01<A, B, T, U, Dim, Dim1, i, j, k>,
    typename promote<T, U>::V, Dim, Dim1, i, k>
  operator*(const Tensor2_Expr<B, U, Dim1, Dim, k, j> &b,
            const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_01<A, B, T, U, Dim, Dim1, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim1, i, k>(
      TensorExpr(a, b));
  }
}
