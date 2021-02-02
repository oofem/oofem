/* This file has all of the declarations for
   Tensor2_symmetric*Tensor2_symmetric.  This includes the double
   contraction A(i,j)*B(i,j) (yielding a typename promote<T,U>::V) as well as
   the more complicated single contraction A(i,j)*B(j,k) (yielding a Tensor2
   expression), and no contractions A(i,j)*B(k,l) -> Ddg */

#pragma once

namespace FTensor
{
  /* Double contraction. */

  /* A(i,j)*B(i,j) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            int Current_Dim0, int Current_Dim1>
  typename promote<T, U>::V
  T2s_times_T2s(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                const Tensor2_symmetric_Expr<B, U, Dim, i, j> &b,
                const Number<Current_Dim0> &, const Number<Current_Dim1> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1)
             * b(Current_Dim0 - 1, Current_Dim1 - 1)
           + T2s_times_T2s(a, b, Number<Current_Dim0 - 1>(),
                           Number<Current_Dim1>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            int Current_Dim1>
  typename promote<T, U>::V
  T2s_times_T2s(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                const Tensor2_symmetric_Expr<B, U, Dim, i, j> &b,
                const Number<1> &, const Number<Current_Dim1> &)
  {
    return a(0, Current_Dim1 - 1) * b(0, Current_Dim1 - 1)
           + T2s_times_T2s(a, b, Number<Dim>(), Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j>
  typename promote<T, U>::V
  T2s_times_T2s(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                const Tensor2_symmetric_Expr<B, U, Dim, i, j> &b,
                const Number<1> &, const Number<1> &)
  {
    return a(0, 0) * b(0, 0);
  }

  template <class A, class B, class T, class U, int Dim, char i, char j>
  typename promote<T, U>::V
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, i, j> &b)
  {
    return T2s_times_T2s(a, b, Number<Dim>(), Number<Dim>());
  }

  /* A(i,j)*B(j,i) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            int Current_Dim0, int Current_Dim1>
  typename promote<T, U>::V
  T2s_times_switched_T2s(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                         const Tensor2_symmetric_Expr<B, U, Dim, j, i> &b,
                         const Number<Current_Dim0> &,
                         const Number<Current_Dim1> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1)
             * b(Current_Dim1 - 1, Current_Dim0 - 1)
           + T2s_times_switched_T2s(a, b, Number<Current_Dim0 - 1>(),
                                    Number<Current_Dim1>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            int Current_Dim1>
  typename promote<T, U>::V
  T2s_times_switched_T2s(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                         const Tensor2_symmetric_Expr<B, U, Dim, j, i> &b,
                         const Number<1> &, const Number<Current_Dim1> &)
  {
    return a(0, Current_Dim1 - 1) * b(Current_Dim1 - 1, 0)
           + T2s_times_switched_T2s(a, b, Number<Dim>(),
                                    Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j>
  typename promote<T, U>::V
  T2s_times_switched_T2s(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
                         const Tensor2_symmetric_Expr<B, U, Dim, j, i> &b,
                         const Number<1> &, const Number<1> &)
  {
    return a(0, 0) * b(0, 0);
  }

  template <class A, class B, class T, class U, int Dim, char i, char j>
  typename promote<T, U>::V
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, j, i> &b)
  {
    return T2s_times_switched_T2s(a, b, Number<Dim>(), Number<Dim>());
  }

  /* Single contraction.  The wrapper class has a different name for
     each possible placing of the indices (e.g. A(i,j)*B(j,k) has the
     number 10 because the contraction indices are on the second and
     first slots (counting from 0). */

  /* A(i,j)*B(j,k) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Tensor2_symmetric_times_Tensor2_symmetric_10
  {
    Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, j, k> iterB;

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
    Tensor2_symmetric_times_Tensor2_symmetric_10(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_symmetric_10<A, B, T, U, Dim, i, j, k>,
    typename promote<T, U>::V, Dim, Dim, i, k>
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_symmetric_10<A, B, T, U, Dim, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim, i, k>(
      TensorExpr(a, b));
  }

  /* A(i,j)*B(k,j) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Tensor2_symmetric_times_Tensor2_symmetric_11
  {
    Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, k, j> iterB;

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
    Tensor2_symmetric_times_Tensor2_symmetric_11(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_symmetric_11<A, B, T, U, Dim, i, j, k>,
    typename promote<T, U>::V, Dim, Dim, i, k>
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, k, j> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_symmetric_11<A, B, T, U, Dim, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim, i, k>(
      TensorExpr(a, b));
  }

  /* A(j,i)*B(j,k) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Tensor2_symmetric_times_Tensor2_symmetric_00
  {
    Tensor2_symmetric_Expr<A, T, Dim, j, i> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, j, k> iterB;

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
    Tensor2_symmetric_times_Tensor2_symmetric_00(
      const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_symmetric_00<A, B, T, U, Dim, i, j, k>,
    typename promote<T, U>::V, Dim, Dim, i, k>
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_symmetric_00<A, B, T, U, Dim, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim, i, k>(
      TensorExpr(a, b));
  }

  /* A(j,i)*B(k,j) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Tensor2_symmetric_times_Tensor2_symmetric_01
  {
    Tensor2_symmetric_Expr<A, T, Dim, j, i> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, k, j> iterB;

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
    Tensor2_symmetric_times_Tensor2_symmetric_01(
      const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor2_Expr<
    Tensor2_symmetric_times_Tensor2_symmetric_01<A, B, T, U, Dim, i, j, k>,
    typename promote<T, U>::V, Dim, Dim, i, k>
  operator*(const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, k, j> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_symmetric_01<A, B, T, U, Dim, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim, Dim, i, k>(
      TensorExpr(a, b));
  }

  /* A(i,j)*B(k,l) -> Ddg */

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j, char k, char l>
  class Tensor2_symmetric_times_Tensor2_symmetric
  {
    Tensor2_symmetric_Expr<A, T, Dim0, i, j> iterA;
    Tensor2_symmetric_Expr<B, U, Dim1, k, l> iterB;

  public:
    Tensor2_symmetric_times_Tensor2_symmetric(
      const Tensor2_symmetric_Expr<A, T, Dim0, i, j> &a,
      const Tensor2_symmetric_Expr<B, U, Dim1, k, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2) * iterB(N3, N4);
    }
  };
  /*
  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j, char k, char l>
  Ddg_Expr<Tensor2_symmetric_times_Tensor2_symmetric<A, B, T, U, Dim0, Dim1, i,
                                                     j, k, l>,
           typename promote<T, U>::V, Dim0, Dim1, i, j, k, l>
  operator*(const Tensor2_symmetric_Expr<A, T, Dim0, i, j> &a,
            const Tensor2_symmetric_Expr<B, U, Dim1, k, l> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_times_Tensor2_symmetric<A, B, T, U, Dim0, Dim1, i, j,
                                                  k, l>;
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i, j, k,
                    l>(TensorExpr(a, b));
		    }*/
}
