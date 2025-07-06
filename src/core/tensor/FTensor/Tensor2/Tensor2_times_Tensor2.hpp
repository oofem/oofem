/* This file has all of the declarations for Tensor2*Tensor2.  This
   includes the double contraction A(i,j)*B(i,j) (yielding a
   typename promote<T,U>) as well as the more complicated single contraction
   A(i,j)*B(j,k) (yielding a Tensor2 expression) and no contraction
   (yielding a Tensor4). */

/* Double contraction. */

/* A(i,j)*B(i,j) */

#pragma once

namespace FTensor
{
  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j, int Current_Dim0, int Current_Dim1>
  typename promote<T, U>::V
  T2_times_T2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
              const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b,
              const Number<Current_Dim0> &, const Number<Current_Dim1> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1)
             * b(Current_Dim0 - 1, Current_Dim1 - 1)
           + T2_times_T2(a, b, Number<Current_Dim0 - 1>(),
                         Number<Current_Dim1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j, int Current_Dim1>
  typename promote<T, U>::V
  T2_times_T2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
              const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b, const Number<1> &,
              const Number<Current_Dim1> &)
  {
    return a(0, Current_Dim1 - 1) * b(0, Current_Dim1 - 1)
           + T2_times_T2(a, b, Number<Dim0>(), Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  typename promote<T, U>::V
  T2_times_T2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
              const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b, const Number<1> &,
              const Number<1> &)
  {
    return a(0, 0) * b(0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  typename promote<T, U>::V
  operator*(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
            const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b)
  {
    return T2_times_T2(a, b, Number<Dim0>(), Number<Dim1>());
  }

  /* A(i,j)*B(j,i) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j, int Current_Dim0, int Current_Dim1>
  typename promote<T, U>::V
  T2_times_switched_T2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                       const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b,
                       const Number<Current_Dim0> &,
                       const Number<Current_Dim1> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1)
             * b(Current_Dim1 - 1, Current_Dim0 - 1)
           + T2_times_switched_T2(a, b, Number<Current_Dim0 - 1>(),
                                  Number<Current_Dim1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j, int Current_Dim1>
  typename promote<T, U>::V
  T2_times_switched_T2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                       const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b,
                       const Number<1> &, const Number<Current_Dim1> &)
  {
    return a(0, Current_Dim1 - 1) * b(Current_Dim1 - 1, 0)
           + T2_times_switched_T2(a, b, Number<Dim0>(),
                                  Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  typename promote<T, U>::V
  T2_times_switched_T2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                       const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b,
                       const Number<1> &, const Number<1> &)
  {
    return a(0, 0) * b(0, 0);
  }

  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  typename promote<T, U>::V
  operator*(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
            const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b)
  {
    return T2_times_switched_T2(a, b, Number<Dim0>(), Number<Dim1>());
  }

  /* Single contraction.  The wrapper class has a different name for
     each possible placing of the indices (e.g. A(i,j)*B(j,k) has the
     number 10 because the contraction indices are on the second and
     first slots (counting from 0). */

  /* A(i,j)*B(j,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor2_times_Tensor2_10
  {
    const Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;
    const Tensor2_Expr<B, U, Dim1, Dim2, j, k> iterB;

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
    Tensor2_times_Tensor2_10(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                             const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  Tensor2_Expr<Tensor2_times_Tensor2_10<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim0, Dim2, i, k>
  operator*(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
            const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
  {
    using TensorExpr
      = Tensor2_times_Tensor2_10<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i,
                        k>(TensorExpr(a, b));
  }

  /* A(i,j)*B(k,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor2_times_Tensor2_11
  {
    Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;
    Tensor2_Expr<B, U, Dim2, Dim1, k, j> iterB;

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
    Tensor2_times_Tensor2_11(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                             const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  Tensor2_Expr<Tensor2_times_Tensor2_11<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim0, Dim2, i, k>
  operator*(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
            const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
  {
    using TensorExpr
      = Tensor2_times_Tensor2_11<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i,
                        k>(TensorExpr(a, b));
  }

  /* A(j,i)*B(j,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor2_times_Tensor2_00
  {
    Tensor2_Expr<A, T, Dim1, Dim0, j, i> iterA;
    Tensor2_Expr<B, U, Dim1, Dim2, j, k> iterB;

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
    Tensor2_times_Tensor2_00(const Tensor2_Expr<A, T, Dim1, Dim0, j, i> &a,
                             const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  Tensor2_Expr<Tensor2_times_Tensor2_00<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim0, Dim2, i, k>
  operator*(const Tensor2_Expr<A, T, Dim1, Dim0, j, i> &a,
            const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
  {
    using TensorExpr
      = Tensor2_times_Tensor2_00<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i,
                        k>(TensorExpr(a, b));
  }

  /* A(j,i)*B(k,j) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor2_times_Tensor2_01
  {
    Tensor2_Expr<A, T, Dim1, Dim0, j, i> iterA;
    Tensor2_Expr<B, U, Dim2, Dim1, k, j> iterB;

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
    Tensor2_times_Tensor2_01(const Tensor2_Expr<A, T, Dim1, Dim0, j, i> &a,
                             const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  Tensor2_Expr<Tensor2_times_Tensor2_01<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim0, Dim2, i, k>
  operator*(const Tensor2_Expr<A, T, Dim1, Dim0, j, i> &a,
            const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
  {
    using TensorExpr
      = Tensor2_times_Tensor2_01<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i,
                        k>(TensorExpr(a, b));
  }

  /* A(i,j)*B(k,l) -> Tensor4 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor2_times_Tensor2
  {
    Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;
    Tensor2_Expr<B, U, Dim2, Dim3, k, l> iterB;

  public:
    Tensor2_times_Tensor2(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                          const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2) * iterB(N3, N4);
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor4_Expr<
    Tensor2_times_Tensor2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>
  operator*(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
            const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b)
  {
    using TensorExpr
      = Tensor2_times_Tensor2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim3, i, j, k, l>(TensorExpr(a, b));
  }
}
