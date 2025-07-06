/* This file has all of the declarations for expressions like
   Tensor4*Tensor2_symmetric and Tensor2_symmetric*Tensor4, yielding a
   Tensor2. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(k,l) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim,
            char i, char j, char k, char l>
  class Tensor4_times_Tensor2_symmetric_23
  {
    Tensor4_Expr<A, T, Dim0, Dim1, Dim, Dim, i, j, k, l> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, k, l> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, 0, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, N2, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, N2, 0, 0) * iterB(0, 0);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>(), Number<Dim>());
    }

    Tensor4_times_Tensor2_symmetric_23(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim, Dim, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, k, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim,
            char i, char j, char k, char l>
  Tensor2_Expr<
    Tensor4_times_Tensor2_symmetric_23<A, B, T, U, Dim0, Dim1, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim1, i, j>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim, Dim, i, j, k, l> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, k, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_symmetric_23<A, B, T, U, Dim0, Dim1, Dim, i, j,
                                           k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, b));
  }

  /* B(k,l)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim,
            char i, char j, char k, char l>
  Tensor2_Expr<
    Tensor4_times_Tensor2_symmetric_23<A, B, T, U, Dim0, Dim1, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim1, i, j>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, k, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim, Dim, i, j, k, l> &a)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_symmetric_23<A, B, T, U, Dim0, Dim1, Dim, i, j,
                                           k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(l,k) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim,
            char i, char j, char k, char l>
  class Tensor4_times_Tensor2_symmetric_32
  {
    Tensor4_Expr<A, T, Dim0, Dim1, Dim, Dim, i, j, k, l> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, l, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(N1, N2, 0, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, N2, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, N2, 0, 0) * iterB(0, 0);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>(), Number<Dim>());
    }

    Tensor4_times_Tensor2_symmetric_32(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim, Dim, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, l, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim,
            char i, char j, char k, char l>
  Tensor2_Expr<
    Tensor4_times_Tensor2_symmetric_32<A, B, T, U, Dim0, Dim1, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim1, i, j>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim, Dim, i, j, k, l> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, l, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_symmetric_32<A, B, T, U, Dim0, Dim1, Dim, i, j,
                                           k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, b));
  }

  /* B(l,k)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim,
            char i, char j, char k, char l>
  Tensor2_Expr<
    Tensor4_times_Tensor2_symmetric_32<A, B, T, U, Dim0, Dim1, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim1, i, j>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, l, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim, Dim, i, j, k, l> &a)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_symmetric_32<A, B, T, U, Dim0, Dim1, Dim, i, j,
                                           k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(i,l) */

  template <class A, class B, class T, class U, int Dim, int Dim1, int Dim2,
            char i, char j, char k, char l>
  class Tensor4_times_Tensor2_symmetric_03
  {
    Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, i, l> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, N1, N2, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(0, N1, N2, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, N2, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(0, N1, N2, 0) * iterB(0, 0);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>(), Number<Dim>());
    }

    Tensor4_times_Tensor2_symmetric_03(
      const Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, i, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, int Dim1, int Dim2,
            char i, char j, char k, char l>
  Tensor2_Expr<
    Tensor4_times_Tensor2_symmetric_03<A, B, T, U, Dim1, Dim2, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim1, Dim2, j, k>
  operator*(const Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, i, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_symmetric_03<A, B, T, U, Dim1, Dim2, Dim, i, j,
                                           k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        k>(TensorExpr(a, b));
  }

  /* B(i,l)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim, int Dim1, int Dim2,
            char i, char j, char k, char l>
  Tensor2_Expr<
    Tensor4_times_Tensor2_symmetric_03<A, B, T, U, Dim1, Dim2, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim1, Dim2, j, k>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, i, l> &b,
            const Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> &a)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_symmetric_03<A, B, T, U, Dim1, Dim2, Dim, i, j,
                                           k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        k>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(l,i) */

  template <class A, class B, class T, class U, int Dim, int Dim1, int Dim2,
            char i, char j, char k, char l>
  class Tensor4_times_Tensor2_symmetric_30
  {
    Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, l, i> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, N1, N2, Current_Dim1 - 1)
               * iterB(Current_Dim1 - 1, Current_Dim0 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(0, N1, N2, Current_Dim1 - 1) * iterB(Current_Dim1 - 1, 0)
             + eval(N1, N2, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(0, N1, N2, 0) * iterB(0, 0);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>(), Number<Dim>());
    }

    Tensor4_times_Tensor2_symmetric_30(
      const Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, l, i> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, int Dim1, int Dim2,
            char i, char j, char k, char l>
  Tensor2_Expr<
    Tensor4_times_Tensor2_symmetric_30<A, B, T, U, Dim1, Dim2, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim1, Dim2, j, k>
  operator*(const Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, l, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_symmetric_30<A, B, T, U, Dim1, Dim2, Dim, i, j,
                                           k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        k>(TensorExpr(a, b));
  }

  /* B(l,i)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim, int Dim1, int Dim2,
            char i, char j, char k, char l>
  Tensor2_Expr<
    Tensor4_times_Tensor2_symmetric_30<A, B, T, U, Dim1, Dim2, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim1, Dim2, j, k>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, l, i> &b,
            const Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> &a)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_symmetric_30<A, B, T, U, Dim1, Dim2, Dim, i, j,
                                           k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        k>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(i,k) */

  template <class A, class B, class T, class U, int Dim, int Dim1, int Dim2,
            char i, char j, char k, char l>
  class Tensor4_times_Tensor2_symmetric_02
  {
    Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, i, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim0> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, N1, Current_Dim1 - 1, N2)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, N2, Number<Current_Dim0 - 1>(),
                    Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &,
         const Number<Current_Dim1> &) const
    {
      return iterA(0, N1, Current_Dim1 - 1, N2) * iterB(0, Current_Dim1 - 1)
             + eval(N1, N2, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(0, N1, 0, N2) * iterB(0, 0);
    }

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>(), Number<Dim>());
    }

    Tensor4_times_Tensor2_symmetric_02(
      const Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, i, k> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim, int Dim1, int Dim2,
            char i, char j, char k, char l>
  Tensor2_Expr<
    Tensor4_times_Tensor2_symmetric_02<A, B, T, U, Dim1, Dim2, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim1, Dim2, j, l>
  operator*(const Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, i, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_symmetric_02<A, B, T, U, Dim1, Dim2, Dim, i, j,
                                           k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        l>(TensorExpr(a, b));
  }

  /* B(i,k)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim, int Dim1, int Dim2,
            char i, char j, char k, char l>
  Tensor2_Expr<
    Tensor4_times_Tensor2_symmetric_02<A, B, T, U, Dim1, Dim2, Dim, i, j, k, l>,
    typename promote<T, U>::V, Dim1, Dim2, j, l>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, i, k> &b,
            const Tensor4_Expr<A, T, Dim, Dim1, Dim2, Dim, i, j, k, l> &a)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_symmetric_02<A, B, T, U, Dim1, Dim2, Dim, i, j,
                                           k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        l>(TensorExpr(a, b));
  }
}
