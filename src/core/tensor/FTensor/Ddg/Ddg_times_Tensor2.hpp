/* This file has all of the declarations for expressions like
   Ddg*Tensor2 and Tensor2*Ddg, yielding a
   Tensor2_symmetric. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_times_Tensor2_23
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor2_Expr<B, U, Dim23, Dim23, k, l> iterB;

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
             + eval(N1, N2, Number<Dim23>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, N2, 0, 0) * iterB(0, 0);
    }

  public:
    Ddg_times_Tensor2_23(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                         const Tensor2_Expr<B, U, Dim23, Dim23, k, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim23>(), Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
    Ddg_times_Tensor2_23<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, i, j>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim23, Dim23, k, l> &b)
  {
    using TensorExpr
      = Ddg_times_Tensor2_23<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim01,
                                  i, j>(TensorExpr(a, b));
  }

  /* B(k,l)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
    Ddg_times_Tensor2_23<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, i, j>
  operator*(const Tensor2_Expr<B, U, Dim23, Dim23, k, l> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    using TensorExpr
      = Ddg_times_Tensor2_23<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim01,
                                  i, j>(TensorExpr(a, b));
  }

  /* A(i,j,k,l)*B(l,k) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_times_Tensor2_32
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor2_Expr<B, U, Dim23, Dim23, l, k> iterB;

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
             + eval(N1, N2, Number<Dim23>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V eval(const int N1, const int N2,
                                   const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, N2, 0, 0) * iterB(0, 0);
    }

  public:
    Ddg_times_Tensor2_32(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                         const Tensor2_Expr<B, U, Dim23, Dim23, l, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim23>(), Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
    Ddg_times_Tensor2_32<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, i, j>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim23, Dim23, l, k> &b)
  {
    using TensorExpr
      = Ddg_times_Tensor2_32<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim01,
                                  i, j>(TensorExpr(a, b));
  }

  /* B(l,k)*A(i,j,k,l) */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Tensor2_symmetric_Expr<
    Ddg_times_Tensor2_32<A, B, T, U, Dim01, Dim23, i, j, k, l>,
    typename promote<T, U>::V, Dim01, i, j>
  operator*(const Tensor2_Expr<B, U, Dim23, Dim23, l, k> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    using TensorExpr
      = Ddg_times_Tensor2_32<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim01,
                                  i, j>(TensorExpr(a, b));
  }
}
