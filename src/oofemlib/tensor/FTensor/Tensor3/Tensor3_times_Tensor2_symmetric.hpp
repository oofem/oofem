/* This file has all of the declarations for expressions like
   Tensor3*Tensor2_symmetric and Tensor2_symmetric*Tensor3, yielding a
   Tensor1. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)*B(j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k>
  class Tensor3_times_Tensor2_symmetric_12
  {
    Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, j, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(N1, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(N1, 0, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, 0, 0) * iterB(0, 0);
    }

  public:
    Tensor3_times_Tensor2_symmetric_12(
      const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim>(), Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k>
  Tensor1_Expr<
    Tensor3_times_Tensor2_symmetric_12<A, B, T, U, Dim0, Dim, i, j, k>,
    typename promote<T, U>::V, Dim0, i>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_symmetric_12<A, B, T, U, Dim0, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(i,j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim0, int Dim, char i,
            char j, char k>
  Tensor1_Expr<
    Tensor3_times_Tensor2_symmetric_12<A, B, T, U, Dim0, Dim, i, j, k>,
    typename promote<T, U>::V, Dim0, i>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b,
            const Tensor3_Expr<A, T, Dim0, Dim, Dim, i, j, k> &a)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_symmetric_12<A, B, T, U, Dim0, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  /* A(j,i,k)*B(j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim1, int Dim, char i,
            char j, char k>
  class Tensor3_times_Tensor2_symmetric_02
  {
    Tensor3_Expr<A, T, Dim, Dim1, Dim, j, i, k> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, j, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, N1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(0, N1, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(0, N1, 0) * iterB(0, 0);
    }

  public:
    Tensor3_times_Tensor2_symmetric_02(
      const Tensor3_Expr<A, T, Dim, Dim1, Dim, j, i, k> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim>(), Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim1, int Dim, char i,
            char j, char k>
  Tensor1_Expr<
    Tensor3_times_Tensor2_symmetric_02<A, B, T, U, Dim1, Dim, i, j, k>,
    typename promote<T, U>::V, Dim1, i>
  operator*(const Tensor3_Expr<A, T, Dim, Dim1, Dim, j, i, k> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_symmetric_02<A, B, T, U, Dim1, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, i>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(j,i,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim1, int Dim, char i,
            char j, char k>
  Tensor1_Expr<
    Tensor3_times_Tensor2_symmetric_02<A, B, T, U, Dim1, Dim, i, j, k>,
    typename promote<T, U>::V, Dim1, i>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b,
            const Tensor3_Expr<A, T, Dim, Dim1, Dim, j, i, k> &a)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_symmetric_02<A, B, T, U, Dim1, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, i>(
      TensorExpr(a, b));
  }

  /* A(j,k,i)*B(j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim2, int Dim, char i,
            char j, char k>
  class Tensor3_times_Tensor2_symmetric_01
  {
    Tensor3_Expr<A, T, Dim, Dim, Dim2, j, k, i> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, j, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, Current_Dim1 - 1, N1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(0, Current_Dim1 - 1, N1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(0, 0, N1) * iterB(0, 0);
    }

  public:
    Tensor3_times_Tensor2_symmetric_01(
      const Tensor3_Expr<A, T, Dim, Dim, Dim2, j, k, i> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim>(), Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim2, int Dim, char i,
            char j, char k>
  Tensor1_Expr<
    Tensor3_times_Tensor2_symmetric_01<A, B, T, U, Dim2, Dim, i, j, k>,
    typename promote<T, U>::V, Dim2, i>
  operator*(const Tensor3_Expr<A, T, Dim, Dim, Dim2, j, k, i> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_symmetric_01<A, B, T, U, Dim2, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, i>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(j,k,i)->Tensor1 */

  template <class A, class B, class T, class U, int Dim2, int Dim, char i,
            char j, char k>
  Tensor1_Expr<
    Tensor3_times_Tensor2_symmetric_01<A, B, T, U, Dim2, Dim, i, j, k>,
    typename promote<T, U>::V, Dim2, i>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b,
            const Tensor3_Expr<A, T, Dim, Dim, Dim2, j, k, i> &a)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_symmetric_01<A, B, T, U, Dim2, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, i>(
      TensorExpr(a, b));
  }
}
