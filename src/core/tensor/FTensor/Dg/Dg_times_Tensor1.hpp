/* This file has all of the declarations for expressions like
   Dg*Tensor1 and Tensor1*Dg, yielding a
   Tensor2_symmetric or Tensor2. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)*B(k)->Tensor2_symmetric */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  class Dg_times_Tensor1_2
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Tensor1_Expr<B, U, Dim2, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, Current_Dim - 1) * iterB(Current_Dim - 1)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(0);
    }

  public:
    Dg_times_Tensor1_2(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
                       const Tensor1_Expr<B, U, Dim2, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Tensor2_symmetric_Expr<Dg_times_Tensor1_2<A, B, T, U, Dim01, Dim2, i, j, k>,
                         typename promote<T, U>::V, Dim01, i, j>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Tensor1_Expr<B, U, Dim2, k> &b)
  {
    using TensorExpr = Dg_times_Tensor1_2<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim01,
                                  i, j>(TensorExpr(a, b));
  }

  /* B(k)*A(i,j,k)->Tensor2_symmetric */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Tensor2_symmetric_Expr<Dg_times_Tensor1_2<A, B, T, U, Dim01, Dim2, i, j, k>,
                         typename promote<T, U>::V, Dim01, i, j>
  operator*(const Tensor1_Expr<B, U, Dim2, k> &b,
            const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a)
  {
    using TensorExpr = Dg_times_Tensor1_2<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim01,
                                  i, j>(TensorExpr(a, b));
  }

  /* A(i,k,j)*B(k)->Tensor2 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  class Dg_times_Tensor1_1
  {
    Dg_Expr<A, T, Dim01, Dim2, i, k, j> iterA;
    Tensor1_Expr<B, U, Dim01, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1, N2) * iterB(Current_Dim - 1)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(N1, 0, N2) * iterB(0);
    }

  public:
    Dg_times_Tensor1_1(const Dg_Expr<A, T, Dim01, Dim2, i, k, j> &a,
                       const Tensor1_Expr<B, U, Dim01, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Tensor2_Expr<Dg_times_Tensor1_1<A, B, T, U, Dim01, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim01, Dim2, i, j>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, i, k, j> &a,
            const Tensor1_Expr<B, U, Dim01, k> &b)
  {
    using TensorExpr = Dg_times_Tensor1_1<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i,
                        j>(TensorExpr(a, b));
  }

  /* B(k)*A(i,k,j)->Tensor2 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Tensor2_Expr<Dg_times_Tensor1_1<A, B, T, U, Dim01, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim01, Dim2, i, j>
  operator*(const Tensor1_Expr<B, U, Dim01, k> &b,
            const Dg_Expr<A, T, Dim01, Dim2, i, k, j> &a)
  {
    using TensorExpr = Dg_times_Tensor1_1<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i,
                        j>(TensorExpr(a, b));
  }

  /* A(k,i,j)*B(k)->Tensor2 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  class Dg_times_Tensor1_0
  {
    Dg_Expr<A, T, Dim01, Dim2, k, i, j> iterA;
    Tensor1_Expr<B, U, Dim01, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1, N2) * iterB(Current_Dim - 1)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(N1, 0, N2) * iterB(0);
    }

  public:
    Dg_times_Tensor1_0(const Dg_Expr<A, T, Dim01, Dim2, k, i, j> &a,
                       const Tensor1_Expr<B, U, Dim01, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Tensor2_Expr<Dg_times_Tensor1_0<A, B, T, U, Dim01, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim01, Dim2, i, j>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, k, i, j> &a,
            const Tensor1_Expr<B, U, Dim01, k> &b)
  {
    using TensorExpr = Dg_times_Tensor1_0<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i,
                        j>(TensorExpr(a, b));
  }

  /* B(k)*A(k,i,j)->Tensor2 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Tensor2_Expr<Dg_times_Tensor1_0<A, B, T, U, Dim01, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim01, Dim2, i, j>
  operator*(const Tensor1_Expr<B, U, Dim01, k> &b,
            const Dg_Expr<A, T, Dim01, Dim2, k, i, j> &a)
  {
    using TensorExpr = Dg_times_Tensor1_0<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i,
                        j>(TensorExpr(a, b));
  }
}
