/* This file has all of the declarations for expressions like
   Ddg*Tensor1 and Tensor1*Ddg, yielding a
   Dg. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(k)->Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_times_Tensor1_2
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor1_Expr<B, U, Dim23, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, Current_Dim - 1, N3) * iterB(Current_Dim - 1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, N2, 0, N3) * iterB(0);
    }

  public:
    Ddg_times_Tensor1_2(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
                        const Tensor1_Expr<B, U, Dim23, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim23>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Dg_Expr<Ddg_times_Tensor1_2<A, B, T, U, Dim01, Dim23, i, j, k, l>,
          typename promote<T, U>::V, Dim01, Dim23, i, j, l>
  operator*(const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
            const Tensor1_Expr<B, U, Dim23, k> &b)
  {
    using TensorExpr
      = Ddg_times_Tensor1_2<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j,
                   l>(TensorExpr(a, b));
  }

  /* B(k)*A(i,j,k,l)->Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  Dg_Expr<Ddg_times_Tensor1_2<A, B, T, U, Dim01, Dim23, i, j, k, l>,
          typename promote<T, U>::V, Dim01, Dim23, i, j, l>
  operator*(const Tensor1_Expr<B, U, Dim23, k> &b,
            const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a)
  {
    using TensorExpr
      = Ddg_times_Tensor1_2<A, B, T, U, Dim01, Dim23, i, j, k, l>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim23, i, j,
                   l>(TensorExpr(a, b));
  }
}
