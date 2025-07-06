/* Fully contract a Ddg with a Ddg. */

#pragma once

namespace FTensor
{
  /* A(i,j,k,l)*B(i,k,j,l) */

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l, int Current_Dim0, int Current_Dim1,
            int Current_Dim2, int Current_Dim3>
  typename promote<T, U>::V
  T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                         const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                         const Number<Current_Dim0> &,
                         const Number<Current_Dim1> &,
                         const Number<Current_Dim2> &,
                         const Number<Current_Dim3> &)
  {
    return a(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1,
             Current_Dim3 - 1)
             * b(Current_Dim0 - 1, Current_Dim2 - 1, Current_Dim1 - 1,
                 Current_Dim3 - 1)
           + T4ddg_times_T4ddg_0213(
               a, b, Number<Current_Dim0 - 1>(), Number<Current_Dim1>(),
               Number<Current_Dim2>(), Number<Current_Dim3>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j, char k,
            char l, int Current_Dim1, int Current_Dim2, int Current_Dim3>
  typename promote<T, U>::V
  T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                         const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                         const Number<1> &, const Number<Current_Dim1> &,
                         const Number<Current_Dim2> &,
                         const Number<Current_Dim3> &)
  {
    return a(0, Current_Dim1 - 1, Current_Dim2 - 1, Current_Dim3 - 1)
             * b(0, Current_Dim2 - 1, Current_Dim1 - 1, Current_Dim3 - 1)
           + T4ddg_times_T4ddg_0213(
               a, b, Number<Dim>(), Number<Current_Dim1 - 1>(),
               Number<Current_Dim2>(), Number<Current_Dim3>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l, int Current_Dim2, int Current_Dim3>
  typename promote<T, U>::V
  T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                         const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                         const Number<1> &, const Number<1> &,
                         const Number<Current_Dim2> &,
                         const Number<Current_Dim3> &)
  {
    return a(0, 0, Current_Dim2 - 1, Current_Dim3 - 1)
             * b(0, Current_Dim2 - 1, 0, Current_Dim3 - 1)
           + T4ddg_times_T4ddg_0213(a, b, Number<Dim>(), Number<Dim>(),
                                    Number<Current_Dim2 - 1>(),
                                    Number<Current_Dim3>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l, int Current_Dim3>
  typename promote<T, U>::V
  T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                         const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                         const Number<1> &, const Number<1> &,
                         const Number<1> &, const Number<Current_Dim3> &)
  {
    return a(0, 0, 0, Current_Dim3 - 1) * b(0, 0, 0, Current_Dim3 - 1)
           + T4ddg_times_T4ddg_0213(a, b, Number<Dim>(), Number<Dim>(),
                                    Number<Dim>(), Number<Current_Dim3 - 1>());
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  typename promote<T, U>::V
  T4ddg_times_T4ddg_0213(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
                         const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b,
                         const Number<1> &, const Number<1> &,
                         const Number<1> &, const Number<1> &)
  {
    return a(0, 0, 0, 0) * b(0, 0, 0, 0);
  }

  template <class A, class B, class T, class U, int Dim, char i, char j,
            char k, char l>
  typename promote<T, U>::V
  operator*(const Ddg_Expr<A, T, Dim, Dim, i, j, k, l> &a,
            const Ddg_Expr<B, U, Dim, Dim, i, k, j, l> &b)
  {
    return T4ddg_times_T4ddg_0213(a, b, Number<Dim>(), Number<Dim>(),
                                  Number<Dim>(), Number<Dim>());
  }
}
