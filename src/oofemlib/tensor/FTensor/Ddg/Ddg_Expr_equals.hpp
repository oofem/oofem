/* Various assignment operators.  I have to explicitly declare the
   second operator= because otherwise the compiler will generate its
   own and not use the template code. */

#pragma once

namespace FTensor
{
  template <class A, class B, class U, int Current_Dim0, int Current_Dim1,
            int Current_Dim2, int Current_Dim3, int Dim01, int Dim23, char i,
            char j, char k, char l>
  void
  T4ddg_equals_T4ddg(A &iter,
                     const Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> &result,
                     const Number<Current_Dim0> &,
                     const Number<Current_Dim1> &,
                     const Number<Current_Dim2> &,
                     const Number<Current_Dim3> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1,
         Current_Dim3 - 1)
      = result(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1,
               Current_Dim3 - 1);
    T4ddg_equals_T4ddg(iter, result, Number<Current_Dim0 - 1>(),
                       Number<Current_Dim1>(), Number<Current_Dim2>(),
                       Number<Current_Dim3>());
  }

  template <class A, class B, class U, int Current_Dim1, int Current_Dim2,
            int Current_Dim3, int Dim01, int Dim23, char i, char j, char k,
            char l>
  void
  T4ddg_equals_T4ddg(A &iter,
                     const Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> &result,
                     const Number<1> &, const Number<Current_Dim1> &,
                     const Number<Current_Dim2> &,
                     const Number<Current_Dim3> &)
  {
    iter(0, Current_Dim1 - 1, Current_Dim2 - 1, Current_Dim3 - 1)
      = result(0, Current_Dim1 - 1, Current_Dim2 - 1, Current_Dim3 - 1);
    T4ddg_equals_T4ddg(iter, result, Number<Current_Dim1 - 1>(),
                       Number<Current_Dim1 - 1>(), Number<Current_Dim2>(),
                       Number<Current_Dim3>());
  }

  template <class A, class B, class U, int Current_Dim2, int Current_Dim3,
            int Dim01, int Dim23, char i, char j, char k, char l>
  void
  T4ddg_equals_T4ddg(A &iter,
                     const Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> &result,
                     const Number<1> &, const Number<1> &,
                     const Number<Current_Dim2> &,
                     const Number<Current_Dim3> &)
  {
    iter(0, 0, Current_Dim2 - 1, Current_Dim3 - 1)
      = result(0, 0, Current_Dim2 - 1, Current_Dim3 - 1);
    T4ddg_equals_T4ddg(iter, result, Number<Dim01>(), Number<Dim01>(),
                       Number<Current_Dim2 - 1>(), Number<Current_Dim3>());
  }

  template <class A, class B, class U, int Current_Dim3, int Dim01, int Dim23,
            char i, char j, char k, char l>
  void
  T4ddg_equals_T4ddg(A &iter,
                     const Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> &result,
                     const Number<1> &, const Number<1> &, const Number<1> &,
                     const Number<Current_Dim3> &)
  {
    iter(0, 0, 0, Current_Dim3 - 1) = result(0, 0, 0, Current_Dim3 - 1);
    T4ddg_equals_T4ddg(iter, result, Number<Dim01>(), Number<Dim01>(),
                       Number<Current_Dim3 - 1>(), Number<Current_Dim3 - 1>());
  }

  template <class A, class B, class U, int Dim01, int Dim23, char i, char j,
            char k, char l>
  void
  T4ddg_equals_T4ddg(A &iter,
                     const Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> &result,
                     const Number<1> &, const Number<1> &, const Number<1> &,
                     const Number<1> &)
  {
    iter(0, 0, 0, 0) = result(0, 0, 0, 0);
  }

  template <class A, class T, int Tensor_Dim01, int Tensor_Dim23, int Dim01,
            int Dim23, char i, char j, char k, char l>
  template <class B, class U>
  Ddg_Expr<Ddg<A, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i, j, k, l> &
  Ddg_Expr<Ddg<A, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i, j, k, l>::
  operator=(const Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> &result)
  {
    T4ddg_equals_T4ddg(iter, result, Number<Dim01>(), Number<Dim01>(),
                       Number<Dim23>(), Number<Dim23>());
    return *this;
  }

  template <class A, class T, int Tensor_Dim01, int Tensor_Dim23, int Dim01,
            int Dim23, char i, char j, char k, char l>
  Ddg_Expr<Ddg<A, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i, j, k, l> &
  Ddg_Expr<Ddg<A, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i, j, k, l>::
  operator=(const Ddg_Expr<Ddg<A, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23,
                           i, j, k, l> &result)
  {
    return operator=<Ddg<A, Tensor_Dim01, Tensor_Dim23>, T>(result);
  }
}
