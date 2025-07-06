/* Various assignment operators, including specializations to
   Christof.  I have to explicitly declare the second
   operator= because otherwise the compiler will generate its own and
   not use the template code. */

#pragma once

namespace FTensor
{
  /* T3dg=T3dg */

  template <class A, class B, class U, int Dim01, int Dim2, char i, char j,
            char k, int Current_Dim0, int Current_Dim1, int Current_Dim2>
  void
  T3dg_equals_T3dg(A &iter, const Dg_Expr<B, U, Dim01, Dim2, i, j, k> &result,
                   const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                   const Number<Current_Dim2> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
      = result(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1);
    T3dg_equals_T3dg(iter, result, Number<Current_Dim0 - 1>(),
                     Number<Current_Dim1>(), Number<Current_Dim2>());
  }

  template <class A, class B, class U, int Dim01, int Dim2, char i, char j,
            char k, int Current_Dim1, int Current_Dim2>
  void
  T3dg_equals_T3dg(A &iter, const Dg_Expr<B, U, Dim01, Dim2, i, j, k> &result,
                   const Number<1> &, const Number<Current_Dim1> &,
                   const Number<Current_Dim2> &)
  {
    iter(0, Current_Dim1 - 1, Current_Dim2 - 1)
      = result(0, Current_Dim1 - 1, Current_Dim2 - 1);
    T3dg_equals_T3dg(iter, result, Number<Current_Dim1 - 1>(),
                     Number<Current_Dim1 - 1>(), Number<Current_Dim2>());
  }

  template <class A, class B, class U, int Dim01, int Dim2, char i, char j,
            char k, int Current_Dim2>
  void
  T3dg_equals_T3dg(A &iter, const Dg_Expr<B, U, Dim01, Dim2, i, j, k> &result,
                   const Number<1> &, const Number<1> &,
                   const Number<Current_Dim2> &)
  {
    iter(0, 0, Current_Dim2 - 1) = result(0, 0, Current_Dim2 - 1);
    T3dg_equals_T3dg(iter, result, Number<Dim01>(), Number<Dim01>(),
                     Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class U, int Dim01, int Dim2, char i, char j,
            char k>
  void
  T3dg_equals_T3dg(A &iter, const Dg_Expr<B, U, Dim01, Dim2, i, j, k> &result,
                   const Number<1> &, const Number<1> &, const Number<1> &)
  {
    iter(0, 0, 0) = result(0, 0, 0);
  }

  template <class A, class T, int Tensor_Dim01, int Tensor_Dim2, int Dim01,
            int Dim2, char i, char j, char k>
  template <class B, class U>
  Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k> &
  Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k>::
  operator=(const Dg_Expr<B, U, Dim01, Dim2, i, j, k> &result)
  {
    T3dg_equals_T3dg(iter, result, Number<Dim01>(), Number<Dim01>(),
                     Number<Dim2>());
    return *this;
  }

  /* T3dg=T3dg_Expr(T3dg) */

  template <class A, class T, int Tensor_Dim01, int Tensor_Dim2, int Dim01,
            int Dim2, char i, char j, char k>
  Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k> &
  Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k>::
  operator=(const Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i,
                          j, k> &result)
  {
    return operator=<Dg<A, Tensor_Dim01, Tensor_Dim2>, T>(result);
  }

  /* T3dg=U */

  template <class A, class U, int Dim01, int Dim2, int Current_Dim0,
            int Current_Dim1, int Current_Dim2>
  void T3dg_equals_generic(A &iter, const U &u, const Number<Current_Dim0> &,
                           const Number<Current_Dim1> &,
                           const Number<Current_Dim2> &, const Number<Dim01> &,
                           const Number<Dim2> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1) = u;
    T3dg_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                        Number<Current_Dim1>(), Number<Current_Dim2>(),
                        Number<Dim01>(), Number<Dim2>());
  }

  template <class A, class U, int Dim01, int Dim2, int Current_Dim1,
            int Current_Dim2>
  void T3dg_equals_generic(A &iter, const U &u, const Number<1> &,
                           const Number<Current_Dim1> &,
                           const Number<Current_Dim2> &, const Number<Dim01> &,
                           const Number<Dim2> &)
  {
    iter(0, Current_Dim1 - 1, Current_Dim2 - 1) = u;
    T3dg_equals_generic(iter, u, Number<Current_Dim1 - 1>(),
                        Number<Current_Dim1 - 1>(), Number<Current_Dim2>(),
                        Number<Dim01>(), Number<Dim2>());
  }

  template <class A, class U, int Dim01, int Dim2, int Current_Dim2>
  void T3dg_equals_generic(A &iter, const U &u, const Number<1> &,
                           const Number<1> &, const Number<Current_Dim2> &,
                           const Number<Dim01> &, const Number<Dim2> &)
  {
    iter(0, 0, Current_Dim2 - 1) = u;
    T3dg_equals_generic(iter, u, Number<Dim01>(), Number<Dim01>(),
                        Number<Current_Dim2 - 1>(), Number<Dim01>(),
                        Number<Dim2>());
  }

  template <class A, class U, int Dim01, int Dim2>
  void T3dg_equals_generic(A &iter, const U &u, const Number<1> &,
                           const Number<1> &, const Number<1> &,
                           const Number<Dim01> &, const Number<Dim2> &)
  {
    iter(0, 0, 0) = u;
  }

  template <class A, class T, int Tensor_Dim01, int Tensor_Dim2, int Dim01,
            int Dim2, char i, char j, char k>
  template <class U>
  Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k> &
  Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k>::
  operator=(const U &u)
  {
    T3dg_equals_generic(iter, u, Number<Dim01>(), Number<Dim01>(),
                        Number<Dim2>(), Number<Dim01>(), Number<Dim2>());
    return *this;
  }

  /* Assignment operators for the specializations to Christof. */

  /* T3ch=T3dg */

  template <class A, class B, class U, int Dim12, int Dim0, char i, char j,
            char k, int Current_Dim0, int Current_Dim1, int Current_Dim2>
  void
  T3ch_equals_T3dg(A &iter, const Dg_Expr<B, U, Dim12, Dim0, i, j, k> &result,
                   const Number<Current_Dim0> &, const Number<Current_Dim1> &,
                   const Number<Current_Dim2> &)
  {
    iter(Current_Dim2 - 1, Current_Dim0 - 1, Current_Dim1 - 1)
      = result(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1);
    T3ch_equals_T3dg(iter, result, Number<Current_Dim0 - 1>(),
                     Number<Current_Dim1>(), Number<Current_Dim2>());
  }

  template <class A, class B, class U, int Dim12, int Dim0, char i, char j,
            char k, int Current_Dim1, int Current_Dim2>
  void
  T3ch_equals_T3dg(A &iter, const Dg_Expr<B, U, Dim12, Dim0, i, j, k> &result,
                   const Number<1> &, const Number<Current_Dim1> &,
                   const Number<Current_Dim2> &)
  {
    iter(Current_Dim2 - 1, 0, Current_Dim1 - 1)
      = result(0, Current_Dim1 - 1, Current_Dim2 - 1);
    T3ch_equals_T3dg(iter, result, Number<Current_Dim1 - 1>(),
                     Number<Current_Dim1 - 1>(), Number<Current_Dim2>());
  }

  template <class A, class B, class U, int Dim12, int Dim0, char i, char j,
            char k, int Current_Dim2>
  void
  T3ch_equals_T3dg(A &iter, const Dg_Expr<B, U, Dim12, Dim0, i, j, k> &result,
                   const Number<1> &, const Number<1> &,
                   const Number<Current_Dim2> &)
  {
    iter(Current_Dim2 - 1, 0, 0) = result(0, 0, Current_Dim2 - 1);
    T3ch_equals_T3dg(iter, result, Number<Dim12>(), Number<Dim12>(),
                     Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class U, int Dim12, int Dim0, char i, char j,
            char k>
  void
  T3ch_equals_T3dg(A &iter, const Dg_Expr<B, U, Dim12, Dim0, i, j, k> &result,
                   const Number<1> &, const Number<1> &, const Number<1> &)
  {
    iter(0, 0, 0) = result(0, 0, 0);
  }

  template <class Tp, class T, int Tensor_Dim0, int Tensor_Dim12, int Dim12,
            int Dim0, char i, char j, char k>
  template <class B, class U>
  Dg_Expr<Christof<Tp, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k> &
  Dg_Expr<Christof<Tp, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k>::
  operator=(const Dg_Expr<B, U, Dim12, Dim0, i, j, k> &result)
  {
    T3ch_equals_T3dg(iter, result, Number<Dim12>(), Number<Dim12>(),
                     Number<Dim0>());
    return *this;
  }

  /* T3ch=T3ch_Expr(T3ch) */

  template <class Tp, class T, int Tensor_Dim0, int Tensor_Dim12, int Dim12,
            int Dim0, char i, char j, char k>
  Dg_Expr<Christof<Tp, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k> &
  Dg_Expr<Christof<Tp, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k>::
  operator=(const Dg_Expr<Christof<Tp, Tensor_Dim0, Tensor_Dim12>, T, Dim12,
                          Dim0, i, j, k> &result)
  {
    return operator=<Christof<Tp, Tensor_Dim0, Tensor_Dim12>, T>(result);
  }

  /* T3ch=U */

  template <class A, class U, int Dim12, int Dim0, int Current_Dim0,
            int Current_Dim1, int Current_Dim2>
  void T3ch_equals_generic(A &iter, const U &u, const Number<Current_Dim0> &,
                           const Number<Current_Dim1> &,
                           const Number<Current_Dim2> &, const Number<Dim12> &,
                           const Number<Dim0> &)
  {
    iter(Current_Dim2 - 1, Current_Dim0 - 1, Current_Dim1 - 1) = u;
    T3ch_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                        Number<Current_Dim1>(), Number<Current_Dim2>(),
                        Number<Dim12>(), Number<Dim0>());
  }

  template <class A, class U, int Dim12, int Dim0, int Current_Dim1,
            int Current_Dim2>
  void T3ch_equals_generic(A &iter, const U &u, const Number<1> &,
                           const Number<Current_Dim1> &,
                           const Number<Current_Dim2> &, const Number<Dim12> &,
                           const Number<Dim0> &)
  {
    iter(Current_Dim2 - 1, 0, Current_Dim1 - 1) = u;
    T3ch_equals_generic(iter, u, Number<Current_Dim1 - 1>(),
                        Number<Current_Dim1 - 1>(), Number<Current_Dim2>(),
                        Number<Dim12>(), Number<Dim0>());
  }

  template <class A, class U, int Dim12, int Dim0, int Current_Dim2>
  void T3ch_equals_generic(A &iter, const U &u, const Number<1> &,
                           const Number<1> &, const Number<Current_Dim2> &,
                           const Number<Dim12> &, const Number<Dim0> &)
  {
    iter(Current_Dim2 - 1, 0, 0) = u;
    T3ch_equals_generic(iter, u, Number<Dim12>(), Number<Dim12>(),
                        Number<Current_Dim2 - 1>(), Number<Dim12>(),
                        Number<Dim0>());
  }

  template <class A, class U, int Dim12, int Dim0>
  void T3ch_equals_generic(A &iter, const U &u, const Number<1> &,
                           const Number<1> &, const Number<1> &,
                           const Number<Dim12> &, const Number<Dim0> &)
  {
    iter(0, 0, 0) = u;
  }

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim12, int Dim12,
            int Dim0, char i, char j, char k>
  template <class U>
  Dg_Expr<Christof<A, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k> &
  Dg_Expr<Christof<A, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k>::
  operator=(const U &u)
  {
    T3ch_equals_generic(iter, u, Number<Dim12>(), Number<Dim12>(),
                        Number<Dim0>(), Number<Dim12>(), Number<Dim0>());
    return *this;
  }

  /* Assignment operators for the specializations to Ddg_number_rhs_0. */

  /* T4ddgrhs0=T3dg */

  template <class A, class B, class U, int Dim23, int Dim1, char i, char j,
            char k, int Current_Dim0, int Current_Dim1, int Current_Dim2, int N>
  void T4ddgrhs0_equals_T3dg(A &iter,
                             const Dg_Expr<B, U, Dim23, Dim1, i, j, k> &result,
                             const Number<Current_Dim0> &,
                             const Number<Current_Dim1> &,
                             const Number<Current_Dim2> &, const Number<N> &)
  {
    iter(N, Current_Dim2 - 1, Current_Dim0 - 1, Current_Dim1 - 1)
      = result(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1);
    T4ddgrhs0_equals_T3dg(iter, result, Number<Current_Dim0 - 1>(),
                          Number<Current_Dim1>(), Number<Current_Dim2>(),
                          Number<N>());
  }

  template <class A, class B, class U, int Dim23, int Dim1, char i, char j,
            char k, int Current_Dim1, int Current_Dim2, int N>
  void T4ddgrhs0_equals_T3dg(A &iter,
                             const Dg_Expr<B, U, Dim23, Dim1, i, j, k> &result,
                             const Number<1> &, const Number<Current_Dim1> &,
                             const Number<Current_Dim2> &, const Number<N> &)
  {
    iter(N, Current_Dim2 - 1, 0, Current_Dim1 - 1)
      = result(0, Current_Dim1 - 1, Current_Dim2 - 1);
    T4ddgrhs0_equals_T3dg(iter, result, Number<Current_Dim1 - 1>(),
                          Number<Current_Dim1 - 1>(), Number<Current_Dim2>(),
                          Number<N>());
  }

  template <class A, class B, class U, int Dim23, int Dim1, char i, char j,
            char k, int Current_Dim2, int N>
  void T4ddgrhs0_equals_T3dg(A &iter,
                             const Dg_Expr<B, U, Dim23, Dim1, i, j, k> &result,
                             const Number<1> &, const Number<1> &,
                             const Number<Current_Dim2> &, const Number<N> &)
  {
    iter(N, Current_Dim2 - 1, 0, 0) = result(0, 0, Current_Dim2 - 1);
    T4ddgrhs0_equals_T3dg(iter, result, Number<Dim23>(), Number<Dim23>(),
                          Number<Current_Dim2 - 1>(), Number<N>());
  }

  template <class A, class B, class U, int Dim23, int Dim1, char i, char j,
            char k, int N>
  void T4ddgrhs0_equals_T3dg(A &iter,
                             const Dg_Expr<B, U, Dim23, Dim1, i, j, k> &result,
                             const Number<1> &, const Number<1> &,
                             const Number<1> &, const Number<N> &)
  {
    iter(N, 0, 0, 0) = result(0, 0, 0);
  }

  template <class A, class T, int Dim23, int Dim1, char i, char j, char k,
            int N>
  template <class B, class U>
  Dg_Expr<Ddg_number_rhs_0<A, T, N>, T, Dim23, Dim1, i, j, k> &
  Dg_Expr<Ddg_number_rhs_0<A, T, N>, T, Dim23, Dim1, i, j, k>::
  operator=(const Dg_Expr<B, U, Dim23, Dim1, i, j, k> &result)
  {
    T4ddgrhs0_equals_T3dg(iter, result, Number<Dim23>(), Number<Dim23>(),
                          Number<Dim1>(), Number<N>());
    return *this;
  }

  /* T4ddgrhs0=T4ddgrhs0 */

  template <class A, class T, int Dim23, int Dim1, char i, char j, char k,
            int N>
  Dg_Expr<Ddg_number_rhs_0<A, T, N>, T, Dim23, Dim1, i, j, k> &
  Dg_Expr<Ddg_number_rhs_0<A, T, N>, T, Dim23, Dim1, i, j, k>::operator=(
    const Dg_Expr<Ddg_number_rhs_0<A, T, N>, T, Dim23, Dim1, i, j, k> &result)
  {
    return operator=<Ddg_number_rhs_0<A, T, N>, T>(result);
  }
}
