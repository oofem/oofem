/* Various assignment operators for generic Tensor1_Expr's as well as
   specializations for Tensor2_number_rhs's.  I have to explicitly
   declare the second operator= because otherwise the compiler will
   generate its own and not use the template code. */

/* =T1_Expr */

#pragma once

namespace FTensor
{
  template <class A, class B, class U, int Dim, char i, int Current_Dim>
  void T1_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim, i> result,
                    const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1) = result(Current_Dim - 1);
    T1_equals_T1(iter, result, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim, char i>
  void T1_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim, i> result,
                    const Number<1> &)
  {
    iter(0) = result(0);
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i>
  template <class B, class U>
  const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
  Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i>::
  operator=(const Tensor1_Expr<B, U, Dim, i> &result)
  {
    T1_equals_T1(iter, result, Number<Dim>());
    return *this;
  }

  /* =T1_Expr(T1) */

  template <class A, class T, int Tensor_Dim, int Dim, char i>
  const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
  Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i>::
  operator=(const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &result)
  {
    return operator=<Tensor1<A, Tensor_Dim>, T>(result);
  }

  /* +=T1 */

  template <class A, class B, class U, int Dim, char i, int Current_Dim>
  void T1_plus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim, i> result,
                         const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1) += result(Current_Dim - 1);
    T1_plus_equals_T1(iter, result, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim, char i>
  void T1_plus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim, i> result,
                         const Number<1> &)
  {
    iter(0) += result(0);
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i>
  template <class B, class U>
  const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
  Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i>::
  operator+=(const Tensor1_Expr<B, U, Dim, i> &result)
  {
    T1_plus_equals_T1(iter, result, Number<Dim>());
    return *this;
  }

  /* -=T1 */

  template <class A, class B, class U, int Dim, char i, int Current_Dim>
  void T1_minus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim, i> result,
                          const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1) -= result(Current_Dim - 1);
    T1_minus_equals_T1(iter, result, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim, char i>
  void T1_minus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim, i> result,
                          const Number<1> &)
  {
    iter(0) -= result(0);
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i>
  template <class B, class U>
  const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
  Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i>::
  operator-=(const Tensor1_Expr<B, U, Dim, i> &result)
  {
    T1_minus_equals_T1(iter, result, Number<Dim>());
    return *this;
  }

  /* General template assignment operators intended mostly for
     doubles (type U), but could be applied to anything you want, like
     complex, etc.  All that is required is that T=U works (or T+=U,
     etc.)  */

  /* =U */

  template <class A, class U, int Current_Dim>
  void T1_equals_generic(A &iter, const U &u, const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1) = u;
    T1_equals_generic(iter, u, Number<Current_Dim - 1>());
  }

  template <class A, class U>
  void T1_equals_generic(A &iter, const U &u, const Number<1> &)
  {
    iter(0) = u;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i>
  template <class U>
  const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
  Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i>::operator=(const U &u)
  {
    T1_equals_generic(iter, u, Number<Dim>());
    return *this;
  }

  /* +=U */

  template <class A, class U, int Current_Dim>
  void T1_plus_equals_generic(A &iter, const U &u, const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1) += u;
    T1_plus_equals_generic(iter, u, Number<Current_Dim - 1>());
  }

  template <class A, class U>
  void T1_plus_equals_generic(A &iter, const U &u, const Number<1> &)
  {
    iter(0) += u;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i>
  template <class U>
  const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
  Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i>::operator+=(const U &u)
  {
    T1_plus_equals_generic(iter, u, Number<Dim>());
    return *this;
  }

  /* -=U */

  template <class A, class U, int Current_Dim>
  void
  T1_minus_equals_generic(A &iter, const U &u, const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1) -= u;
    T1_minus_equals_generic(iter, u, Number<Current_Dim - 1>());
  }

  template <class A, class U>
  void T1_minus_equals_generic(A &iter, const U &u, const Number<1> &)
  {
    iter(0) -= u;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i>
  template <class U>
  const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
  Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i>::operator-=(const U &u)
  {
    T1_minus_equals_generic(iter, u, Number<Dim>());
    return *this;
  }

  /* *=U */

  template <class A, class U, int Current_Dim>
  void
  T1_times_equals_generic(A &iter, const U &u, const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1) *= u;
    T1_times_equals_generic(iter, u, Number<Current_Dim - 1>());
  }

  template <class A, class U>
  void T1_times_equals_generic(A &iter, const U &u, const Number<1> &)
  {
    iter(0) *= u;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i>
  template <class U>
  const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
  Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i>::operator*=(const U &u)
  {
    T1_times_equals_generic(iter, u, Number<Dim>());
    return *this;
  }

  /* /=U */

  template <class A, class U, int Current_Dim>
  void
  T1_divide_equals_generic(A &iter, const U &u, const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1) /= u;
    T1_divide_equals_generic(iter, u, Number<Current_Dim - 1>());
  }

  template <class A, class U>
  void T1_divide_equals_generic(A &iter, const U &u, const Number<1> &)
  {
    iter(0) /= u;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i>
  template <class U>
  const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
  Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i>::operator/=(const U &u)
  {
    T1_divide_equals_generic(iter, u, Number<Dim>());
    return *this;
  }

  /* Specializations for Tensor2_number_rhs_0's */

  /* =T1 */

  template <class A, class B, class U, int Dim1, char i, int N, int Current_Dim>
  void T2rhs0_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                        const Number<N> &N1, const Number<Current_Dim> &)
  {
    iter(N, Current_Dim - 1) = result(Current_Dim - 1);
    T2rhs0_equals_T1(iter, result, N1, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim1, char i, int N>
  void T2rhs0_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                        const Number<N> &, const Number<1> &)
  {
    iter(N, 0) = result(0);
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class B, class U>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::
  operator=(const Tensor1_Expr<B, U, Dim1, i> &result)
  {
    T2rhs0_equals_T1(iter, result, Number<N>(), Number<Dim1>());
    return *this;
  }

  template <class A, class T, int Dim1, char i, int N>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::operator=(
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &result)
  {
    return operator=<Tensor2_number_rhs_0<A, T, N>, T>(result);
  }

  /* +=T1 */

  template <class A, class B, class U, int Dim1, char i, int N, int Current_Dim>
  void
  T2rhs0_plus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                        const Number<N> &N1, const Number<Current_Dim> &)
  {
    iter(N, Current_Dim - 1) += result(Current_Dim - 1);
    T2rhs0_plus_equals_T1(iter, result, N1, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim1, char i, int N>
  void
  T2rhs0_plus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                        const Number<N> &, const Number<1> &)
  {
    iter(N, 0) += result(0);
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class B, class U>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::
  operator+=(const Tensor1_Expr<B, U, Dim1, i> &result)
  {
    T2rhs0_plus_equals_T1(iter, result, Number<N>(), Number<Dim1>());
    return *this;
  }

  template <class A, class T, int Dim1, char i, int N>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::operator+=(
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &result)
  {
    return operator+=<Tensor2_number_rhs_0<A, T, N>, T>(result);
  }

  /* -=T1 */

  template <class A, class B, class U, int Dim1, char i, int N, int Current_Dim>
  void
  T2rhs0_minus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                         const Number<N> &N1, const Number<Current_Dim> &)
  {
    iter(N, Current_Dim - 1) -= result(Current_Dim - 1);
    T2rhs0_minus_equals_T1(iter, result, N1, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim1, char i, int N>
  void
  T2rhs0_minus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                         const Number<N> &, const Number<1> &)
  {
    iter(N, 0) -= result(0);
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class B, class U>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::
  operator-=(const Tensor1_Expr<B, U, Dim1, i> &result)
  {
    T2rhs0_minus_equals_T1(iter, result, Number<N>(), Number<Dim1>());
    return *this;
  }

  template <class A, class T, int Dim1, char i, int N>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::operator-=(
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &result)
  {
    return operator-=<Tensor2_number_rhs_0<A, T, N>, T>(result);
  }

  /* =U */

  template <class A, class U, int N, int Current_Dim>
  void T2rhs0_equals_generic(A &iter, const U &u, const Number<N> &N1,
                             const Number<Current_Dim> &)
  {
    iter(N, Current_Dim - 1) = u;
    T2rhs0_equals_generic(iter, u, N1, Number<Current_Dim - 1>());
  }

  template <class A, class U, int N>
  void T2rhs0_equals_generic(A &iter, const U &u, const Number<N> &,
                             const Number<1> &)
  {
    iter(N, 0) = u;
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class U>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::
  operator=(const U &u)
  {
    T2rhs0_equals_generic(iter, u, Number<N>(), Number<Dim1>());
    return *this;
  }

  /* +=U */

  template <class A, class U, int N, int Current_Dim>
  void T2rhs0_plus_equals_generic(A &iter, const U &u, const Number<N> &N1,
                                  const Number<Current_Dim> &)
  {
    iter(N, Current_Dim - 1) += u;
    T2rhs0_plus_equals_generic(iter, u, N1, Number<Current_Dim - 1>());
  }

  template <class A, class U, int N>
  void T2rhs0_plus_equals_generic(A &iter, const U &u, const Number<N> &,
                                  const Number<1> &)
  {
    iter(N, 0) += u;
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class U>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::
  operator+=(const U &u)
  {
    T2rhs0_plus_equals_generic(iter, u, Number<N>(), Number<Dim1>());
    return *this;
  }

  /* -=U */

  template <class A, class U, int N, int Current_Dim>
  void T2rhs0_minus_equals_generic(A &iter, const U &u, const Number<N> &N1,
                                   const Number<Current_Dim> &)
  {
    iter(N, Current_Dim - 1) -= u;
    T2rhs0_minus_equals_generic(iter, u, N1, Number<Current_Dim - 1>());
  }

  template <class A, class U, int N>
  void T2rhs0_minus_equals_generic(A &iter, const U &u, const Number<N> &,
                                   const Number<1> &)
  {
    iter(N, 0) -= u;
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class U>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::
  operator-=(const U &u)
  {
    T2rhs0_minus_equals_generic(iter, u, Number<N>(), Number<Dim1>());
    return *this;
  }

  /* *=U */

  template <class A, class U, int N, int Current_Dim>
  void T2rhs0_times_equals_generic(A &iter, const U &u, const Number<N> &N1,
                                   const Number<Current_Dim> &)
  {
    iter(N, Current_Dim - 1) *= u;
    T2rhs0_times_equals_generic(iter, u, N1, Number<Current_Dim - 1>());
  }

  template <class A, class U, int N>
  void T2rhs0_times_equals_generic(A &iter, const U &u, const Number<N> &,
                                   const Number<1> &)
  {
    iter(N, 0) *= u;
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class U>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::
  operator*=(const U &u)
  {
    T2rhs0_times_equals_generic(iter, u, Number<N>(), Number<Dim1>());
    return *this;
  }

  /* /=U */

  template <class A, class U, int N, int Current_Dim>
  void T2rhs0_divide_equals_generic(A &iter, const U &u, const Number<N> &N1,
                                    const Number<Current_Dim> &)
  {
    iter(N, Current_Dim - 1) /= u;
    T2rhs0_divide_equals_generic(iter, u, N1, Number<Current_Dim - 1>());
  }

  template <class A, class U, int N>
  void T2rhs0_divide_equals_generic(A &iter, const U &u, const Number<N> &,
                                    const Number<1> &)
  {
    iter(N, 0) /= u;
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class U>
  const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>::
  operator/=(const U &u)
  {
    T2rhs0_divide_equals_generic(iter, u, Number<N>(), Number<Dim1>());
    return *this;
  }

  /* Specializations for Tensor2_number_rhs_1's */

  /* =T1 */

  template <class A, class B, class U, int Dim1, char i, int N, int Current_Dim>
  void T2rhs1_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                        const Number<N> &N1, const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1, N) = result(Current_Dim - 1);
    T2rhs1_equals_T1(iter, result, N1, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim1, char i, int N>
  void T2rhs1_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                        const Number<N> &, const Number<1> &)
  {
    iter(0, N) = result(0);
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class B, class U>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::
  operator=(const Tensor1_Expr<B, U, Dim1, i> &result)
  {
    T2rhs1_equals_T1(iter, result, Number<N>(), Number<Dim1>());
    return *this;
  }

  template <class A, class T, int Dim1, char i, int N>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::operator=(
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &result)
  {
    return operator=<Tensor2_number_rhs_1<A, T, N>, T>(result);
  }

  /* +=T1 */

  template <class A, class B, class U, int Dim1, char i, int N, int Current_Dim>
  void
  T2rhs1_plus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                        const Number<N> &N1, const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1, N) += result(Current_Dim - 1);
    T2rhs1_plus_equals_T1(iter, result, N1, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim1, char i, int N>
  void
  T2rhs1_plus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                        const Number<N> &, const Number<1> &)
  {
    iter(0, N) += result(0);
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class B, class U>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::
  operator+=(const Tensor1_Expr<B, U, Dim1, i> &result)
  {
    T2rhs1_plus_equals_T1(iter, result, Number<N>(), Number<Dim1>());
    return *this;
  }

  template <class A, class T, int Dim1, char i, int N>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::operator+=(
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &result)
  {
    return operator+=<Tensor2_number_rhs_1<A, T, N>, T>(result);
  }

  /* -=T1 */

  template <class A, class B, class U, int Dim1, char i, int N, int Current_Dim>
  void
  T2rhs1_minus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                         const Number<N> &N1, const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1, N) -= result(Current_Dim - 1);
    T2rhs1_minus_equals_T1(iter, result, N1, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim1, char i, int N>
  void
  T2rhs1_minus_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim1, i> &result,
                         const Number<N> &, const Number<1> &)
  {
    iter(0, N) -= result(0);
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class B, class U>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::
  operator-=(const Tensor1_Expr<B, U, Dim1, i> &result)
  {
    T2rhs1_minus_equals_T1(iter, result, Number<N>(), Number<Dim1>());
    return *this;
  }

  template <class A, class T, int Dim1, char i, int N>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::operator-=(
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &result)
  {
    return operator-=<Tensor2_number_rhs_1<A, T, N>, T>(result);
  }

  /* =U */

  template <class A, class U, int N, int Current_Dim>
  void T2rhs1_equals_generic(A &iter, const U &u, const Number<N> &N1,
                             const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1, N) = u;
    T2rhs1_equals_generic(iter, u, N1, Number<Current_Dim - 1>());
  }

  template <class A, class U, int N>
  void T2rhs1_equals_generic(A &iter, const U &u, const Number<N> &,
                             const Number<1> &)
  {
    iter(0, N) = u;
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class U>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::
  operator=(const U &u)
  {
    T2rhs1_equals_generic(iter, u, Number<N>(), Number<Dim1>());
    return *this;
  }

  /* +=U */

  template <class A, class U, int N, int Current_Dim>
  void T2rhs1_plus_equals_generic(A &iter, const U &u, const Number<N> &N1,
                                  const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1, N) += u;
    T2rhs1_plus_equals_generic(iter, u, N1, Number<Current_Dim - 1>());
  }

  template <class A, class U, int N>
  void T2rhs1_plus_equals_generic(A &iter, const U &u, const Number<N> &,
                                  const Number<1> &)
  {
    iter(0, N) += u;
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class U>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::
  operator+=(const U &u)
  {
    T2rhs1_plus_equals_generic(iter, u, Number<N>(), Number<Dim1>());
    return *this;
  }

  /* -=U */

  template <class A, class U, int N, int Current_Dim>
  void T2rhs1_minus_equals_generic(A &iter, const U &u, const Number<N> &N1,
                                   const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1, N) -= u;
    T2rhs1_minus_equals_generic(iter, u, N1, Number<Current_Dim - 1>());
  }

  template <class A, class U, int N>
  void T2rhs1_minus_equals_generic(A &iter, const U &u, const Number<N> &,
                                   const Number<1> &)
  {
    iter(0, N) -= u;
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class U>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::
  operator-=(const U &u)
  {
    T2rhs1_minus_equals_generic(iter, u, Number<N>(), Number<Dim1>());
    return *this;
  }

  /* *=U */

  template <class A, class U, int N, int Current_Dim>
  void T2rhs1_times_equals_generic(A &iter, const U &u, const Number<N> &N1,
                                   const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1, N) *= u;
    T2rhs1_times_equals_generic(iter, u, N1, Number<Current_Dim - 1>());
  }

  template <class A, class U, int N>
  void T2rhs1_times_equals_generic(A &iter, const U &u, const Number<N> &,
                                   const Number<1> &)
  {
    iter(0, N) *= u;
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class U>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::
  operator*=(const U &u)
  {
    T2rhs1_times_equals_generic(iter, u, Number<N>(), Number<Dim1>());
    return *this;
  }

  /* /=U */

  template <class A, class U, int N, int Current_Dim>
  void T2rhs1_divide_equals_generic(A &iter, const U &u, const Number<N> &N1,
                                    const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1, N) /= u;
    T2rhs1_divide_equals_generic(iter, u, N1, Number<Current_Dim - 1>());
  }

  template <class A, class U, int N>
  void T2rhs1_divide_equals_generic(A &iter, const U &u, const Number<N> &,
                                    const Number<1> &)
  {
    iter(0, N) /= u;
  }

  template <class A, class T, int Dim1, char i, int N>
  template <class U>
  const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
  Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>::
  operator/=(const U &u)
  {
    T2rhs1_divide_equals_generic(iter, u, Number<N>(), Number<Dim1>());
    return *this;
  }

  /* Specializations for Dg_number_rhs_12's */

  template <class A, class B, class U, int Dim, char i, int N1, int N2,
            int Current_Dim>
  void T3dgrhs12_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim, i> &result,
                           const Number<N1> &NN1, const Number<N2> &NN2,
                           const Number<Current_Dim> &)
  {
    iter(Current_Dim - 1, N1, N2) = result(Current_Dim - 1);
    T3dgrhs12_equals_T1(iter, result, NN1, NN2, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim, char i, int N1, int N2>
  void T3dgrhs12_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim, i> &result,
                           const Number<N1> &, const Number<N2> &,
                           const Number<1> &)
  {
    iter(0, N1, N2) = result(0);
  }

  template <class A, class T, int Dim, char i, int N1, int N2>
  template <class B, class U>
  const Tensor1_Expr<Dg_number_rhs_12<A, T, N1, N2>, T, Dim, i> &
  Tensor1_Expr<Dg_number_rhs_12<A, T, N1, N2>, T, Dim, i>::
  operator=(const Tensor1_Expr<B, U, Dim, i> &result)
  {
    T3dgrhs12_equals_T1(iter, result, Number<N1>(), Number<N2>(),
                        Number<Dim>());
    return *this;
  }

  template <class A, class T, int Dim, char i, int N1, int N2>
  const Tensor1_Expr<Dg_number_rhs_12<A, T, N1, N2>, T, Dim, i> &
  Tensor1_Expr<Dg_number_rhs_12<A, T, N1, N2>, T, Dim, i>::operator=(
    const Tensor1_Expr<Dg_number_rhs_12<A, T, N1, N2>, T, Dim, i> &result)
  {
    return operator=<Dg_number_rhs_12<A, T, N1, N2>, T>(result);
  }

  /* Specializations for Dg_number_rhs_01's */

  template <class A, class B, class U, int Dim, char i, int N1, int N2,
            int Current_Dim>
  void T3dgrhs01_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim, i> &result,
                           const Number<N1> &NN1, const Number<N2> &NN2,
                           const Number<Current_Dim> &)
  {
    iter(N1, N2, Current_Dim - 1) = result(Current_Dim - 1);
    T3dgrhs01_equals_T1(iter, result, NN1, NN2, Number<Current_Dim - 1>());
  }

  template <class A, class B, class U, int Dim, char i, int N1, int N2>
  void T3dgrhs01_equals_T1(A &iter, const Tensor1_Expr<B, U, Dim, i> &result,
                           const Number<N1> &, const Number<N2> &,
                           const Number<1> &)
  {
    iter(N1, N2, 0) = result(0);
  }

  template <class A, class T, int Dim, char i, int N1, int N2>
  template <class B, class U>
  const Tensor1_Expr<Dg_number_rhs_01<A, T, N1, N2>, T, Dim, i> &
  Tensor1_Expr<Dg_number_rhs_01<A, T, N1, N2>, T, Dim, i>::
  operator=(const Tensor1_Expr<B, U, Dim, i> &result)
  {
    T3dgrhs01_equals_T1(iter, result, Number<N1>(), Number<N2>(),
                        Number<Dim>());
    return *this;
  }

  template <class A, class T, int Dim, char i, int N1, int N2>
  const Tensor1_Expr<Dg_number_rhs_01<A, T, N1, N2>, T, Dim, i> &
  Tensor1_Expr<Dg_number_rhs_01<A, T, N1, N2>, T, Dim, i>::operator=(
    const Tensor1_Expr<Dg_number_rhs_01<A, T, N1, N2>, T, Dim, i> &result)
  {
    return operator=<Dg_number_rhs_01<A, T, N1, N2>, T>(result);
  }
}
