/* Various assignment operators.  I have to explicitly declare the
   second operator= because otherwise the compiler will generate its
   own and not use the template code. */

/* T2s=T2s */

#pragma once

namespace FTensor
{
  template <class A, class B, class U, int Dim, char i, char j,
            int Current_Dim0, int Current_Dim1>
  void
  T2s_equals_T2s(A &iter,
                 const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                 const Number<Current_Dim0> &, const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      = result(Current_Dim0 - 1, Current_Dim1 - 1);
    T2s_equals_T2s(iter, result, Number<Current_Dim0 - 1>(),
                   Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim, char i, char j,
            int Current_Dim1>
  void T2s_equals_T2s(A &iter,
                      const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                      const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) = result(0, Current_Dim1 - 1);
    T2s_equals_T2s(iter, result, Number<Current_Dim1 - 1>(),
                   Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim, char i, char j>
  void T2s_equals_T2s(A &iter,
                      const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                      const Number<1> &, const Number<1> &)
  {
    iter(0, 0) = result(0, 0);
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class B, class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result)
  {
    T2s_equals_T2s(iter, result, Number<Dim>(), Number<Dim>());
    return *this;
  }

  /* T2s=T2s_Expr(T2s) */

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator=(const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T,
                                         Dim, i, j> &result)
  {
    return operator=<Tensor2_symmetric<A, Tensor_Dim>, T>(result);
  }

  /* T2s+=T2s */

  template <class A, class B, class U, int Dim, char i, char j,
            int Current_Dim0, int Current_Dim1>
  void
  T2s_plus_equals_T2s(A &iter,
                      const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                      const Number<Current_Dim0> &,
                      const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      += result(Current_Dim0 - 1, Current_Dim1 - 1);
    T2s_plus_equals_T2s(iter, result, Number<Current_Dim0 - 1>(),
                        Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim, char i, char j,
            int Current_Dim1>
  void
  T2s_plus_equals_T2s(A &iter,
                      const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                      const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) += result(0, Current_Dim1 - 1);
    T2s_plus_equals_T2s(iter, result, Number<Current_Dim1 - 1>(),
                        Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim, char i, char j>
  void
  T2s_plus_equals_T2s(A &iter,
                      const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                      const Number<1> &, const Number<1> &)
  {
    iter(0, 0) += result(0, 0);
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class B, class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator+=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result)
  {
    T2s_plus_equals_T2s(iter, result, Number<Dim>(), Number<Dim>());
    return *this;
  }

  /* T2s-=T2s */

  template <class A, class B, class U, int Dim, char i, char j,
            int Current_Dim0, int Current_Dim1>
  void
  T2s_minus_equals_T2s(A &iter,
                       const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                       const Number<Current_Dim0> &,
                       const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      -= result(Current_Dim0 - 1, Current_Dim1 - 1);
    T2s_minus_equals_T2s(iter, result, Number<Current_Dim0 - 1>(),
                         Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim, char i, char j,
            int Current_Dim1>
  void
  T2s_minus_equals_T2s(A &iter,
                       const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                       const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) -= result(0, Current_Dim1 - 1);
    T2s_minus_equals_T2s(iter, result, Number<Current_Dim1 - 1>(),
                         Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim, char i, char j>
  void
  T2s_minus_equals_T2s(A &iter,
                       const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                       const Number<1> &, const Number<1> &)
  {
    iter(0, 0) -= result(0, 0);
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class B, class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator-=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result)
  {
    T2s_minus_equals_T2s(iter, result, Number<Dim>(), Number<Dim>());
    return *this;
  }

  /* T2s&=T2s */

  template <class A, class B, class U, int Dim, char i, char j,
            int Current_Dim0, int Current_Dim1>
  void
  T2s_and_equals_T2s(A &iter,
                     const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                     const Number<Current_Dim0> &,
                     const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1)
      *= result(Current_Dim0 - 1, Current_Dim1 - 1);
    T2s_and_equals_T2s(iter, result, Number<Current_Dim0 - 1>(),
                       Number<Current_Dim1>());
  }

  template <class A, class B, class U, int Dim, char i, char j,
            int Current_Dim1>
  void
  T2s_and_equals_T2s(A &iter,
                     const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                     const Number<1> &, const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) *= result(0, Current_Dim1 - 1);
    T2s_and_equals_T2s(iter, result, Number<Current_Dim1 - 1>(),
                       Number<Current_Dim1 - 1>());
  }

  template <class A, class B, class U, int Dim, char i, char j>
  void
  T2s_and_equals_T2s(A &iter,
                     const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result,
                     const Number<1> &, const Number<1> &)
  {
    iter(0, 0) *= result(0, 0);
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class B, class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator&=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result)
  {
    T2s_and_equals_T2s(iter, result, Number<Dim>(), Number<Dim>());
    return *this;
  }

  /* This is for when the indices are switched (i,j) -> (j,i). */

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class B, class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator=(const Tensor2_symmetric_Expr<B, U, Dim, j, i> &result)
  {
    T2s_equals_T2s(iter, result, Number<Dim>(), Number<Dim>());
    return *this;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class B, class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator+=(const Tensor2_symmetric_Expr<B, U, Dim, j, i> &result)
  {
    T2s_plus_equals_T2s(iter, result, Number<Dim>(), Number<Dim>());
    return *this;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class B, class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator-=(const Tensor2_symmetric_Expr<B, U, Dim, j, i> &result)
  {
    T2s_minus_equals_T2s(iter, result, Number<Dim>(), Number<Dim>());
    return *this;
  }

  /* Operations with just generics. */

  /* T2s=U */

  template <class A, class U, int Current_Dim0, int Current_Dim1>
  void T2s_equals_generic(A &iter, const U &u, const Number<Current_Dim0> &,
                          const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1) = u;
    T2s_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                       Number<Current_Dim1>());
  }

  template <class A, class U, int Current_Dim1>
  void T2s_equals_generic(A &iter, const U &u, const Number<1> &,
                          const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) = u;
    T2s_equals_generic(iter, u, Number<Current_Dim1 - 1>(),
                       Number<Current_Dim1 - 1>());
  }

  template <class A, class U>
  void
  T2s_equals_generic(A &iter, const U &u, const Number<1> &, const Number<1> &)
  {
    iter(0, 0) = u;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator=(const U &u)
  {
    T2s_equals_generic(iter, u, Number<Dim>(), Number<Dim>());
    return *this;
  }

  /* T2s+=U */

  template <class A, class U, int Current_Dim0, int Current_Dim1>
  void
  T2s_plus_equals_generic(A &iter, const U &u, const Number<Current_Dim0> &,
                          const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1) += u;
    T2s_plus_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                            Number<Current_Dim1>());
  }

  template <class A, class U, int Current_Dim1>
  void T2s_plus_equals_generic(A &iter, const U &u, const Number<1> &,
                               const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) += u;
    T2s_plus_equals_generic(iter, u, Number<Current_Dim1 - 1>(),
                            Number<Current_Dim1 - 1>());
  }

  template <class A, class U>
  void T2s_plus_equals_generic(A &iter, const U &u, const Number<1> &,
                               const Number<1> &)
  {
    iter(0, 0) += u;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator+=(const U &u)
  {
    T2s_plus_equals_generic(iter, u, Number<Dim>(), Number<Dim>());
    return *this;
  }

  /* T2s-=U */

  template <class A, class U, int Current_Dim0, int Current_Dim1>
  void
  T2s_minus_equals_generic(A &iter, const U &u, const Number<Current_Dim0> &,
                           const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1) -= u;
    T2s_minus_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                             Number<Current_Dim1>());
  }

  template <class A, class U, int Current_Dim1>
  void T2s_minus_equals_generic(A &iter, const U &u, const Number<1> &,
                                const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) -= u;
    T2s_minus_equals_generic(iter, u, Number<Current_Dim1 - 1>(),
                             Number<Current_Dim1 - 1>());
  }

  template <class A, class U>
  void T2s_minus_equals_generic(A &iter, const U &u, const Number<1> &,
                                const Number<1> &)
  {
    iter(0, 0) -= u;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator-=(const U &u)
  {
    T2s_minus_equals_generic(iter, u, Number<Dim>(), Number<Dim>());
    return *this;
  }

  /* T2s*=U */

  template <class A, class U, int Current_Dim0, int Current_Dim1>
  void
  T2s_times_equals_generic(A &iter, const U &u, const Number<Current_Dim0> &,
                           const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1) *= u;
    T2s_times_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                             Number<Current_Dim1>());
  }

  template <class A, class U, int Current_Dim1>
  void T2s_times_equals_generic(A &iter, const U &u, const Number<1> &,
                                const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) *= u;
    T2s_times_equals_generic(iter, u, Number<Current_Dim1 - 1>(),
                             Number<Current_Dim1 - 1>());
  }

  template <class A, class U>
  void T2s_times_equals_generic(A &iter, const U &u, const Number<1> &,
                                const Number<1> &)
  {
    iter(0, 0) *= u;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator*=(const U &u)
  {
    T2s_times_equals_generic(iter, u, Number<Dim>(), Number<Dim>());
    return *this;
  }

  /* T2s/=U */

  template <class A, class U, int Current_Dim0, int Current_Dim1>
  void
  T2s_divide_equals_generic(A &iter, const U &u, const Number<Current_Dim0> &,
                            const Number<Current_Dim1> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1) /= u;
    T2s_divide_equals_generic(iter, u, Number<Current_Dim0 - 1>(),
                              Number<Current_Dim1>());
  }

  template <class A, class U, int Current_Dim1>
  void T2s_divide_equals_generic(A &iter, const U &u, const Number<1> &,
                                 const Number<Current_Dim1> &)
  {
    iter(0, Current_Dim1 - 1) /= u;
    T2s_divide_equals_generic(iter, u, Number<Current_Dim1 - 1>(),
                              Number<Current_Dim1 - 1>());
  }

  template <class A, class U>
  void T2s_divide_equals_generic(A &iter, const U &u, const Number<1> &,
                                 const Number<1> &)
  {
    iter(0, 0) /= u;
  }

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  template <class U>
  const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>::
  operator/=(const U &u)
  {
    T2s_divide_equals_generic(iter, u, Number<Dim>(), Number<Dim>());
    return *this;
  }

  /* Various assignment operators for Dg_number_rhs_2.  I have
     to explicitly declare the second operator= because otherwise the
     compiler will generate its own and not use the template code. */

  template <class A, class B, class U, int Dim01, char i, char j, int N,
            int Current_Dim0, int Current_Dim1>
  void
  T3dgrhs2_equals_T2s(A &iter,
                      const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &result,
                      const Number<Current_Dim0> &,
                      const Number<Current_Dim1> &, const Number<N> &)
  {
    iter(Current_Dim0 - 1, Current_Dim1 - 1, N)
      = result(Current_Dim0 - 1, Current_Dim1 - 1);
    T3dgrhs2_equals_T2s(iter, result, Number<Current_Dim0 - 1>(),
                        Number<Current_Dim1>(), Number<N>());
  }

  template <class A, class B, class U, int Dim01, char i, char j, int N,
            int Current_Dim1>
  void
  T3dgrhs2_equals_T2s(A &iter,
                      const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &result,
                      const Number<1> &, const Number<Current_Dim1> &,
                      const Number<N> &)
  {
    iter(0, Current_Dim1 - 1, N) = result(0, Current_Dim1 - 1);
    T3dgrhs2_equals_T2s(iter, result, Number<Current_Dim1 - 1>(),
                        Number<Current_Dim1 - 1>(), Number<N>());
  }

  template <class A, class B, class U, int Dim01, char i, char j, int N>
  void
  T3dgrhs2_equals_T2s(A &iter,
                      const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &result,
                      const Number<1> &, const Number<1> &, const Number<N> &)
  {
    iter(0, 0, N) = result(0, 0);
  }

  template <class A, class T, int Dim, char i, char j, int N>
  template <class B, class U>
  const Tensor2_symmetric_Expr<Dg_number_rhs_2<A, T, N>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Dg_number_rhs_2<A, T, N>, T, Dim, i, j>::
  operator=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result)
  {
    T3dgrhs2_equals_T2s(iter, result, Number<Dim>(), Number<Dim>(),
                        Number<N>());
    return *this;
  }

  template <class A, class T, int Dim, char i, char j, int N>
  const Tensor2_symmetric_Expr<Dg_number_rhs_2<A, T, N>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Dg_number_rhs_2<A, T, N>, T, Dim, i, j>::
  operator=(const Tensor2_symmetric_Expr<Dg_number_rhs_2<A, T, N>, T, Dim, i,
                                         j> &result)
  {
    return operator=<Dg_number_rhs_2<A, T, N>, T>(result);
  }

  /* Various assignment operators for Ddg_number_rhs_01.  I have
     to explicitly declare the second operator= because otherwise the
     compiler will generate its own and not use the template code. */

  template <class A, class B, class U, int Dim01, char i, char j, int N0,
            int N1, int Current_Dim0, int Current_Dim1>
  void T4ddgrhs01_equals_T2s(
    A &iter, const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &result,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
    const Number<N0> &, const Number<N1> &)
  {
    iter(N0, N1, Current_Dim0 - 1, Current_Dim1 - 1)
      = result(Current_Dim0 - 1, Current_Dim1 - 1);
    T4ddgrhs01_equals_T2s(iter, result, Number<Current_Dim0 - 1>(),
                          Number<Current_Dim1>(), Number<N0>(), Number<N1>());
  }

  template <class A, class B, class U, int Dim01, char i, char j, int N0,
            int N1, int Current_Dim1>
  void T4ddgrhs01_equals_T2s(
    A &iter, const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &result,
    const Number<1> &, const Number<Current_Dim1> &, const Number<N0> &,
    const Number<N1> &)
  {
    iter(N0, N1, 0, Current_Dim1 - 1) = result(0, Current_Dim1 - 1);
    T4ddgrhs01_equals_T2s(iter, result, Number<Current_Dim1 - 1>(),
                          Number<Current_Dim1 - 1>(), Number<N0>(),
                          Number<N1>());
  }

  template <class A, class B, class U, int Dim01, char i, char j, int N0,
            int N1>
  void T4ddgrhs01_equals_T2s(
    A &iter, const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &result,
    const Number<1> &, const Number<1> &, const Number<N0> &,
    const Number<N1> &)
  {
    iter(N0, N1, 0, 0) = result(0, 0);
  }
  /*
  template <class A, class T, int Dim, char i, char j, int N0, int N1>
  template <class B, class U>
  const Tensor2_symmetric_Expr<Ddg_number_rhs_01<A, T, N0, N1>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Ddg_number_rhs_01<A, T, N0, N1>, T, Dim, i, j>::
  operator=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result)
  {
    T4ddgrhs01_equals_T2s(iter, result, Number<Dim>(), Number<Dim>(),
                          Number<N0>(), Number<N1>());
    return *this;
  }

  template <class A, class T, int Dim, char i, char j, int N0, int N1>
  const Tensor2_symmetric_Expr<Ddg_number_rhs_01<A, T, N0, N1>, T, Dim, i, j> &
  Tensor2_symmetric_Expr<Ddg_number_rhs_01<A, T, N0, N1>, T, Dim, i, j>::
  operator=(const Tensor2_symmetric_Expr<Ddg_number_rhs_01<A, T, N0, N1>, T,
                                         Dim, i, j> &result)
  {
    return operator=<Ddg_number_rhs_01<A, T, N0, N1>, T>(result);
    }*/
}
