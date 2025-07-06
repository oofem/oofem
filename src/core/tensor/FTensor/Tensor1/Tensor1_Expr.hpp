/* Declares a wrapper class for rank 1 Tensor expressions. Note that
   Tensor1_Expr_equals is included at the end, since it needs the
   definition of the class in order to compile. */

#pragma once

#include "Tensor1_and_Tensor1.hpp"
#include "Tensor1_carat_Tensor1.hpp"
#include "Tensor1_divide_generic.hpp"
#include "Tensor1_minus_Tensor1.hpp"
#include "Tensor1_minus_generic.hpp"
#include "Tensor1_or_Tensor1.hpp"
#include "Tensor1_plus_Tensor1.hpp"
#include "Tensor1_plus_generic.hpp"
#include "Tensor1_times_Tensor1.hpp"
#include "Tensor1_times_generic.hpp"
#include "dTensor1.hpp"
#include "d_one_sided_Tensor1.hpp"
#include "ddTensor1.hpp"
#include "diffusion_Tensor1.hpp"
#include "generic_minus_Tensor1.hpp"
#include "interpolate_Tensor1.hpp"
#include "minus_Tensor1.hpp"

namespace FTensor
{
  template <class A, class T, int Dim, char i> class Tensor1_Expr
  {
    A iter;

  public:
    Tensor1_Expr(const A &a) : iter(a) {}
    T operator()(const int N) const { return iter(N); }
  };

  template <class A, class T, int Tensor_Dim, int Dim, char i>
  class Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i>
  {
    Tensor1<A, Tensor_Dim> &iter;

  public:
    Tensor1_Expr(Tensor1<A, Tensor_Dim> &a) : iter(a) {}

    T &operator()(const int N) { return iter(N); }
    T operator()(const int N) const { return iter(N); }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
    operator=(const Tensor1_Expr<B, U, Dim, i> &result);

    const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
    operator=(const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &result);

    template <class B, class U>
    const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
    operator+=(const Tensor1_Expr<B, U, Dim, i> &result);

    template <class B, class U>
    const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
    operator-=(const Tensor1_Expr<B, U, Dim, i> &result);

    /* General template assignment operators intended mostly for
       doubles (type T), but could be applied to anything you want, like
       complex, etc.  All that is required is that T=B works (or T+=B,
       etc.)  */

    template <class B>
    const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
    operator=(const B &d);
    template <class B>
    const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
    operator+=(const B &d);
    template <class B>
    const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
    operator-=(const B &d);
    template <class B>
    const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
    operator*=(const B &d);
    template <class B>
    const Tensor1_Expr<Tensor1<A, Tensor_Dim>, T, Dim, i> &
    operator/=(const B &d);
  };

  /* Specialized for Tensor2_number_rhs_0 (Tensor2{_symmetric} with the
     first index explicitly given). */

  template <class A, class T, int Dim1, char i, int N>
  class Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i>
  {
    A &iter;

  public:
    Tensor1_Expr(A &a) : iter(a) {}
    T &operator()(const int N1) { return iter(N, N1); }
    T operator()(const int N1) const { return iter(N, N1); }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
    operator=(const Tensor1_Expr<B, U, Dim1, i> &result);

    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &operator=(
      const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &result);

    template <class B, class U>
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
    operator+=(const Tensor1_Expr<B, U, Dim1, i> &result);

    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &operator+=(
      const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &result);

    template <class B, class U>
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
    operator-=(const Tensor1_Expr<B, U, Dim1, i> &result);

    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &operator-=(
      const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &result);

    template <class B>
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
    operator=(const B &result);

    template <class B>
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
    operator+=(const B &result);

    template <class B>
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
    operator-=(const B &result);

    template <class B>
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
    operator*=(const B &result);

    template <class B>
    const Tensor1_Expr<Tensor2_number_rhs_0<A, T, N>, T, Dim1, i> &
    operator/=(const B &result);
  };

  /* Specialized for Tensor2_number_rhs_1 (Tensor2{_symmetric} with the
     second index explicitly given). */

  template <class A, class T, int Dim1, char i, int N>
  class Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i>
  {
    A &iter;

  public:
    Tensor1_Expr(A &a) : iter(a) {}
    T &operator()(const int N1) { return iter(N1, N); }
    T operator()(const int N1) const { return iter(N1, N); }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
    operator=(const Tensor1_Expr<B, U, Dim1, i> &result);

    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &operator=(
      const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &result);

    template <class B, class U>
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
    operator+=(const Tensor1_Expr<B, U, Dim1, i> &result);

    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &operator+=(
      const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &result);

    template <class B, class U>
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
    operator-=(const Tensor1_Expr<B, U, Dim1, i> &result);

    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &operator-=(
      const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &result);

    template <class B>
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
    operator=(const B &result);

    template <class B>
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
    operator+=(const B &result);

    template <class B>
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
    operator-=(const B &result);

    template <class B>
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
    operator*=(const B &result);

    template <class B>
    const Tensor1_Expr<Tensor2_number_rhs_1<A, T, N>, T, Dim1, i> &
    operator/=(const B &result);
  };

  /* Specialized for Dg_number_rhs_12 (A Dg with
     explicit numbers in the (first or second) and third slots). */

  template <class A, class T, int Dim, char i, int N1, int N2>
  class Tensor1_Expr<Dg_number_rhs_12<A, T, N1, N2>, T, Dim, i>
  {
    A &iter;

  public:
    Tensor1_Expr(A &a) : iter(a) {}
    T &operator()(const int N) { return iter(N, N1, N2); }
    T operator()(const int N) const { return iter(N, N1, N2); }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    const Tensor1_Expr<Dg_number_rhs_12<A, T, N1, N2>, T, Dim, i> &
    operator=(const Tensor1_Expr<B, U, Dim, i> &result);

    const Tensor1_Expr<Dg_number_rhs_12<A, T, N1, N2>, T, Dim, i> &operator=(
      const Tensor1_Expr<Dg_number_rhs_12<A, T, N1, N2>, T, Dim, i> &result);
  };

  /* Specialized for Dg_number_rhs_01 (A Dg with
     explicit numbers in the first and second slots). */

  template <class A, class T, int Dim, char i, int N1, int N2>
  class Tensor1_Expr<Dg_number_rhs_01<A, T, N1, N2>, T, Dim, i>
  {
    A &iter;

  public:
    Tensor1_Expr(A &a) : iter(a) {}
    T &operator()(const int N) { return iter(N1, N2, N); }
    T operator()(const int N) const { return iter(N1, N2, N); }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    const Tensor1_Expr<Dg_number_rhs_01<A, T, N1, N2>, T, Dim, i> &
    operator=(const Tensor1_Expr<B, U, Dim, i> &result);

    const Tensor1_Expr<Dg_number_rhs_01<A, T, N1, N2>, T, Dim, i> &operator=(
      const Tensor1_Expr<Dg_number_rhs_01<A, T, N1, N2>, T, Dim, i> &result);
  };
}

#include "Tensor1_Expr_equals.hpp"
