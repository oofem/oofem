/* Declares a wrapper class for symmetric rank 2 Tensor expressions.
   It is specialized for Tensor3_number_rhs_2.  Note that
   Tensor2_symmetric_Expr_equals.hpp is included at the end, because it
   needs the definition of Tensor2_symmetric_Expr. */

#pragma once

#include "Tensor2_symmetric_and_Tensor2_symmetric.hpp"
#include "Tensor2_symmetric_carat_Tensor2.hpp"
#include "Tensor2_symmetric_divide_generic.hpp"
#include "Tensor2_symmetric_minus_Tensor2.hpp"
#include "Tensor2_symmetric_minus_Tensor2_symmetric.hpp"
#include "Tensor2_symmetric_minus_generic.hpp"
#include "Tensor2_symmetric_mod_Tensor2_symmetric.hpp"
#include "Tensor2_symmetric_plus_Tensor2.hpp"
#include "Tensor2_symmetric_plus_Tensor2_symmetric.hpp"
#include "Tensor2_symmetric_plus_generic.hpp"
#include "Tensor2_symmetric_times_Tensor1.hpp"
#include "Tensor2_symmetric_times_Tensor2.hpp"
#include "Tensor2_symmetric_times_Tensor2_symmetric.hpp"
#include "Tensor2_symmetric_times_generic.hpp"
#include "dTensor2_symmetric.hpp"
#include "d_boundary_Tensor2_symmetric.hpp"
#include "d_one_sided_Tensor2_symmetric.hpp"
#include "ddTensor2_symmetric.hpp"
#include "dd_boundary_Tensor2_symmetric.hpp"
#include "diffusion_Tensor2_symmetric.hpp"
#include "generic_minus_Tensor2_symmetric.hpp"
#include "interpolate_Tensor2_symmetric.hpp"
#include "minus_Tensor2_symmetric.hpp"

namespace FTensor
{
  template <class A, class T, int Dim, char i, char j>
  class Tensor2_symmetric_Expr
  {
    A iter;

  public:
    Tensor2_symmetric_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2) const { return iter(N1, N2); }
  };

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  class Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i, j>
  {
    Tensor2_symmetric<A, Tensor_Dim> &iter;

  public:
    Tensor2_symmetric_Expr(Tensor2_symmetric<A, Tensor_Dim> &a) : iter(a) {}
    T &operator()(const int N1, const int N2) { return iter(N1, N2); }
    T operator()(const int N1, const int N2) const { return iter(N1, N2); }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result);

    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator=(const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T,
                                           Dim, i, j> &result);

    template <class B, class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator+=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result);

    template <class B, class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator-=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result);

    template <class B, class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator&=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result);

    /* This is for when the indices are switched (i,j) -> (j,i). */

    template <class B, class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator=(const Tensor2_symmetric_Expr<B, U, Dim, j, i> &result);

    template <class B, class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator+=(const Tensor2_symmetric_Expr<B, U, Dim, j, i> &result);

    template <class B, class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator-=(const Tensor2_symmetric_Expr<B, U, Dim, j, i> &result);

    /* Operations with just generics. */

    template <class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator=(const U &d);
    template <class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator+=(const U &d);
    template <class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator-=(const U &d);
    template <class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator*=(const U &d);
    template <class U>
    const Tensor2_symmetric_Expr<Tensor2_symmetric<A, Tensor_Dim>, T, Dim, i,
                                 j> &
    operator/=(const U &d);
  };

  /* Specialized for Dg_number_rhs_2 (Dg with the
     second index explicitly given). */

  template <class A, class T, int Dim, char i, char j, int N>
  class Tensor2_symmetric_Expr<Dg_number_rhs_2<A, T, N>, T, Dim, i, j>
  {
    A &iter;

  public:
    Tensor2_symmetric_Expr(A &a) : iter(a) {}
    T &operator()(const int N1, const int N2) { return iter(N1, N2, N); }
    T operator()(const int N1, const int N2) const { return iter(N1, N2, N); }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    const Tensor2_symmetric_Expr<Dg_number_rhs_2<A, T, N>, T, Dim, i, j> &
    operator=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result);

    const Tensor2_symmetric_Expr<Dg_number_rhs_2<A, T, N>, T, Dim, i, j> &
    operator=(const Tensor2_symmetric_Expr<Dg_number_rhs_2<A, T, N>, T, Dim, i,
                                           j> &result);
  };

  /* Specialized for Ddg_number_rhs_01 (Ddg with the
     first and second index explicitly given). */
  /*
  template <class A, class T, int Dim, char i, char j, int N0, int N1>
  class Tensor2_symmetric_Expr<Ddg_number_rhs_01<A, T, N0, N1>, T, Dim, i, j>
  {
    A &iter;

  public:
    Tensor2_symmetric_Expr(A &a) : iter(a) {}
    T &operator()(const int N2, const int N3) { return iter(N0, N1, N2, N3); }
    T operator()(const int N2, const int N3) const
    {
      return iter(N0, N1, N2, N3);
      }*/

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */
  /*
    template <class B, class U>
    const Tensor2_symmetric_Expr<Ddg_number_rhs_01<A, T, N0, N1>, T, Dim, i, j>
      &operator=(const Tensor2_symmetric_Expr<B, U, Dim, i, j> &result);

    const Tensor2_symmetric_Expr<Ddg_number_rhs_01<A, T, N0, N1>, T, Dim, i, j>
      &operator=(const Tensor2_symmetric_Expr<Ddg_number_rhs_01<A, T, N0, N1>,
                                              T, Dim, i, j> &result);
  };*/
}

#include "Tensor2_symmetric_Expr_equals.hpp"
