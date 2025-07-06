/* Declares a wrapper class for symmetric rank 3 Tensor expressions
   that are symmetric on the first two indices.  It is also used for
   Christof's but with some games played with the
   indices. */

#pragma once

#include "Dg_and_Dg.hpp"
#include "Dg_and_Tensor1.hpp"
#include "Dg_and_Tensor2_symmetric.hpp"
#include "Dg_divide_generic.hpp"
#include "Dg_minus_Dg.hpp"
#include "Dg_or_Dg.hpp"
#include "Dg_plus_Dg.hpp"
#include "Dg_times_Dg.hpp"
#include "Dg_times_Tensor1.hpp"
#include "Dg_times_Tensor2.hpp"
#include "Dg_times_Tensor2_symmetric.hpp"
#include "Dg_times_generic.hpp"
#include "minus_Dg.hpp"

namespace FTensor
{
  template <class A, class T, int Dim01, int Dim2, char i, char j, char k>
  class Dg_Expr
  {
    A iter;

  public:
    Dg_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3) const
    {
      return iter(N1, N2, N3);
    }
  };

  template <class A, class T, int Tensor_Dim01, int Tensor_Dim2, int Dim01,
            int Dim2, char i, char j, char k>
  class Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k>
  {
    Dg<A, Tensor_Dim01, Tensor_Dim2> &iter;

  public:
    Dg_Expr(Dg<A, Tensor_Dim01, Tensor_Dim2> &a) : iter(a) {}
    T &operator()(const int N1, const int N2, const int N3)
    {
      return iter(N1, N2, N3);
    }
    T operator()(const int N1, const int N2, const int N3) const
    {
      return iter(N1, N2, N3);
    }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k> &
    operator=(const Dg_Expr<B, U, Dim01, Dim2, i, j, k> &result);

    Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k> &
    operator=(const Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2,
                            i, j, k> &result);

    template <class U>
    Dg_Expr<Dg<A, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k> &
    operator=(const U &d);
  };

  /* I need a version for const and non-const Christof,
     otherwise it will use the default and the index order will get all
     messed up. The class A is either T or T*, depending on whether
     Christof has an array of T's or T *'s. */

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim12, int Dim12,
            int Dim0, char i, char j, char k>
  class Dg_Expr<const Christof<A, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0,
                i, j, k>
  {
    Christof<A, Tensor_Dim0, Tensor_Dim12> iter;

  public:
    Dg_Expr(const Christof<A, Tensor_Dim0, Tensor_Dim12> &a) : iter(a) {}

    /* Need to switch the index order because a christof tensor is just
       a dg tensor with the indices switched. */

    T operator()(const int N1, const int N2, const int N3) const
    {
      return iter(N3, N1, N2);
    }
  };

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim12, int Dim12,
            int Dim0, char i, char j, char k>
  class Dg_Expr<Christof<A, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k>
  {
    Christof<A, Tensor_Dim0, Tensor_Dim12> &iter;

  public:
    Dg_Expr(Christof<A, Tensor_Dim0, Tensor_Dim12> &a) : iter(a) {}

    /* Need to switch the index order because a christof tensor is just
       a dg tensor with the indices switched.  I have to explicitly
       declare the second operator= because otherwise the compiler will
       generate its own and not use the template code. */

    T operator()(const int N1, const int N2, const int N3) const
    {
      return iter(N3, N1, N2);
    }
    template <class B, class U>
    Dg_Expr<Christof<A, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k> &
    operator=(const Dg_Expr<B, U, Dim12, Dim0, i, j, k> &result);

    Dg_Expr<Christof<A, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k> &
    operator=(const Dg_Expr<Christof<A, Tensor_Dim0, Tensor_Dim12>, T, Dim12,
                            Dim0, i, j, k> &result);

    template <class U>
    Dg_Expr<Christof<A, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k> &
    operator=(const U &d);
  };

  /* Specialized for Ddg_number_rhs_0 (Ddg with the
     first index explicitly given). */

  template <class A, class T, int Dim23, int Dim1, char i, char j, char k,
            int N0>
  class Dg_Expr<Ddg_number_rhs_0<A, T, N0>, T, Dim23, Dim1, i, j, k>
  {
    A &iter;

  public:
    Dg_Expr(A &a) : iter(a) {}
    T &operator()(const int N1, const int N2, const int N3)
    {
      return iter(N0, N1, N2, N3);
    }
    T operator()(const int N1, const int N2, const int N3) const
    {
      return iter(N0, N1, N2, N3);
    }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    Dg_Expr<Ddg_number_rhs_0<A, T, N0>, T, Dim23, Dim1, i, j, k> &
    operator=(const Dg_Expr<B, U, Dim23, Dim1, i, j, k> &result);

    Dg_Expr<Ddg_number_rhs_0<A, T, N0>, T, Dim23, Dim1, i, j, k> &
    operator=(const Dg_Expr<Ddg_number_rhs_0<A, T, N0>, T, Dim23, Dim1, i, j,
                            k> &result);
  };
}

#include "Dg_Expr_equals.hpp"
