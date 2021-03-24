/* Declares a wrapper class for rank 2 Tensor expressions.  I
   specialize it for when I wrap a simple Tensor2 or Tensor2_ptr so
   that it has a reference to the Tensor2(_ptr) and not a copy.
   Otherwise assignment wouldn't work. */

#pragma once

#include "Tensor2_and_Tensor1.hpp"
#include "Tensor2_carat_Tensor2.hpp"
#include "Tensor2_divide_generic.hpp"
#include "Tensor2_minus_Tensor2.hpp"
#include "Tensor2_or_Tensor2.hpp"
#include "Tensor2_plus_Tensor2.hpp"
#include "Tensor2_times_Tensor1.hpp"
#include "Tensor2_times_Tensor2.hpp"
#include "Tensor2_times_generic.hpp"
#include "Tensor2_transform.hpp"
#include "conj_Tensor2.hpp"
#include "minus_Tensor2.hpp"

#include "../permute.hpp"

namespace FTensor
{
  template <class A, class T, int Dim0, int Dim1, char i, char j>
  class Tensor2_Expr
  {
    A iter;

  public:
    Tensor2_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2) const { return iter(N1, N2); }
  };

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Dim0,
            int Dim1, char i, char j>
  class Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1>, T, Dim0, Dim1, i, j>
  {
    Tensor2<A, Tensor_Dim0, Tensor_Dim1> &iter;

  public:
    Tensor2_Expr(Tensor2<A, Tensor_Dim0, Tensor_Dim1> &a) : iter(a) {}
    T &operator()(const int N1, const int N2) { return iter(N1, N2); }
    T operator()(const int N1, const int N2) const { return iter(N1, N2); }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U, int Dim1_0, int Dim1_1, char i_1, char j_1>
    auto &equals(const Tensor2_Expr<B, U, Dim1_0, Dim1_1, i_1, j_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          {
            iter(ii, jj) = permute(*this, rhs, ii, jj);
          }
      return *this;
    }

    template <class B, class U, int Dim1_0, int Dim1_1, char i_1, char j_1>
    auto &
    operator=(const FTensor::Tensor2_Expr<B, U, Dim1_0, Dim1_1, i_1, j_1> &rhs)
    {
      return equals(rhs);
    }

    auto &operator=(const Tensor2_Expr<Tensor2<A, Tensor_Dim0, Tensor_Dim1>, T,
                                       Dim0, Dim1, i, j> &rhs)
    {
      return equals(rhs);
    }

    template <class B, class U, int Dim1_0, int Dim1_1, char i_1, char j_1>
    auto &operator+=(const Tensor2_Expr<B, U, Dim1_0, Dim1_1, i_1, j_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          {
            iter(ii, jj) += permute(*this, rhs, ii, jj);
          }
      return *this;
    }

    template <class B, class U, int Dim1_0, int Dim1_1, char i_1, char j_1>
    auto &operator-=(const Tensor2_Expr<B, U, Dim1_0, Dim1_1, i_1, j_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          {
            iter(ii, jj) -= permute(*this, rhs, ii, jj);
          }
      return *this;
    }

    /* This is for int's, double's, etc. */

    template <class U> auto &operator=(const U &u)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          {
            iter(ii, jj) = u;
          }
      return *this;
    }

    template <class U> auto &operator+=(const U &u)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          {
            iter(ii, jj) += u;
          }
      return *this;
    }

    template <class U> auto &operator-=(const U &u)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          {
            iter(ii, jj) -= u;
          }
      return *this;
    }

    template <class U> auto &operator*=(const U &u)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          {
            iter(ii, jj) *= u;
          }
      return *this;
    }

    template <class U> auto &operator/=(const U &u)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          {
            iter(ii, jj) /= u;
          }
      return *this;
    }
  };

  /* Specialized for Dg_number_rhs_0 (Dg with the
     first or second index explicitly given). */

  template <class A, class T, int Dim0, int Dim1, char i, char j, int N>
  class Tensor2_Expr<Dg_number_rhs_0<A, T, N>, T, Dim0, Dim1, i, j>
  {
    A &iter;

  public:
    Tensor2_Expr(A &a) : iter(a) {}
    T &operator()(const int N1, const int N2) { return iter(N, N1, N2); }
    T operator()(const int N1, const int N2) const { return iter(N, N1, N2); }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U, int Dim1_0, int Dim1_1, char i_1, char j_1>
    auto &operator=(const Tensor2_Expr<B, U, Dim1_0, Dim1_1, i_1, j_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          {
            iter(ii, jj) = permute(*this, rhs, ii, jj);
          }
      return *this;
    }

    auto &operator=(
      const Tensor2_Expr<Dg_number_rhs_0<A, T, N>, T, Dim0, Dim1, i, j> &result)
    {
      return operator=<Dg_number_rhs_0<A, T, N>, T>(result);
    }
  };
}
