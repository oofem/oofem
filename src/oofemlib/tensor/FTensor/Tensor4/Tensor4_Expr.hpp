/* Declare a wrapper class for generic rank 4 Tensor expressions. */

#pragma once

#include "Tensor4_minus_Tensor4.hpp"
#include "Tensor4_plus_Tensor4.hpp"
#include "Tensor4_times_Tensor2_single.hpp"
#include "Tensor4_times_Tensor2_double.hpp"
#include "Tensor4_times_Tensor2_symmetric.hpp"
#include "Tensor4_times_Tensor3_triple.hpp"
#include "Tensor4_times_generic.hpp"

#include "../permute.hpp"

namespace FTensor
{
  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, char i,
            char j, char k, char l>
  class Tensor4_Expr
  {
    A iter;

  public:
    Tensor4_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iter(N1, N2, N3, N4);
    }
  };

  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, char i,
            char j, char k, char l>
  class Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0, Dim1, Dim2,
                     Dim3, i, j, k, l>
  {
    Tensor4<A, Dim0, Dim1, Dim2, Dim3> &iter;

  public:
    Tensor4_Expr(Tensor4<A, Dim0, Dim1, Dim2, Dim3> &a) : iter(a) {}
    T &operator()(const int N1, const int N2, const int N3, const int N4)
    {
      return iter(N1, N2, N3, N4);
    }
    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iter(N1, N2, N3, N4);
    }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3,
              char i_1, char j_1, char k_1, char l_1>
    auto &equals(const Tensor4_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, i_1,
                                    j_1, k_1, l_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            for(int ll = 0; ll < Dim3; ++ll)
              {
                iter(ii, jj, kk, ll) = permute(*this, rhs, ii, jj, kk, ll);
              }
      return *this;
    }

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3,
              char i_1, char j_1, char k_1, char l_1>
    auto &operator=(const Tensor4_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3,
                                       i_1, j_1, k_1, l_1> &rhs)
    {
      return equals(rhs);
    }

    auto &
    operator=(const Tensor4_Expr<Tensor4<A, Dim0, Dim1, Dim2, Dim3>, T, Dim0,
                                 Dim1, Dim2, Dim3, i, j, k, l> &rhs)
    {
      return equals(rhs);
    }

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3,
              char i_1, char j_1, char k_1, char l_1>
    auto &operator+=(const Tensor4_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3,
                                        i_1, j_1, k_1, l_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            for(int ll = 0; ll < Dim3; ++ll)
              {
                iter(ii, jj, kk, ll) += permute(*this, rhs, ii, jj, kk, ll);
              }
      return *this;
    }

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3,
              char i_1, char j_1, char k_1, char l_1>
    auto &operator-=(const Tensor4_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3,
                                        i_1, j_1, k_1, l_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            for(int ll = 0; ll < Dim3; ++ll)
              {
                iter(ii, jj, kk, ll) += permute(*this, rhs, ii, jj, kk, ll);
              }
      return *this;
    }
  };
}
