/* Declare a wrapper class for generic rank 3 Tensor expressions.
   There isn't a Tensor3 class yet.  I only use Tensor3_Expr as an
   intermediate expression which immediately get contracted with
   something to make a Tensor2 or Tensor1. */

#pragma once

#include "Tensor3_minus_Tensor3.hpp"
#include "Tensor3_or_Tensor3.hpp"
#include "Tensor3_plus_Tensor3.hpp"
#include "Tensor3_times_Dg.hpp"
#include "Tensor3_times_Tensor1.hpp"
#include "Tensor3_times_Tensor2.hpp"
#include "Tensor3_times_Tensor2_symmetric.hpp"
#include "Tensor3_times_Tensor3.hpp"
#include "Tensor3_times_generic.hpp"

namespace FTensor
{
  template <class A, class T, int Dim0, int Dim1, int Dim2, char i, char j,
            char k>
  class Tensor3_Expr
  {
    A iter;

  public:
    Tensor3_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3) const
    {
      return iter(N1, N2, N3);
    }
  };

  template <class A, class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Dim0, int Dim1, int Dim2, char i, char j, char k>
  class Tensor3_Expr<Tensor3<A, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T,
                     Dim0, Dim1, Dim2, i, j, k>
  {
    Tensor3<A, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &iter;

  public:
    Tensor3_Expr(Tensor3<A, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &a)
        : iter(a)
    {}
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

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, char i_1,
              char j_1, char k_1>
    auto &
    equals(const Tensor3_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, i_1, j_1, k_1> &rhs)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            {
              iter(ii, jj, kk) = permute(*this, rhs, ii, jj, kk);
            }
      return *this;
    }

    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, char i_1,
              char j_1, char k_1>
    auto &operator=(
      const Tensor3_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, i_1, j_1, k_1> &rhs)
    {
      return equals(rhs);
    }

    auto &operator=(
      const Tensor3_Expr<Tensor3<A, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T,
                         Dim0, Dim1, Dim2, i, j, k> &rhs)
    {
      return equals(rhs);
    }

    template <class U>
    Tensor3_Expr<Tensor3<A, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim0,
                 Dim1, Dim2, i, j, k> &
    operator=(const U &u)
    {
      for(int ii = 0; ii < Dim0; ++ii)
        for(int jj = 0; jj < Dim1; ++jj)
          for(int kk = 0; kk < Dim2; ++kk)
            {
              iter(ii, jj, kk) = u;
            }
      return *this;
    }
  };
}
