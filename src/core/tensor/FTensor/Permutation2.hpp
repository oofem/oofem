#pragma once

#include "error_when_instantiated.hpp"

namespace FTensor
{
  template <int Dim0, int Dim1, char i, char j> struct Permutation2
  {
    Permutation2() = default;

    template <class B, class U>
    U eval(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &rhs, const int N0,
           const int N1)
    {
      return rhs(N0, N1);
    }

    template <class B, class U>
    U eval(const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &rhs, const int N0,
           const int N1)
    {
      return rhs(N1, N0);
    }

    template <class B, class U, int Dim1_0, int Dim1_1, char i1, char j1>
    U eval(const Tensor2_Expr<B, U, Dim1_0, Dim1_1, i1, j1> &rhs, const int N0,
           const int N1)
    {
      static_assert(error_when_instantiated<B>(), "Incompatible indices");
      return rhs(0, 0);
    }
  };
}
