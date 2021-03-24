#pragma once

#include "error_when_instantiated.hpp"

namespace FTensor
{
  template <int Dim0, int Dim1, int Dim2, char i, char j, char k>
  struct Permutation3
  {
    Permutation3() = default;

    template <class B, class U>
    U eval(const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &rhs,
           const int N0, const int N1, const int N2)
    {
      return rhs(N0, N1, N2);
    }

    template <class B, class U>
    U eval(const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &rhs,
           const int N0, const int N1, const int N2)
    {
      return rhs(N0, N2, N1);
    }

    template <class B, class U>
    U eval(const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &rhs,
           const int N0, const int N1, const int N2)
    {
      return rhs(N1, N0, N2);
    }

    template <class B, class U>
    U eval(const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &rhs,
           const int N0, const int N1, const int N2)
    {
      return rhs(N1, N2, N0);
    }

    template <class B, class U>
    U eval(const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &rhs,
           const int N0, const int N1, const int N2)
    {
      return rhs(N2, N0, N1);
    }

    template <class B, class U>
    U eval(const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &rhs,
           const int N0, const int N1, const int N2)
    {
      return rhs(N2, N1, N0);
    }

    // Catch-all version for incompatible indices
    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, char i1,
              char j1, char k1>
    U eval(const Tensor3_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, i1, j1, k1> &rhs,
           const int N0, const int N1, const int N2)
    {
      static_assert(error_when_instantiated<B>(), "Incompatible indices");
      return rhs(0, 0, 0);
    }
  };
}
