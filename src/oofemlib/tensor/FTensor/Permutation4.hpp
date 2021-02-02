#pragma once

#include "error_when_instantiated.hpp"

namespace FTensor
{
  template <int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k,
            char l>
  struct Permutation4
  {
    Permutation4() = default;

    // (i,...)
    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N0, N1, N2, N3);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim0, Dim1, Dim3, Dim2, i, j, l, k> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N0, N1, N3, N2);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim0, Dim2, Dim1, Dim3, i, k, j, l> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N0, N2, N1, N3);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim0, Dim2, Dim3, Dim1, i, k, l, j> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N0, N2, N3, N1);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim0, Dim3, Dim1, Dim2, i, l, j, k> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N0, N3, N1, N2);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim0, Dim3, Dim2, Dim1, i, l, k, j> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N0, N3, N2, N1);
    }

    // (j,...)
    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim1, Dim0, Dim2, Dim3, j, i, k, l> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N1, N0, N2, N3);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim1, Dim0, Dim3, Dim2, j, i, l, k> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N1, N0, N3, N2);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim1, Dim2, Dim0, Dim3, j, k, i, l> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N1, N2, N0, N3);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim1, Dim2, Dim3, Dim0, j, k, l, i> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N1, N2, N3, N0);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim1, Dim3, Dim0, Dim2, j, l, i, k> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N1, N3, N0, N2);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim1, Dim3, Dim2, Dim0, j, l, k, i> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N1, N3, N2, N0);
    }

    // (k,...)
    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim2, Dim1, Dim0, Dim3, k, j, i, l> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N2, N1, N0, N3);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim2, Dim1, Dim3, Dim0, k, j, l, i> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N2, N1, N3, N0);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim2, Dim0, Dim1, Dim3, k, i, j, l> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N2, N0, N1, N3);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim2, Dim0, Dim3, Dim1, k, i, l, j> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N2, N0, N3, N1);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim2, Dim3, Dim1, Dim0, k, l, j, i> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N2, N3, N1, N0);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim2, Dim3, Dim0, Dim1, k, l, i, j> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N2, N3, N0, N1);
    }

    // (l,...)
    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim3, Dim1, Dim2, Dim0, l, j, k, i> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N3, N1, N2, N0);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim3, Dim1, Dim0, Dim2, l, j, i, k> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N3, N1, N0, N2);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim3, Dim2, Dim1, Dim0, l, k, j, i> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N3, N2, N1, N0);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim3, Dim2, Dim0, Dim1, l, k, i, j> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N3, N2, N0, N1);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim3, Dim0, Dim1, Dim2, l, i, j, k> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N3, N0, N1, N2);
    }

    template <class B, class U>
    U eval(const Tensor4_Expr<B, U, Dim3, Dim0, Dim2, Dim1, l, i, k, j> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      return rhs(N3, N0, N2, N1);
    }

    // Catch-all version for incompatible indices
    template <class B, class U, int Dim1_0, int Dim1_1, int Dim1_2, int Dim1_3,
              char i1, char j1, char k1, char l1>
    U eval(const Tensor4_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, i1, j1, k1,
                              l1> &rhs,
           const int N0, const int N1, const int N2, const int N3)
    {
      static_assert(error_when_instantiated<B>(), "Incompatible indices");
      return rhs(0, 0, 0, 0);
    }
  };
}
