#pragma once

#include "Permutation2.hpp"
#include "Permutation3.hpp"
#include "Permutation4.hpp"

namespace FTensor
{
  template <class A, class T, class B, class U, int Dim0_0, int Dim0_1,
            int Dim1_0, int Dim1_1, char i0, char j0, char i1, char j1>
  U permute(const Tensor2_Expr<A, T, Dim0_0, Dim0_1, i0, j0> &,
            const Tensor2_Expr<B, U, Dim1_0, Dim1_1, i1, j1> &rhs,
            const int N0, const int N1)
  {
    return Permutation2<Dim0_0, Dim0_1, i0, j0>().eval(rhs, N0, N1);
  }

  template <class A, class T, class B, class U, int Dim0_0, int Dim0_1,
            int Dim0_2, int Dim1_0, int Dim1_1, int Dim1_2, char i0, char j0,
            char k0, char i1, char j1, char k1>
  U permute(const Tensor3_Expr<A, T, Dim0_0, Dim0_1, Dim0_2, i0, j0, k0> &,
            const Tensor3_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, i1, j1, k1> &rhs,
            const int N0, const int N1, const int N2)
  {
    return Permutation3<Dim0_0, Dim0_1, Dim0_2, i0, j0, k0>().eval(rhs, N0, N1,
                                                                   N2);
  }

  template <class A, class T, class B, class U, int Dim0_0, int Dim0_1,
            int Dim0_2, int Dim0_3, int Dim1_0, int Dim1_1, int Dim1_2,
            int Dim1_3, char i0, char j0, char k0, char l0, char i1, char j1,
            char k1, char l1>
  U permute(
    const Tensor4_Expr<A, T, Dim0_0, Dim0_1, Dim0_2, Dim0_3, i0, j0, k0, l0> &,
    const Tensor4_Expr<B, U, Dim1_0, Dim1_1, Dim1_2, Dim1_3, i1, j1, k1, l1>
      &rhs,
    const int N0, const int N1, const int N2, const int N3)
  {
    return Permutation4<Dim0_0, Dim0_1, Dim0_2, Dim0_3, i0, j0, k0, l0>().eval(
      rhs, N0, N1, N2, N3);
  }
}
