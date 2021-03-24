/* This is for expressions where a number is used for one slot, and
   an index for the other, yielding a Tensor2_symmetric_Expr. */

#pragma once

namespace FTensor
{
  template <class A, class T, int N> class Christof_number_0
  {
    const A iterA;

  public:
    T operator()(const int N1, const int N2) const { return iterA(N, N1, N2); }
    Christof_number_0(const A &a) : iterA(a) {}
  };
}
