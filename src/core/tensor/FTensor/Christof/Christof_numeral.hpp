/* This is for expressions where an int is used for one slot, and an
   Index for the others, yielding a Tensor2_symmetric_Expr or
   Tensor2_Expr. */

#pragma once

namespace FTensor
{
  template <class A, class T> class Christof_numeral_0
  {
    const A iterA;
    const int N;

  public:
    T operator()(const int N1, const int N2) const { return iterA(N, N1, N2); }
    Christof_numeral_0(const A &a, const int NN) : iterA(a), N(NN) {}
  };

  /* We can use this version for when there is an int in either of the
     last slots, because Christof is symmetric. */

  template <class A, class T> class Christof_numeral_1
  {
    const A iterA;
    const int N;

  public:
    T operator()(const int N1, const int N2) const { return iterA(N1, N, N2); }
    Christof_numeral_1(const A &a, const int NN) : iterA(a), N(NN) {}
  };
}
