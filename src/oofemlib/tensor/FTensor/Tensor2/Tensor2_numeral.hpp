/* This is for expressions where a number is used for one slot, and
   an index for another, yielding a Tensor1_Expr. */

#pragma once

namespace FTensor
{
  template <class A, class T> class Tensor2_numeral_1
  {
    A iterA;
    const int N;

  public:
    T operator()(const int N1) const { return iterA(N1, N); }
    Tensor2_numeral_1(A &a, const int NN) : iterA(a), N(NN) {}
  };

  template <class A, class T> class Tensor2_numeral_0
  {
    A iterA;
    const int N;

  public:
    T operator()(const int N1) const { return iterA(N, N1); }
    Tensor2_numeral_0(A &a, const int NN) : iterA(a), N(NN) {}
  };
}
