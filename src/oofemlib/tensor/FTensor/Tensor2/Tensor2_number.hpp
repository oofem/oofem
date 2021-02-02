/* This is for expressions where a number is used for one slot, and
   an index for another, yielding a Tensor1_Expr. */

#pragma once

namespace FTensor
{
  template <class A, class T, int N> class Tensor2_number_1
  {
    A iterA;

  public:
    T operator()(const int N1) const { return iterA(N1, N); }
    Tensor2_number_1(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N> class Tensor2_number_0
  {
    A iterA;

  public:
    T operator()(const int N1) const { return iterA(N, N1); }
    Tensor2_number_0(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N> class Tensor2_number_rhs_0
  {};

  template <class A, class T, int N> class Tensor2_number_rhs_1
  {};
}
