/* This is for expressions where a number is used for one or two
   slots, and an index for the others, yielding a Tensor1_Expr,
   Tensor2_Expr, or Tensor2_symmetric_Expr. */

#pragma once

namespace FTensor
{
  /* First or second slot. */

  template <class A, class T> class Dg_numeral_0
  {
    A iterA;
    int N;

  public:
    T operator()(const int N1, const int N2) const { return iterA(N, N1, N2); }
    Dg_numeral_0(const A &a, const int NN) : iterA(a), N(NN) {}
  };

  /* Third slot. */

  template <class A, class T> class Dg_numeral_2
  {
    A iterA;
    int N;

  public:
    T operator()(const int N1, const int N2) const { return iterA(N1, N2, N); }
    Dg_numeral_2(const A &a, const int NN) : iterA(a), N(NN) {}
  };

  /* Second and third slot (or first and third slot). */

  template <class A, class T> class Dg_numeral_12
  {
    A iterA;
    int N1, N2;

  public:
    T operator()(const int N) const { return iterA(N, N1, N2); }
    Dg_numeral_12(const A &a, const int NN1, const int NN2)
        : iterA(a), N1(NN1), N2(NN2)
    {}
  };

  /* First and second slot. */

  template <class A, class T> class Dg_numeral_01
  {
    A iterA;
    int N1, N2;

  public:
    T operator()(const int N) const { return iterA(N1, N2, N); }
    Dg_numeral_01(const A &a, const int NN1, const int NN2)
        : iterA(a), N1(NN1), N2(NN2)
    {}
  };
}
