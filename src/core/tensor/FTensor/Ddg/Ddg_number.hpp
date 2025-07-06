/* This is for expressions where a number is used for one or more
   slots, and an index for the others, yielding a
   Tensor2_symmetric_Expr or Dg_Expr. */

#pragma once

namespace FTensor
{
  /* First and second slot. */

  template <class A, class T, int N1, int N2> class Ddg_number_01
  {
    A iterA;

  public:
    T operator()(const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4);
    }
    Ddg_number_01(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N0, int N1> class Ddg_number_rhs_01
  {};

  /* Only first slot */

  template <class A, class T, int N1> class Ddg_number_0
  {
    A iterA;

  public:
    T operator()(const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4);
    }
    Ddg_number_0(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N0> class Ddg_number_rhs_0
  {};
}
