/* This is for expressions where a number is used for one or two
   slots, and an index for the others, yielding a Tensor1_Expr,
   Tensor2_Expr, or Tensor2_symmetric_Expr. */

#pragma once

namespace FTensor
{
  /* First or second slot. */

  template <class A, class T, int N> class Dg_number_0
  {
    A iterA;

  public:
    T operator()(const int N1, const int N2) const { return iterA(N, N1, N2); }
    Dg_number_0(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N> class Dg_number_rhs_0
  {};

  /* Third slot. */

  template <class A, class T, int N> class Dg_number_2
  {
    A iterA;

  public:
    T operator()(const int N1, const int N2) const { return iterA(N1, N2, N); }
    Dg_number_2(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N> class Dg_number_rhs_2
  {};

  /* Second and third slot (or first and third slot). */

  template <class A, class T, int N1, int N2> class Dg_number_12
  {
    A iterA;

  public:
    T operator()(const int N) const { return iterA(N, N1, N2); }
    Dg_number_12(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N1, int N2> class Dg_number_rhs_12
  {};

  /* First and second slot. */

  template <class A, class T, int N1, int N2> class Dg_number_01
  {
    A iterA;

  public:
    T operator()(const int N) const { return iterA(N1, N2, N); }
    Dg_number_01(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N1, int N2> class Dg_number_rhs_01
  {};
}
