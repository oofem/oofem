/* This is for expressions where a number is used for one or two
   slots, and an index for the others, yielding a Tensor1_Expr or
   Tensor2_Expr. */

#pragma once

namespace FTensor
{
  /* First slot. */

  template <class A, class T, int N> class Tensor3_number_0
  {
    A iterA;

  public:
    T operator()(const int N1, const int N2) const { return iterA(N, N1, N2); }
    Tensor3_number_0(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N> class Tensor3_number_rhs_0
  {};

  /* Second slot. */

  template <class A, class T, int N> class Tensor3_number_1
  {
    A iterA;

  public:
    T operator()(const int N1, const int N2) const { return iterA(N1, N, N2); }
    Tensor3_number_1(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N> class Tensor3_number_rhs_1
  {};

  /* Third slot. */

  template <class A, class T, int N> class Tensor3_number_2
  {
    A iterA;

  public:
    T operator()(const int N1, const int N2) const { return iterA(N1, N2, N); }
    Tensor3_number_2(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N> class Tensor3_number_rhs_2
  {};

  /* First and second slot. */

  template <class A, class T, int N1, int N2> class Tensor3_number_01
  {
    A iterA;

  public:
    T operator()(const int N) const { return iterA(N1, N2, N); }
    Tensor3_number_01(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N1, int N2> class Tensor3_number_rhs_01
  {};

  /* First and third slot. */

  template <class A, class T, int N1, int N2> class Tensor3_number_02
  {
    A iterA;

  public:
    T operator()(const int N) const { return iterA(N1, N, N2); }
    Tensor3_number_02(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N1, int N2> class Tensor3_number_rhs_02
  {};

  /* Second and third slot. */

  template <class A, class T, int N1, int N2> class Tensor3_number_12
  {
    A iterA;

  public:
    T operator()(const int N) const { return iterA(N, N1, N2); }
    Tensor3_number_12(const A &a) : iterA(a) {}
  };

  template <class A, class T, int N1, int N2> class Tensor3_number_rhs_12
  {};
}
