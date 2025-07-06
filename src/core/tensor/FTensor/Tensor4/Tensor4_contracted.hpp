/* A helper class that automatically returns the internal contraction
   of a rank 4 tensor. */

#pragma once

namespace FTensor
{
  template <class A, class T, int Dim> class Tensor4_contracted_23
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(N0, N1, ii, ii);
        }
      return result;
    }
    Tensor4_contracted_23(const A &a) : iterA(a) {}
  };

  template <class A, class T, int Dim> class Tensor4_contracted_13
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(N0, ii, N1, ii);
        }
      return result;
    }
    Tensor4_contracted_13(const A &a) : iterA(a) {}
  };

  template <class A, class T, int Dim> class Tensor4_contracted_12
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(N0, ii, ii, N1);
        }
      return result;
    }
    Tensor4_contracted_12(const A &a) : iterA(a) {}
  };

  template <class A, class T, int Dim> class Tensor4_contracted_03
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(ii, N0, N1, ii);
        }
      return result;
    }
    Tensor4_contracted_03(const A &a) : iterA(a) {}
  };

  template <class A, class T, int Dim> class Tensor4_contracted_02
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(ii, N0, ii, N1);
        }
      return result;
    }
    Tensor4_contracted_02(const A &a) : iterA(a) {}
  };

  template <class A, class T, int Dim> class Tensor4_contracted_01
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(ii, ii, N0, N1);
        }
      return result;
    }
    Tensor4_contracted_01(const A &a) : iterA(a) {}
  };
}
