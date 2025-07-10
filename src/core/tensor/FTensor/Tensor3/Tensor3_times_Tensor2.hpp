// This file has all of the declarations for expressions like Tensor3*Tensor2
// and Tensor2*Tensor3, yielding a Tensor3 or Tensor1.

#pragma once

#include "../permute.hpp"

namespace FTensor
{
  // A(i,j,k) * B(l,m) double contraction
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m,
            int DimA, int DimX, int DimY, char a, char x, char y>
  class Tensor3_times_Tensor2_double
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    Tensor2_Expr<B, U, Dim3, Dim4, l, m> iterB;

  public:
    Tensor3_times_Tensor2_double(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &Itera,
      const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &Iterb)
        : iterA(Itera), iterB(Iterb)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      typename promote<T, U>::V result(0);

      for(int xx = 0; xx < DimX; ++xx)
        for(int yy = 0; yy < DimY; ++yy)
          {
            result += Permutation3<DimA, DimX, DimY, a, x, y>().eval(iterA, N1,
                                                                     xx, yy)
                      * Permutation2<DimX, DimY, x, y>().eval(iterB, xx, yy);
          }
      return result;
    }
  };

  // A(i,j,k)*B(j,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                 const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim1, Dim2,
                                     i, j, k, j, k, Dim0, Dim1, Dim2, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  // B(j,k)*A(i,j,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b,
                 const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    return a * b;
  }

  // A(i,j,k)*B(k,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                 const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim2, Dim1,
                                     i, j, k, k, j, Dim0, Dim1, Dim2, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  // B(k,j)*A(i,j,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b,
                 const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    return a * b;
  }

  // A(i,j,k)*B(i,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                 const Tensor2_Expr<B, U, Dim0, Dim2, i, k> &b)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim0, Dim2,
                                     i, j, k, i, k, Dim1, Dim0, Dim2, j, i, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, j>(
      TensorExpr(a, b));
  }

  // A(i,j,k)*B(i,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor2_Expr<B, U, Dim0, Dim2, i, k> &b,
                 const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    return a * b;
  }

  // A(i,j,k)*B(k,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                 const Tensor2_Expr<B, U, Dim2, Dim0, k, i> &b)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim2, Dim0,
                                     i, j, k, k, i, Dim1, Dim0, Dim2, j, i, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, j>(
      TensorExpr(a, b));
  }

  // B(k,i)*A(i,j,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor2_Expr<B, U, Dim2, Dim0, k, i> &b,
                 const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    return a * b;
  }

  // A(i,j,k)*B(i,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                 const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim0, Dim1,
                                     i, j, k, i, j, Dim2, Dim0, Dim1, k, i, j>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, k>(
      TensorExpr(a, b));
  }

  // B(i,j)*A(i,j,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b,
                 const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    return a * b;
  }

  // A(i,j,k)*B(j,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
                 const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b)
  {
    using TensorExpr
      = Tensor3_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim1, Dim0,
                                     i, j, k, j, i, Dim2, Dim0, Dim1, k, i, j>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, k>(
      TensorExpr(a, b));
  }

  // B(j,i)*A(i,j,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  auto operator*(const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b,
                 const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    return a * b;
  }

  /* A(i,j,k)*B(k,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_2_01
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    Tensor2_Expr<B, U, Dim2, Dim3, k, l> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, Current_Dim - 1) * iterB(Current_Dim - 1, N3)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(0, N3);
    }

  public:
    Tensor3_times_Tensor2_2_01(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_2_01<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j, l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b)
  {
    using TensorExpr = Tensor3_times_Tensor2_2_01<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim3, i, j, l>(TensorExpr(a, b));
  }

  /* B(k,l)*A(i,j,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_2_01<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j, l>
  operator*(const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    using TensorExpr = Tensor3_times_Tensor2_2_01<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim3, i, j, l>(TensorExpr(a, b));
  }

  /* A(i,j,k)*B(l,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_2_10
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> iterA;
    Tensor2_Expr<B, U, Dim3, Dim2, l, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, Current_Dim - 1) * iterB(N3, Current_Dim - 1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(N3, 0);
    }

  public:
    Tensor3_times_Tensor2_2_10(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
      const Tensor2_Expr<B, U, Dim2, Dim3, l, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_2_10<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j, l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a,
            const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b)
  {
    using TensorExpr = Tensor3_times_Tensor2_2_10<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim3, i, j, l>(TensorExpr(a, b));
  }

  /* B(l,k)*A(i,j,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_2_10<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j, l>
  operator*(const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim2, i, j, k> &a)
  {
    using TensorExpr = Tensor3_times_Tensor2_2_10<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim3, i, j, l>(TensorExpr(a, b));
  }

  /* A(i,j,l)*B(j,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_1_01
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> iterA;
    Tensor2_Expr<B, U, Dim1, Dim2, j, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1, N3) * iterB(Current_Dim - 1, N2)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, 0, N3) * iterB(0, N2);
    }

  public:
    Tensor3_times_Tensor2_1_01(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a,
      const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_1_01<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k, l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a,
            const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
  {
    using TensorExpr = Tensor3_times_Tensor2_1_01<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2,
                        Dim3, i, k, l>(TensorExpr(a, b));
  }

  /* B(j,k)*A(i,j,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_1_01<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k, l>
  operator*(const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a)
  {
    using TensorExpr = Tensor3_times_Tensor2_1_01<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2,
                        Dim3, i, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,l)*B(k,j)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_1_10
  {
    Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> iterA;
    Tensor2_Expr<B, U, Dim2, Dim1, k, j> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1, N3) * iterB(N2, Current_Dim - 1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, 0, N3) * iterB(N2, 0);
    }

  public:
    Tensor3_times_Tensor2_1_10(
      const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a,
      const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_1_10<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k, l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a,
            const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
  {
    using TensorExpr = Tensor3_times_Tensor2_1_10<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2,
                        Dim3, i, k, l>(TensorExpr(a, b));
  }

  /* B(k,j)*A(i,j,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_1_10<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k, l>
  operator*(const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b,
            const Tensor3_Expr<A, T, Dim0, Dim1, Dim3, i, j, l> &a)
  {
    using TensorExpr = Tensor3_times_Tensor2_1_10<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2,
                        Dim3, i, k, l>(TensorExpr(a, b));
  }

  /* A(i,k,l)*B(i,j)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_0_01
  {
    Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> iterA;
    Tensor2_Expr<B, U, Dim0, Dim1, i, j> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(Current_Dim - 1, N2, N3) * iterB(Current_Dim - 1, N1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(0, N2, N3) * iterB(0, N1);
    }

  public:
    Tensor3_times_Tensor2_0_01(
      const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a,
      const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim0>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_0_01<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim1, Dim2, Dim3, j, k, l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a,
            const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b)
  {
    using TensorExpr = Tensor3_times_Tensor2_0_01<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2,
                        Dim3, j, k, l>(TensorExpr(a, b));
  }

  /* B(i,j)*A(i,k,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_0_01<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim1, Dim2, Dim3, j, k, l>
  operator*(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b,
            const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a)
  {
    using TensorExpr = Tensor3_times_Tensor2_0_01<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2,
                        Dim3, j, k, l>(TensorExpr(a, b));
  }

  /* A(i,k,l)*B(j,i)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  class Tensor3_times_Tensor2_0_10
  {
    Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> iterA;
    Tensor2_Expr<B, U, Dim1, Dim0, j, i> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(Current_Dim - 1, N2, N3) * iterB(N1, Current_Dim - 1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(0, N2, N3) * iterB(N1, 0);
    }

  public:
    Tensor3_times_Tensor2_0_10(
      const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a,
      const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim0>());
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_0_10<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim1, Dim2, Dim3, j, k, l>
  operator*(const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a,
            const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b)
  {
    using TensorExpr = Tensor3_times_Tensor2_0_10<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2,
                        Dim3, j, k, l>(TensorExpr(a, b));
  }

  /* B(j,i)*A(i,k,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  Tensor3_Expr<
    Tensor3_times_Tensor2_0_10<A, B, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
    typename promote<T, U>::V, Dim1, Dim2, Dim3, j, k, l>
  operator*(const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b,
            const Tensor3_Expr<A, T, Dim0, Dim2, Dim3, i, k, l> &a)
  {
    using TensorExpr = Tensor3_times_Tensor2_0_10<A, B, T, U, Dim0, Dim1, Dim2,
                                                  Dim3, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2,
                        Dim3, j, k, l>(TensorExpr(a, b));
  }
}
