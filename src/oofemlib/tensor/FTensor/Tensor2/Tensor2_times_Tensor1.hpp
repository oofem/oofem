// This file has all of the declarations for expressions like
// Tensor2*Tensor1 and Tensor1*Tensor2, yielding a Tensor1 or Tensor3.

#pragma once

#include "../permute.hpp"

namespace FTensor
{
  // A(i,j) * B(k) single contraction
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k, int DimA, int DimX, char a, char x>
  class Tensor2_times_Tensor1_single
  {
    const Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;
    const Tensor1_Expr<B, U, Dim2, k> iterB;

  public:
    Tensor2_times_Tensor1_single(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &iter_a,
                                 const Tensor1_Expr<B, U, Dim2, k> &iter_b)
        : iterA(iter_a), iterB(iter_b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      typename promote<T, U>::V result(0);
      for(int xx = 0; xx < DimX; ++xx)
        {
          result += iterB(xx)
                    * Permutation2<DimA, DimX, a, x>().eval(iterA, N1, xx);
        }
      return result;
    }
  };

  // A(i,j)*B(j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  auto operator*(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                 const Tensor1_Expr<B, U, Dim1, j> &b)
  {
    using TensorExpr
      = Tensor2_times_Tensor1_single<A, B, T, U, Dim0, Dim1, Dim1, i, j, j,
                                     Dim0, Dim1, i, j>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  // B(j)*A(i,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  auto operator*(const Tensor1_Expr<B, U, Dim1, j> &b,
                 const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a)
  {
    return a * b;
  }

  // A(i,j)*B(i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  auto operator*(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                 const Tensor1_Expr<B, U, Dim0, i> &b)
  {
    using TensorExpr
      = Tensor2_times_Tensor1_single<A, B, T, U, Dim0, Dim1, Dim0, i, j, i,
                                     Dim1, Dim0, j, i>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, j>(
      TensorExpr(a, b));
  }

  // B(i)*A(i,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, char i,
            char j>
  auto operator*(const Tensor1_Expr<B, U, Dim0, i> &b,
                 const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a)
  {
    return a * b;
  }

  /* A(i,j)*B(k) -> Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  class Tensor2_times_Tensor1
  {
    const Tensor2_Expr<A, T, Dim0, Dim1, i, j> iterA;
    const Tensor1_Expr<B, U, Dim2, k> iterB;

  public:
    Tensor2_times_Tensor1(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
                          const Tensor1_Expr<B, U, Dim2, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return iterA(N1, N2) * iterB(N3);
    }
  };

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  Tensor3_Expr<Tensor2_times_Tensor1<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim0, Dim1, Dim2, i, j, k>
  operator*(const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a,
            const Tensor1_Expr<B, U, Dim2, k> &b)
  {
    using TensorExpr
      = Tensor2_times_Tensor1<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, i, j, k>(TensorExpr(a, b));
  }

  /* B(k)*A(i,j) -> Tensor3 */

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            char i, char j, char k>
  Tensor3_Expr<Tensor2_times_Tensor1<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim0, Dim1, Dim2, i, j, k>
  operator*(const Tensor1_Expr<B, U, Dim2, k> &b,
            const Tensor2_Expr<A, T, Dim0, Dim1, i, j> &a)
  {
    using TensorExpr
      = Tensor2_times_Tensor1<A, B, T, U, Dim0, Dim1, Dim2, i, j, k>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, i, j, k>(TensorExpr(a, b));
  }
}
