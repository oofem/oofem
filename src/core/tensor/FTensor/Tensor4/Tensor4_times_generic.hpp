/* Multiplies a Tensor4 with a generic, yielding a Tensor2. */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
  class Tensor4_times_generic
  {
    const Tensor4_Expr<A, T, Dim0, Dim1,  Dim2, Dim3, i, j, k, l> iterA;
    const U d;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) * d;
    }

    Tensor4_times_generic(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
                          const U &d0)
        : iterA(a), d(d0)
    {}
  };

  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
  Tensor4_Expr<Tensor4_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
               typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3,  i, j, k, l>
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a, const U &d0)
  {
    using TensorExpr = Tensor4_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i,
                        j, k, l>(TensorExpr(a, d0));
  }

  template <class A, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, char i, char j, char k, char l>
  Tensor4_Expr<Tensor4_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>,
               typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>
  operator*(const U &d0, const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    using TensorExpr = Tensor4_times_generic<A, T, U, Dim0, Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i,
                        j, k, l>(TensorExpr(a, d0));
  }
}
