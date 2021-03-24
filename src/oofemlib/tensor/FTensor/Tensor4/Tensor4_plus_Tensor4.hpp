// Adds Tensor4+Tensor4 -> Tensor4

#pragma once

namespace FTensor
{
  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim2_0, int Dim3_0, int Dim0_1, int Dim1_1, int Dim2_1,
            int Dim3_1, char i0, char j0, char k0, char l0, char i1, char j1,
            char k1, char l1>
  class Tensor4_plus_Tensor4
  {
    Tensor4_Expr<A, T, Dim0_0, Dim1_0, Dim2_0, Dim3_0, i0, j0, k0, l0> iterA;
    Tensor4_Expr<B, U, Dim0_1, Dim1_1, Dim2_1, Dim3_1, i1, j1, k1, l1> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) + permute(iterA, iterB, N1, N2, N3, N4);
    }

    Tensor4_plus_Tensor4(const Tensor4_Expr<A, T, Dim0_0, Dim1_0, Dim2_0,
                                            Dim3_0, i0, j0, k0, l0> &a,
                         const Tensor4_Expr<B, U, Dim0_1, Dim1_1, Dim2_1,
                                            Dim3_1, i1, j1, k1, l1> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim2_0, int Dim3_0, int Dim0_1, int Dim1_1, int Dim2_1,
            int Dim3_1, char i0, char j0, char k0, char l0, char i1, char j1,
            char k1, char l1>
  auto operator+(
    const Tensor4_Expr<A, T, Dim0_0, Dim1_0, Dim2_0, Dim3_0, i0, j0, k0, l0> &a,
    const Tensor4_Expr<B, U, Dim0_1, Dim1_1, Dim2_1, Dim3_1, i1, j1, k1, l1> &b)
  {
    using TensorExpr
      = Tensor4_plus_Tensor4<A, B, T, U, Dim0_0, Dim1_0, Dim2_0, Dim3_0,
                             Dim0_1, Dim1_1, Dim2_1, Dim3_1, i0, j0, k0, l0,
                             i1, j1, k1, l1>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0_0, Dim1_0,
                        Dim2_0, Dim3_0, i0, j0, k0, l0>(TensorExpr(a, b));
  }
}
