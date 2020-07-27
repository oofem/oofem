/* Takes a derivative of a Tensor0<T*> yielding a Tensor1. */

#pragma once

namespace FTensor
{
  template <class T, int Dim, char i> class dTensor0
  {
    const Tensor0<T *> &a;
    const Tensor1<int, Dim> &d_ijk;
    const Tensor1<double, Dim> &d_xyz;

  public:
    typename promote<T, double>::V operator()(const int N) const
    {
      return (*(&a + d_ijk(N)) - *(&a - d_ijk(N))) * d_xyz(N) * 0.5;
    }
    dTensor0(const Tensor0<T *> &A, const Tensor1<int, Dim> &D_ijk,
             const Tensor1<double, Dim> &D_xyz)
        : a(A), d_ijk(D_ijk), d_xyz(D_xyz)
    {}
  };

  template <class T, int Dim, char i>
  const Tensor1_Expr<const dTensor0<T, Dim, i>, typename promote<T, double>::V,
                     Dim, i>
  d(const Tensor0<T *> &a, const Index<i, Dim> index,
    const Tensor1<int, Dim> &d_ijk, const Tensor1<double, Dim> &d_xyz)
  {
    using Tensor_Expr = dTensor0<T, Dim, i>;
    return Tensor1_Expr<Tensor_Expr, typename promote<T, double>::V, Dim, i>(
      Tensor_Expr(a, d_ijk, d_xyz));
  }
}
