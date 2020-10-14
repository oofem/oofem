/* Takes a second derivative of a Tensor0 yielding a Tensor2_symmetric. */

#pragma once

namespace FTensor
{
  template <class T, int Dim, char i, char j> class ddTensor0
  {
    const Tensor0<T *> &a;
    const Tensor1<int, Dim> &d_ijk;
    const Tensor1<double, Dim> &d_xyz;

  public:
    typename promote<T, double>::V operator()(const int N1, const int N2) const
    {
      return N1 == N2
               ? (*(&a + d_ijk(N1)) - 2 * a + *(&a - d_ijk(N1))) * d_xyz(N1)
                   * d_xyz(N1)
               : (*(&a + d_ijk(N1) + d_ijk(N2)) - *(&a - d_ijk(N1) + d_ijk(N2))
                  - *(&a + d_ijk(N1) - d_ijk(N2))
                  + *(&a - d_ijk(N1) - d_ijk(N2)))
                   * d_xyz(N1) * d_xyz(N2) * 0.25;
    }
    ddTensor0(const Tensor0<T *> &A, const Tensor1<int, Dim> &D_ijk,
              const Tensor1<double, Dim> &D_xyz)
        : a(A), d_ijk(D_ijk), d_xyz(D_xyz)
    {}
  };

  template <class T, int Dim, char i, char j>
  const Tensor2_symmetric_Expr<const ddTensor0<T, Dim, i, j>,
                               typename promote<T, double>::V, Dim, i, j>
  dd(const Tensor0<T *> &a, const Index<i, Dim> index1,
     const Index<j, Dim> index2, const Tensor1<int, Dim> &d_ijk,
     const Tensor1<double, Dim> &d_xyz)
  {
    using Tensor_Expr = ddTensor0<T, Dim, i, j>;
    return Tensor2_symmetric_Expr<Tensor_Expr, typename promote<T, double>::V,
                                  Dim, i, j>(Tensor_Expr(a, d_ijk, d_xyz));
  }
}
