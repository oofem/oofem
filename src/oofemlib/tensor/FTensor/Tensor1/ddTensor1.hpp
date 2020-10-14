/* Takes a second derivative of a Tensor1, yielding a Dg, but
   with it's indices switched around so that it has the symmetries of
   a Christof. */

#pragma once

namespace FTensor
{
  template <class T, int Dim0, int Dim12, char i, char j, char k>
  class ddTensor1
  {
    const Tensor1<T *, Dim0> &a;
    const Tensor1<int, Dim12> &d_ijk;
    const Tensor1<double, Dim12> &d_xyz;

  public:
    /* Note the indices switch here. */

    typename promote<T, double>::V
    operator()(const int N2, const int N3, const int N1) const
    {
      return N2 == N3 ? (*(a.ptr(N1) + d_ijk(N2)) - 2 * a(N1)
                         + *(a.ptr(N1) - d_ijk(N2)))
                          * d_xyz(N2) * d_xyz(N2)
                      : (*(a.ptr(N1) + d_ijk(N2) + d_ijk(N3))
                         - *(a.ptr(N1) - d_ijk(N2) + d_ijk(N3))
                         - *(a.ptr(N1) + d_ijk(N2) - d_ijk(N3))
                         + *(a.ptr(N1) - d_ijk(N2) - d_ijk(N3)))
                          * d_xyz(N2) * d_xyz(N3) * 0.25;
    }
    ddTensor1(const Tensor1<T *, Dim0> &A, const Tensor1<int, Dim12> &D_ijk,
              const Tensor1<double, Dim12> &D_xyz)
        : a(A), d_ijk(D_ijk), d_xyz(D_xyz)
    {}
  };

  /* Note that the indices are mixed up here to convert between being
     symmetric on the last two indices to the first two indices. */

  template <class T, int Dim0, int Dim12, char i, char j, char k>
  const Dg_Expr<const ddTensor1<T, Dim0, Dim12, i, j, k>,
                typename promote<T, double>::V, Dim0, Dim12, i, j, k>
  dd(const Tensor1<T *, Dim0> &a, const Index<k, Dim0> index0,
     const Index<i, Dim12> index1, const Index<j, Dim12> index2,
     const Tensor1<int, Dim12> &d_ijk, const Tensor1<double, Dim12> &d_xyz)
  {
    using TensorExpr = ddTensor1<T, Dim0, Dim12, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, double>::V, Dim0, Dim12, i,
                   j, k>(TensorExpr(a, d_ijk, d_xyz));
  }
}
