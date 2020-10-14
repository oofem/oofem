/* Takes a second derivative of a Tensor2_symmetric, yielding a
   Ddg.  This is mostly useful near boundaries where you might
   have to take one-sided derivatives. */

#pragma once

namespace FTensor
{
  template <class T, int Dim01, int Dim23, char i, char j, char k, char l>
  class dd_boundary_Tensor2_symmetric
  {
    const Tensor2_symmetric<T *, Dim01> &a;
    const Tensor1<int, Dim23> &d_ijk;
    const Tensor1<double, Dim23> &d_xyz;
    const Tensor2<bool, Dim23, 2> &boundary;

  public:
    typename promote<T, double>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return N3 == N4
               ? (boundary(N3, 0)
                    ? (*(a.ptr(N1, N2) + 2 * d_ijk(N3))
                       - 2 * *(a.ptr(N1, N2) + d_ijk(N3)) + a(N1, N2))
                    : (boundary(N3, 1)
                         ? (a(N1, N2) - 2 * *(a.ptr(N1, N2) - d_ijk(N3))
                            + *(a.ptr(N1, N2) - 2 * d_ijk(N3)))
                         : (*(a.ptr(N1, N2) + d_ijk(N3)) - 2 * a(N1, N2)
                            + *(a.ptr(N1, N2) - d_ijk(N3)))))
                   * d_xyz(N3) * d_xyz(N3)
               : (boundary(N3, 0)
                    ? (boundary(N4, 0)
                         ? ((*(a.ptr(N1, N2) + d_ijk(N3) + d_ijk(N4))
                             - *(a.ptr(N1, N2) + d_ijk(N4))
                             - *(a.ptr(N1, N2) + d_ijk(N3)) + a(N1, N2))
                            * d_xyz(N3) * d_xyz(N4))
                         : (boundary(N4, 1)
                              ? ((*(a.ptr(N1, N2) + d_ijk(N3)) - a(N1, N2)
                                  - *(a.ptr(N1, N2) + d_ijk(N3) - d_ijk(N4))
                                  + *(a.ptr(N1, N2) - d_ijk(N4)))
                                 * d_xyz(N3) * d_xyz(N4))
                              : (*(a.ptr(N1, N2) + d_ijk(N3) + d_ijk(N4))
                                 - *(a.ptr(N1, N2) + d_ijk(N4))
                                 - *(a.ptr(N1, N2) + d_ijk(N3) - d_ijk(N4))
                                 + *(a.ptr(N1, N2) - d_ijk(N4)))
                                  * d_xyz(N3) * d_xyz(N4) * 0.5))
                    : (boundary(N3, 1)
                         ? (boundary(N4, 0)
                              ? ((*(a.ptr(N1, N2) + d_ijk(N4))
                                  - *(a.ptr(N1, N2) - d_ijk(N3) + d_ijk(N4))
                                  - a(N1, N2) + *(a.ptr(N1, N2) - d_ijk(N3)))
                                 * d_xyz(N3) * d_xyz(N4))
                              : (boundary(N4, 1)
                                   ? ((a(N1, N2) - *(a.ptr(N1, N2) - d_ijk(N3))
                                       - *(a.ptr(N1, N2) - d_ijk(N4))
                                       + *(a.ptr(N1, N2) - d_ijk(N3)
                                           - d_ijk(N4)))
                                      * d_xyz(N3) * d_xyz(N4))
                                   : (*(a.ptr(N1, N2) + d_ijk(N4))
                                      - *(a.ptr(N1, N2) - d_ijk(N3) + d_ijk(N4))
                                      - *(a.ptr(N1, N2) - d_ijk(N4))
                                      + *(a.ptr(N1, N2) - d_ijk(N3)
                                          - d_ijk(N4)))
                                       * d_xyz(N3) * d_xyz(N4) * 0.5))
                         : (boundary(N4, 0)
                              ? ((*(a.ptr(N1, N2) + d_ijk(N3) + d_ijk(N4))
                                  - *(a.ptr(N1, N2) - d_ijk(N3) + d_ijk(N4))
                                  - *(a.ptr(N1, N2) + d_ijk(N3))
                                  + *(a.ptr(N1, N2) - d_ijk(N3)))
                                 * d_xyz(N3) * d_xyz(N4) * 0.5)
                              : (boundary(N4, 1)
                                   ? ((*(a.ptr(N1, N2) + d_ijk(N3))
                                       - *(a.ptr(N1, N2) - d_ijk(N3))
                                       - *(a.ptr(N1, N2) + d_ijk(N3)
                                           - d_ijk(N4))
                                       + *(a.ptr(N1, N2) - d_ijk(N3)
                                           - d_ijk(N4)))
                                      * d_xyz(N3) * d_xyz(N4) * 0.5)
                                   : ((*(a.ptr(N1, N2) + d_ijk(N3) + d_ijk(N4))
                                       - *(a.ptr(N1, N2) - d_ijk(N3)
                                           + d_ijk(N4))
                                       - *(a.ptr(N1, N2) + d_ijk(N3)
                                           - d_ijk(N4))
                                       + *(a.ptr(N1, N2) - d_ijk(N3)
                                           - d_ijk(N4)))
                                      * d_xyz(N3) * d_xyz(N4) * 0.25)))));
    }
    dd_boundary_Tensor2_symmetric(const Tensor2_symmetric<T *, Dim01> &A,
                                  const Tensor1<int, Dim23> &D_ijk,
                                  const Tensor1<double, Dim23> &D_xyz,
                                  const Tensor2<bool, Dim23, 2> &Boundary)
        : a(A), d_ijk(D_ijk), d_xyz(D_xyz), boundary(Boundary)
    {}
  };
  /*
  template <class T, int Dim01, int Dim23, char i, char j, char k, char l>
  const Ddg_Expr<
    const dd_boundary_Tensor2_symmetric<T, Dim01, Dim23, i, j, k, l>,
    typename promote<T, double>::V, Dim01, Dim23, i, j, k, l>
  dd_boundary(const Tensor2_symmetric<T *, Dim01> &a,
              const Index<i, Dim01> index1, const Index<j, Dim01> index2,
              const Index<k, Dim23> index3, const Index<l, Dim23> index4,
              const Tensor1<int, Dim23> &d_ijk,
              const Tensor1<double, Dim23> &d_xyz,
              const Tensor2<bool, Dim23, 2> &boundary)
  {
    using Tensor_Expr
      = dd_boundary_Tensor2_symmetric<T, Dim01, Dim23, i, j, k, l>;
    return Ddg_Expr<Tensor_Expr, typename promote<T, double>::V, Dim01, Dim23,
                    i, j, k, l>(Tensor_Expr(a, d_ijk, d_xyz, boundary));
		    }*/
}
