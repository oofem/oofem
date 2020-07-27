/* Takes a second derivative of a Tensor0<T*> yielding a Tensor1.
   This is primarily useful at boundaries, where you have to take
   one-sided derivatives. */

#pragma once

namespace FTensor
{
  template <class T, int Dim, char i, char j> class dd_boundary_Tensor0
  {
    const Tensor0<T *> &a;
    const Tensor1<int, Dim> &d_ijk;
    const Tensor1<double, Dim> &d_xyz;
    const Tensor2<bool, Dim, 2> &boundary;

  public:
    typename promote<T, double>::V operator()(const int N1, const int N2) const
    {
      return N1 == N2
               ? (boundary(N1, 0)
                    ? (*(&a + 2 * d_ijk(N1)) - 2 * *(&a + d_ijk(N1)) + a)
                    : (boundary(N1, 1)
                         ? (a - 2 * *(&a - d_ijk(N1)) + *(&a - 2 * d_ijk(N1)))
                         : (*(&a + d_ijk(N1)) - 2 * a + *(&a - d_ijk(N1)))))
                   * d_xyz(N1) * d_xyz(N1)
               : (boundary(N1, 0)
                    ? (boundary(N2, 0)
                         ? ((*(&a + d_ijk(N1) + d_ijk(N2)) - *(&a + d_ijk(N2))
                             - *(&a + d_ijk(N1)) + a)
                            * d_xyz(N1) * d_xyz(N2))
                         : (boundary(N2, 1) ? ((*(&a + d_ijk(N1)) - a
                                                - *(&a + d_ijk(N1) - d_ijk(N2))
                                                + *(&a - d_ijk(N2)))
                                               * d_xyz(N1) * d_xyz(N2))
                                            : (*(&a + d_ijk(N1) + d_ijk(N2))
                                               - *(&a + d_ijk(N2))
                                               - *(&a + d_ijk(N1) - d_ijk(N2))
                                               + *(&a - d_ijk(N2)))
                                                * d_xyz(N1) * d_xyz(N2) * 0.5))
                    : (boundary(N1, 1)
                         ? (boundary(N2, 0)
                              ? ((*(&a + d_ijk(N2))
                                  - *(&a - d_ijk(N1) + d_ijk(N2)) - a
                                  + *(&a - d_ijk(N1)))
                                 * d_xyz(N1) * d_xyz(N2))
                              : (boundary(N2, 1)
                                   ? ((a - *(&a - d_ijk(N1)) - *(&a - d_ijk(N2))
                                       + *(&a - d_ijk(N1) - d_ijk(N2)))
                                      * d_xyz(N1) * d_xyz(N2))
                                   : (*(&a + d_ijk(N2))
                                      - *(&a - d_ijk(N1) + d_ijk(N2))
                                      - *(&a - d_ijk(N2))
                                      + *(&a - d_ijk(N1) - d_ijk(N2)))
                                       * d_xyz(N1) * d_xyz(N2) * 0.5))
                         : (boundary(N2, 0)
                              ? ((*(&a + d_ijk(N1) + d_ijk(N2))
                                  - *(&a - d_ijk(N1) + d_ijk(N2))
                                  - *(&a + d_ijk(N1)) + *(&a - d_ijk(N1)))
                                 * d_xyz(N1) * d_xyz(N2) * 0.5)
                              : (boundary(N2, 1)
                                   ? ((*(&a + d_ijk(N1)) - *(&a - d_ijk(N1))
                                       - *(&a + d_ijk(N1) - d_ijk(N2))
                                       + *(&a - d_ijk(N1) - d_ijk(N2)))
                                      * d_xyz(N1) * d_xyz(N2) * 0.5)
                                   : ((*(&a + d_ijk(N1) + d_ijk(N2))
                                       - *(&a - d_ijk(N1) + d_ijk(N2))
                                       - *(&a + d_ijk(N1) - d_ijk(N2))
                                       + *(&a - d_ijk(N1) - d_ijk(N2)))
                                      * d_xyz(N1) * d_xyz(N2) * 0.25)))));
    }
    dd_boundary_Tensor0(const Tensor0<T *> &A, const Tensor1<int, Dim> &D_ijk,
                        const Tensor1<double, Dim> &D_xyz,
                        const Tensor2<bool, Dim, 2> &Boundary)
        : a(A), d_ijk(D_ijk), d_xyz(D_xyz), boundary(Boundary)
    {}
  };

  template <class T, int Dim, char i, char j>
  const Tensor2_symmetric_Expr<const dd_boundary_Tensor0<T, Dim, i, j>,
                               typename promote<T, double>::V, Dim, i, j>
  dd_boundary(const Tensor0<T *> &a, const Index<i, Dim> index3,
              const Index<j, Dim> index4, const Tensor1<int, Dim> &d_ijk,
              const Tensor1<double, Dim> &d_xyz,
              const Tensor2<bool, Dim, 2> &boundary)
  {
    using Tensor_Expr = dd_boundary_Tensor0<T, Dim, i, j>;
    return Tensor2_symmetric_Expr<Tensor_Expr, typename promote<T, double>::V,
                                  Dim, i, j>(
      Tensor_Expr(a, d_ijk, d_xyz, boundary));
  }
}
