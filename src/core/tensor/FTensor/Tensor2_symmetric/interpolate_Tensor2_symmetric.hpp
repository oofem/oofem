/* Interpolates to (i0+distance[0], j0+distance[1], k0+distance[2])
   yielding a Tensor2_symmetric.  (i0,j0,k0) are 3D array coordinates,
   conjugate==1-distance, and (di,dj,dk) are the stride of the array
   coordinates.  It is assumed that the Tensor2_symmetric<T*,Dim> is zero
   centered. */

#pragma once

namespace FTensor
{
  template <class T, int Dim, char i, char j>
  class interpolate_Tensor2_symmetric
  {
    const Tensor2_symmetric<T *, Dim> &a;
    const int di, dj, dk, i0, j0, k0;
    const double *distance, *conjugate;

  public:
    typename promote<T, double>::V operator()(const int N1, const int N2) const
    {
      return conjugate[0] * conjugate[1] * conjugate[2]
               * (*(a.ptr(N1, N2) + di * i0 + dj * j0 + dk * k0))
             + distance[0] * conjugate[1] * conjugate[2]
                 * (*(a.ptr(N1, N2) + di * (i0 + 1) + dj * (j0) + dk * (k0)))
             + conjugate[0] * distance[1] * conjugate[2]
                 * (*(a.ptr(N1, N2) + di * (i0) + dj * (j0 + 1) + dk * (k0)))
             + distance[0] * distance[1] * conjugate[2]
                 * (*(a.ptr(N1, N2) + di * (i0 + 1) + dj * (j0 + 1)
                      + dk * (k0)))
             + conjugate[0] * conjugate[1] * distance[2]
                 * (*(a.ptr(N1, N2) + di * (i0) + dj * (j0) + dk * (k0 + 1)))
             + distance[0] * conjugate[1] * distance[2]
                 * (*(a.ptr(N1, N2) + di * (i0 + 1) + dj * (j0)
                      + dk * (k0 + 1)))
             + conjugate[0] * distance[1] * distance[2]
                 * (*(a.ptr(N1, N2) + di * (i0) + dj * (j0 + 1)
                      + dk * (k0 + 1)))
             + distance[0] * distance[1] * distance[2]
                 * (*(a.ptr(N1, N2) + di * (i0 + 1) + dj * (j0 + 1)
                      + dk * (k0 + 1)));
    }
    interpolate_Tensor2_symmetric(const Tensor2_symmetric<T *, Dim> &A,
                                  const int Di, const int Dj, const int Dk,
                                  const int I0, const int J0, const int K0,
                                  const double Distance[3],
                                  const double Conjugate[3])
        : a(A), di(Di), dj(Dj), dk(Dk), i0(I0), j0(J0), k0(K0),
          distance(Distance), conjugate(Conjugate)
    {}
  };

  template <class T, int Dim, char i, char j>
  const Tensor2_symmetric_Expr<const interpolate_Tensor2_symmetric<T, Dim, i, j>,
                               typename promote<T, double>::V, Dim, i, j>
  interpolate(const Tensor2_symmetric<T *, Dim> &a, const Index<i, Dim> index1,
              const Index<j, Dim> index2, const int &di, const int &dj,
              const int &dk, const int &i0, const int &j0, const int &k0,
              const double distance[3], const double conjugate[3])
  {
    using Tensor_Expr = interpolate_Tensor2_symmetric<T, Dim, i, j>;
    return Tensor2_symmetric_Expr<Tensor_Expr, typename promote<T, double>::V,
                                  Dim, i, j>(
      Tensor_Expr(a, di, dj, dk, i0, j0, k0, distance, conjugate));
  }
}
