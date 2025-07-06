#pragma once

namespace FTensor
{
  // A(i,j,k,l) * B(m,n) double contraction
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, int Dim5, char i, char j, char k, char l,
            char m, char n, int DimA, int DimB, int DimX, int DimY, char a,
            char b, char x, char y>
  class Tensor4_times_Tensor2_double
  {
    Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    Tensor2_Expr<B, U, Dim4, Dim5, m, n> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      typename promote<T, U>::V result(0);
      for(int xx = 0; xx < DimX; ++xx)
        for(int yy = 0; yy < DimY; ++yy)
          {
            // Permutation is where the indices get checked.
            result += Permutation4<DimA, DimB, DimX, DimY, a, b, x, y>().eval(
                        iterA, N1, N2, xx, yy)
                      * Permutation2<DimX, DimY, x, y>().eval(iterB, xx, yy);
          }
      return result;
    }

    Tensor4_times_Tensor2_double(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &iter_a,
      const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &iter_b)
        : iterA(iter_a), iterB(iter_b)
    {}
  };

  // A(i,j,k,l)*B(k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
                                     Dim3, i, j, k, l, k, l, Dim0, Dim1, Dim2,
                                     Dim3, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, b));
  }

  // B(k,l)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(l,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
                                     Dim2, i, j, k, l, l, k, Dim0, Dim1, Dim2,
                                     Dim3, i, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i,
                        j>(TensorExpr(a, b));
  }

  // B(l,k)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(j,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim1, Dim3, j, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
                                     Dim3, i, j, k, l, j, l, Dim0, Dim2, Dim1,
                                     Dim3, i, k, j, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i,
                        k>(TensorExpr(a, b));
  }

  // B(j,l)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim1, Dim3, j, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(l,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim3, Dim1, l, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
                                     Dim1, i, j, k, l, l, j, Dim0, Dim2, Dim1,
                                     Dim3, i, k, j, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i,
                        k>(TensorExpr(a, b));
  }

  // B(l,j)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim3, Dim1, l, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(j,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
                                     Dim2, i, j, k, l, j, k, Dim0, Dim3, Dim1,
                                     Dim2, i, l, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, i,
                        l>(TensorExpr(a, b));
  }

  // B(j,k)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(k,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
                                     Dim1, i, j, k, l, k, j, Dim0, Dim3, Dim1,
                                     Dim2, i, l, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, i,
                        l>(TensorExpr(a, b));
  }

  // A(i,j,k,l)*B(k,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(i,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
                                     Dim3, i, j, k, l, i, l, Dim1, Dim2, Dim0,
                                     Dim3, j, k, i, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        k>(TensorExpr(a, b));
  }

  // B(i,l)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(l,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
                                     Dim0, i, j, k, l, l, i, Dim1, Dim2, Dim0,
                                     Dim3, j, k, i, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        k>(TensorExpr(a, b));
  }

  // B(l,i)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(i,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim0, Dim2, i, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
                                     Dim2, i, j, k, l, i, k, Dim1, Dim3, Dim0,
                                     Dim2, j, l, i, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        l>(TensorExpr(a, b));
  }

  // B(i,k)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim0, Dim2, i, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(k,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim2, Dim0, k, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
                                     Dim0, i, j, k, l, k, i, Dim1, Dim3, Dim0,
                                     Dim2, j, l, i, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j,
                        l>(TensorExpr(a, b));
  }

  // A(i,j,k,l)*B(k,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim2, Dim0, k, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(i,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
                                     Dim1, i, j, k, l, i, j, Dim2, Dim3, Dim0,
                                     Dim1, k, l, i, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, k,
                        l>(TensorExpr(a, b));
  }

  // A(i,j,k,l)*B(i,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(j,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
                                     Dim0, i, j, k, l, j, i, Dim2, Dim3, Dim0,
                                     Dim1, k, l, i, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, k,
                        l>(TensorExpr(a, b));
  }

  // A(i,j,k,l)*B(j,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }
}
