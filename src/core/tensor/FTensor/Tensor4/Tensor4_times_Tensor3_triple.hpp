#pragma once

namespace FTensor
{
  // A(i,j,k,l) * B(m,n,o) triple contraction
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k,
            char l, char m, char n, char o, int DimA, int DimX, int DimY,
            int DimZ, char a, char x, char y, char z>
  class Tensor4_times_Tensor3_triple
  {
    Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1) const
    {
      typename promote<T, U>::V result(0);
      for(int xx = 0; xx < DimX; ++xx)
        for(int yy = 0; yy < DimY; ++yy)
          for(int zz = 0; zz < DimZ; ++zz)
            {
              // Permutation is where the indices get checked.
              result
                += Permutation4<DimA, DimX, DimY, DimZ, a, x, y, z>().eval(
                     iterA, N1, xx, yy, zz)
                   * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx,
                                                                    yy, zz);
            }
      return result;
    }

    Tensor4_times_Tensor3_triple(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &iter_a,
      const Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> &iter_b)
        : iterA(iter_a), iterB(iter_b)
    {}
  };

  // A(i,j,k,l)*B(j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim1, Dim2, Dim3, j, k, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
                                     Dim2, Dim3, i, j, k, l, j, k, l, Dim0,
                                     Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  // B(j,k,l)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim3, j, k, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(j,l,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim1, Dim3, Dim2, j, l, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
                                     Dim3, Dim2, i, j, k, l, j, l, k, Dim0,
                                     Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  // B(j,l,k)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim2, j, l, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(k,j,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim2, Dim1, Dim3, k, j, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
                                     Dim1, Dim3, i, j, k, l, k, j, l, Dim0,
                                     Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  // B(k,j,l)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim3, k, j, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(k,l,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim2, Dim3, Dim1, k, l, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
                                     Dim3, Dim1, i, j, k, l, k, l, j, Dim0,
                                     Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  // B(k,l,j)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim1, k, l, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(l,j,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim3, Dim1, Dim2, l, j, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
                                     Dim1, Dim2, i, j, k, l, l, j, k, Dim0,
                                     Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  // B(l,j,k)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim2, l, j, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(l,k,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim3, Dim2, Dim1, l, k, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
                                     Dim2, Dim1, i, j, k, l, l, k, j, Dim0,
                                     Dim1, Dim2, Dim3, i, j, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim0, i>(
      TensorExpr(a, b));
  }

  // B(l,k,j)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim1, l, k, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(i,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim0, Dim2, Dim3, i, k, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
                                     Dim2, Dim3, i, j, k, l, i, k, l, Dim1,
                                     Dim0, Dim2, Dim3, j, i, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, j>(
      TensorExpr(a, b));
  }

  // B(i,k,l)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim3, i, k, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(i,l,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim0, Dim3, Dim2, i, l, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
                                     Dim3, Dim2, i, j, k, l, i, l, k, Dim1,
                                     Dim0, Dim2, Dim3, j, i, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, j>(
      TensorExpr(a, b));
  }

  // B(i,l,k)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim2, i, l, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(k,i,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim2, Dim0, Dim3, k, i, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
                                     Dim0, Dim3, i, j, k, l, k, i, l, Dim1,
                                     Dim0, Dim2, Dim3, j, i, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, j>(
      TensorExpr(a, b));
  }

  // B(k,i,l)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim3, k, i, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(k,l,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim2, Dim3, Dim0, k, l, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
                                     Dim3, Dim0, i, j, k, l, k, l, i, Dim1,
                                     Dim0, Dim2, Dim3, j, i, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, j>(
      TensorExpr(a, b));
  }

  // B(k,l,i)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim2, Dim3, Dim0, k, l, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(l,i,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim3, Dim0, Dim2, l, i, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
                                     Dim0, Dim2, i, j, k, l, l, i, k, Dim1,
                                     Dim0, Dim2, Dim3, j, i, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, j>(
      TensorExpr(a, b));
  }

  // B(l,i,k)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim2, l, i, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(l,k,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim3, Dim2, Dim0, l, k, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
                                     Dim2, Dim0, i, j, k, l, l, k, i, Dim1,
                                     Dim0, Dim2, Dim3, j, i, k, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim1, j>(
      TensorExpr(a, b));
  }

  // B(l,k,i)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim3, Dim2, Dim0, l, k, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(i,j,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim0, Dim1, Dim3, i, j, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
                                     Dim1, Dim3, i, j, k, l, i, j, l, Dim2,
                                     Dim0, Dim1, Dim3, k, i, j, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, k>(
      TensorExpr(a, b));
  }

  // B(i,j,l)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim3, i, j, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(i,l,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim0, Dim3, Dim1, i, l, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
                                     Dim3, Dim1, i, j, k, l, i, l, j, Dim2,
                                     Dim0, Dim1, Dim3, k, i, j, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, k>(
      TensorExpr(a, b));
  }

  // B(i,l,j)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim0, Dim3, Dim1, i, l, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(j,i,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim1, Dim0, Dim3, j, i, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
                                     Dim0, Dim3, i, j, k, l, j, i, l, Dim2,
                                     Dim0, Dim1, Dim3, k, i, j, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, k>(
      TensorExpr(a, b));
  }

  // B(j,i,l)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim3, j, i, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(j,l,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim1, Dim3, Dim0, j, l, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
                                     Dim3, Dim0, i, j, k, l, j, l, i, Dim2,
                                     Dim0, Dim1, Dim3, k, i, j, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, k>(
      TensorExpr(a, b));
  }

  // B(j,l,i)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim1, Dim3, Dim0, j, l, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(l,i,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim3, Dim0, Dim1, l, i, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
                                     Dim0, Dim1, i, j, k, l, l, i, j, Dim2,
                                     Dim0, Dim1, Dim3, k, i, j, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, k>(
      TensorExpr(a, b));
  }

  // B(l,i,j)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim3, Dim0, Dim1, l, i, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(l,j,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim3, Dim1, Dim0, l, j, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
                                     Dim1, Dim0, i, j, k, l, l, j, i, Dim2,
                                     Dim0, Dim1, Dim3, k, i, j, l>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, k>(
      TensorExpr(a, b));
  }

  // B(l,j,i)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim3, Dim1, Dim0, l, j, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(i,j,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
                                     Dim1, Dim2, i, j, k, l, i, j, k, Dim3,
                                     Dim0, Dim1, Dim2, l, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim3, l>(
      TensorExpr(a, b));
  }

  // B(i,j,k)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(i,k,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
                                     Dim2, Dim1, i, j, k, l, i, k, j, Dim3,
                                     Dim0, Dim1, Dim2, l, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim3, l>(
      TensorExpr(a, b));
  }

  // B(i,k,j)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(j,i,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
                                     Dim0, Dim2, i, j, k, l, j, i, k, Dim3,
                                     Dim0, Dim1, Dim2, l, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim3, l>(
      TensorExpr(a, b));
  }

  // B(j,i,k)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(j,k,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
                                     Dim2, Dim0, i, j, k, l, j, k, i, Dim3,
                                     Dim0, Dim1, Dim2, l, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim3, l>(
      TensorExpr(a, b));
  }

  // B(j,k,i)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(k,i,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
                                     Dim0, Dim1, i, j, k, l, k, i, j, Dim3,
                                     Dim0, Dim1, Dim2, l, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim3, l>(
      TensorExpr(a, b));
  }

  // B(k,i,j)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(k,j,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
                                     Dim1, Dim0, i, j, k, l, k, j, i, Dim3,
                                     Dim0, Dim1, Dim2, l, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim3, l>(
      TensorExpr(a, b));
  }

  // B(k,j,i)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, char i, char j, char k, char l>
  auto
  operator*(const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }
}
