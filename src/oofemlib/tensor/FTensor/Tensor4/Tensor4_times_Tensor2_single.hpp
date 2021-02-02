#pragma once

namespace FTensor
{
  // A(i,j,k,l) * B(m,n) single contraction
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, int Dim5, char i, char j, char k, char l,
            char m, char n, int DimA, int DimB, int DimC, int DimX, int DimY,
            char a, char b, char c, char x, char y>
  class Tensor4_times_Tensor2_single
  {
    Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
    Tensor2_Expr<B, U, Dim4, Dim5, m, n> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      typename promote<T, U>::V result(0);
      for(int xx = 0; xx < DimX; ++xx)
        {
          // Permutation is where the indices get checked.
          result += Permutation4<DimA, DimB, DimC, DimX, a, b, c, x>().eval(
                      iterA, N1, N2, N3, xx)
                    * Permutation2<DimX, DimY, x, y>().eval(iterB, xx, N4);
        }
      return result;
    }

    Tensor4_times_Tensor2_single(
      const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &iter_a,
      const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &iter_b)
        : iterA(iter_a), iterB(iter_b)
    {}
  };

  // A(i,j,k,l)*B(l,m)

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3,
                                     Dim4, i, j, k, l, l, m, Dim0, Dim1, Dim2,
                                     Dim3, Dim4, i, j, k, l, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
  }

  // B(l,m)*A(i,j,k,l)

  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(m,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                     Dim3, i, j, k, l, m, l, Dim0, Dim1, Dim2,
                                     Dim3, Dim4, i, j, k, l, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim2, Dim4, i, j, k, m>(TensorExpr(a, b));
  }

  // B(m,l)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(k,m)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2,
                                     Dim4, i, j, k, l, k, m, Dim0, Dim1, Dim3,
                                     Dim2, Dim4, i, j, l, k, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
  }

  // B(k,m)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(m,k)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                     Dim2, i, j, k, l, m, k, Dim0, Dim1, Dim3,
                                     Dim2, Dim4, i, j, l, k, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1,
                        Dim3, Dim4, i, j, l, m>(TensorExpr(a, b));
  }

  // B(m,k)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(j,m)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1,
                                     Dim4, i, j, k, l, j, m, Dim0, Dim2, Dim3,
                                     Dim1, Dim4, i, k, l, j, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2,
                        Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
  }

  // B(j,m)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(m,j)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                     Dim1, i, j, k, l, m, j, Dim0, Dim2, Dim3,
                                     Dim1, Dim4, i, k, l, j, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2,
                        Dim3, Dim4, i, k, l, m>(TensorExpr(a, b));
  }

  // B(m,j)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(i,m)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0,
                                     Dim4, i, j, k, l, i, m, Dim1, Dim2, Dim3,
                                     Dim0, Dim4, j, k, l, i, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2,
                        Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
  }

  // B(i,m)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }

  // A(i,j,k,l)*B(m,i)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
            const Tensor2_Expr<B, U, Dim4, Dim0, m, i> &b)
  {
    using TensorExpr
      = Tensor4_times_Tensor2_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4,
                                     Dim0, i, j, k, l, m, i, Dim1, Dim2, Dim3,
                                     Dim0, Dim4, j, k, l, i, m>;
    return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2,
                        Dim3, Dim4, j, k, l, m>(TensorExpr(a, b));
  }

  // B(m,i)*A(i,j,k,l)
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, char i, char j, char k, char l, char m>
  auto
  operator*(const Tensor2_Expr<B, U, Dim4, Dim0, m, i> &b,
            const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
  {
    return a * b;
  }
}
