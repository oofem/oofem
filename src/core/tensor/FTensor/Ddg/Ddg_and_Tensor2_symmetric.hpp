/* Multiply a Tensor2_symmetric and a Ddg together but don't
   contract, yielding a Ddg. */

#pragma once

namespace FTensor
{
  /* Base Template */
  template <class A, class B, class T, class U, int Dim01_0, int Dim23_0,
            int Dim_1, char i0, char j0, char k0, char l0, char i1, char j1>
  class Ddg_and_Tensor2_symmetric
  {};

  /* A(i,j,k,l) & B(i,j) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_and_Tensor2_symmetric<A, B, T, U, Dim01, Dim23, Dim01, i, j, k, l,
                                  i, j>
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor2_symmetric_Expr<B, U, Dim01, i, j> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) * iterB(N1, N2);
    }

    Ddg_and_Tensor2_symmetric(
      const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim01, i, j> &b)
        : iterA(a), iterB(b)
    {}
  };

  /* B(i,j) & A(i,j,k,l) -> Ddg */

  /* A(i,j,k,l) & B(k,l) -> Ddg */

  template <class A, class B, class T, class U, int Dim01, int Dim23, char i,
            char j, char k, char l>
  class Ddg_and_Tensor2_symmetric<A, B, T, U, Dim01, Dim23, Dim23, i, j, k, l,
                                  k, l>
  {
    Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> iterA;
    Tensor2_symmetric_Expr<B, U, Dim23, k, l> iterB;

  public:
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iterA(N1, N2, N3, N4) * iterB(N3, N4);
    }

    Ddg_and_Tensor2_symmetric(
      const Ddg_Expr<A, T, Dim01, Dim23, i, j, k, l> &a,
      const Tensor2_symmetric_Expr<B, U, Dim23, k, l> &b)
        : iterA(a), iterB(b)
    {}
  };

  template <class A, class B, class T, class U, int Dim01_0, int Dim23_0,
            int Dim_1, char i0, char j0, char k0, char l0, char i1, char j1>
  Ddg_Expr<Ddg_and_Tensor2_symmetric<A, B, T, U, Dim01_0, Dim23_0, Dim_1, i0,
                                     j0, k0, l0, i1, j1>,
           typename promote<T, U>::V, Dim01_0, Dim23_0, i0, j0, k0, l0>
  operator&(const Ddg_Expr<A, T, Dim01_0, Dim23_0, i0, j0, k0, l0> &a,
            const Tensor2_symmetric_Expr<B, U, Dim_1, i1, j1> &b)
  {
    using TensorExpr
      = Ddg_and_Tensor2_symmetric<A, B, T, U, Dim01_0, Dim23_0, Dim_1, i0, j0,
                                  k0, l0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the & operator");
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01_0, Dim23_0,
                    i0, j0, k0, l0>(TensorExpr(a, b));
  }

  /* B(k,l) & A(i,j,k,l) -> Ddg */

  template <class A, class B, class T, class U, int Dim01_0, int Dim23_0,
            int Dim_1, char i0, char j0, char k0, char l0, char i1, char j1>
  Ddg_Expr<Ddg_and_Tensor2_symmetric<A, B, T, U, Dim01_0, Dim23_0, Dim_1, i0,
                                     j0, k0, l0, i1, j1>,
           typename promote<T, U>::V, Dim01_0, Dim23_0, i0, j0, k0, l0>
  operator&(const Tensor2_symmetric_Expr<B, U, Dim_1, i1, j1> &b,
            const Ddg_Expr<A, T, Dim01_0, Dim23_0, i0, j0, k0, l0> &a)
  {
    using TensorExpr
      = Ddg_and_Tensor2_symmetric<A, B, T, U, Dim01_0, Dim23_0, Dim_1, i0, j0,
                                  k0, l0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the & operator");
    return Ddg_Expr<TensorExpr, typename promote<T, U>::V, Dim01_0, Dim23_0,
                    i0, j0, k0, l0>(TensorExpr(a, b));
  }

  /* I originally put these declarations for unknown reasons, but they
     won't work because the result is not a Ddg.  The
     multiplication messes up the symmetries. */

  //  /* A(i,j,k,l) & B(j,l) -> Ddg */

  //  template<class A, class B, class T, class U, int Dim,
  //    char i, char j, char k, char l>
  //  class Ddg_and_Tensor2_symmetric_13
  //  {
  //    const Ddg_Expr<A,T,Dim,Dim,i,j,k,l> iterA;
  //    const Tensor2_symmetric_Expr<B,U,Dim,j,l> iterB;
  //  public:
  //    typename promote<T,U>::V operator()(const int N1, const int N2, const
  //    int N3,
  //  		    const int N4) const
  //    {
  //      return iterA(N1,N2,N3,N4)*iterB(N2,N4);
  //    }

  //    Ddg_and_Tensor2_symmetric_13
  //    (const Ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
  //     const Tensor2_symmetric_Expr<B,U,Dim,j,l> &b): iterA(a), iterB(b) {}
  //  };

  //  template<class A, class B, class T, class U, int Dim,
  //    char i, char j, char k, char l>
  //  const Ddg_Expr
  //  <const Ddg_and_Tensor2_symmetric_13<A,B,T,U,Dim,i,j,k,l>,
  //    typename promote<T,U>::V,Dim,Dim,i,j,k,l>
  //  operator&(const Ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
  //  	  const Tensor2_symmetric_Expr<B,U,Dim,j,l> &b)
  //  {
  //    typedef Ddg_and_Tensor2_symmetric_13<A,B,T,U,Dim,i,j,k,l>
  //      TensorExpr;
  //    return Ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,j,k,l>
  //      (TensorExpr(a,b));
  //  }

  //  /* B(j,l) & A(i,j,k,l) -> Ddg */

  //  template<class A, class B, class T, class U, int Dim,
  //    char i, char j, char k, char l>
  //  const Ddg_Expr
  //  <const Ddg_and_Tensor2_symmetric_13<A,B,T,U,Dim,i,j,k,l>,
  //    typename promote<T,U>::V,Dim,Dim,i,j,k,l>
  //  operator&(const Tensor2_symmetric_Expr<B,U,Dim,j,l> &b,
  //  	  const Ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a)
  //  {
  //    typedef Ddg_and_Tensor2_symmetric_13<A,B,T,U,Dim,i,j,k,l>
  //      TensorExpr;
  //    return Ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,j,k,l>
  //      (TensorExpr(a,b));
  //  }

  //  /* A(i,j,k,l) & B(l,j) -> Ddg */

  //  template<class A, class B, class T, class U, int Dim,
  //    char i, char j, char k, char l>
  //  class Ddg_and_Tensor2_symmetric_31
  //  {
  //    const Ddg_Expr<A,T,Dim,Dim,i,j,k,l> iterA;
  //    const Tensor2_symmetric_Expr<B,U,Dim,l,j> iterB;
  //  public:
  //    typename promote<T,U>::V operator()(const int N1, const int N2, const
  //    int N3,
  //  		    const int N4) const
  //    {
  //      return iterA(N1,N2,N3,N4)*iterB(N2,N4);
  //    }

  //    Ddg_and_Tensor2_symmetric_31
  //    (const Ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
  //     const Tensor2_symmetric_Expr<B,U,Dim,l,j> &b): iterA(a), iterB(b) {}
  //  };

  //  template<class A, class B, class T, class U, int Dim,
  //    char i, char j, char k, char l>
  //  const Ddg_Expr
  //  <const Ddg_and_Tensor2_symmetric_31<A,B,T,U,Dim,i,j,k,l>,
  //    typename promote<T,U>::V,Dim,Dim,i,j,k,l>
  //  operator&(const Ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
  //  	  const Tensor2_symmetric_Expr<B,U,Dim,l,j> &b)
  //  {
  //    typedef Ddg_and_Tensor2_symmetric_31<A,B,T,U,Dim,i,j,k,l>
  //      TensorExpr;
  //    return Ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,j,k,l>
  //      (TensorExpr(a,b));
  //  }

  //  /* B(l,j) & A(i,j,k,l) -> Ddg */

  //  template<class A, class B, class T, class U, int Dim,
  //    char i, char j, char k, char l>
  //  const Ddg_Expr
  //  <const Ddg_and_Tensor2_symmetric_31<A,B,T,U,Dim,i,j,k,l>,
  //    typename promote<T,U>::V,Dim,Dim,i,j,k,l>
  //  operator&(const Tensor2_symmetric_Expr<B,U,Dim,l,j> &b,
  //  	  const Ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a)
  //  {
  //    typedef Ddg_and_Tensor2_symmetric_31<A,B,T,U,Dim,i,j,k,l>
  //      TensorExpr;
  //    return Ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,j,k,l>
  //      (TensorExpr(a,b));
  //  }
}
