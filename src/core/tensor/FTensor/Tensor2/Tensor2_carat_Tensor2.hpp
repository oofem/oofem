/* Creates a Tensor2_symmetric expression by contracting two Tensor2's
   together. There are different versions, depending on where the
   contracting indices are located (i.e. whether it is A(i,j)^B(j,k)
   or A(i,j)^B(k,j)).  The classes are numbered to differentiate
   between these.  Thus, A(i,j)^B(j,k) has 10 appended to the name
   because I count from 0. */

#pragma once

namespace FTensor
{
  /* Base Template */
  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim0_1, int Dim1_1, char i0, char j0, char i1, char j1>
  class Tensor2_carat_Tensor2
  {};

  /* A(i,j)^B(j,k) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  class Tensor2_carat_Tensor2<A, B, T, U, Dim, Dim1, Dim1, Dim, i, j, j, k>
  {
    const Tensor2_Expr<A, T, Dim, Dim1, i, j> iterA;
    const Tensor2_Expr<B, U, Dim1, Dim, j, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1) * iterB(Current_Dim - 1, N2)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(N1, 0) * iterB(0, N2);
    }

  public:
    Tensor2_carat_Tensor2(const Tensor2_Expr<A, T, Dim, Dim1, i, j> &a,
                          const Tensor2_Expr<B, U, Dim1, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim1>());
    }
  };

  /* A(i,j)^B(k,j) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  class Tensor2_carat_Tensor2<A, B, T, U, Dim, Dim1, Dim, Dim1, i, j, k, j>
  {
    const Tensor2_Expr<A, T, Dim, Dim1, i, j> iterA;
    const Tensor2_Expr<B, U, Dim, Dim1, k, j> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1) * iterB(N2, Current_Dim - 1)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(N1, 0) * iterB(N2, 0);
    }

  public:
    Tensor2_carat_Tensor2(const Tensor2_Expr<A, T, Dim, Dim1, i, j> &a,
                          const Tensor2_Expr<B, U, Dim, Dim1, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim1>());
    }
  };

  /* A(j,i)^B(j,k) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  class Tensor2_carat_Tensor2<A, B, T, U, Dim1, Dim, Dim1, Dim, j, i, j, k>
  {
    const Tensor2_Expr<A, T, Dim1, Dim, j, i> iterA;
    const Tensor2_Expr<B, U, Dim1, Dim, j, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(Current_Dim - 1, N1) * iterB(Current_Dim - 1, N2)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(0, N1) * iterB(0, N2);
    }

  public:
    Tensor2_carat_Tensor2(const Tensor2_Expr<A, T, Dim1, Dim, j, i> &a,
                          const Tensor2_Expr<B, U, Dim1, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim1>());
    }
  };

  /* A(j,i)^B(k,j) */

  template <class A, class B, class T, class U, int Dim, int Dim1, char i,
            char j, char k>
  class Tensor2_carat_Tensor2<A, B, T, U, Dim1, Dim, Dim, Dim1, j, i, k, j>
  {
    const Tensor2_Expr<A, T, Dim1, Dim, j, i> iterA;
    const Tensor2_Expr<B, U, Dim, Dim1, k, j> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<Current_Dim> &) const
    {
      return iterA(Current_Dim - 1, N1) * iterB(N2, Current_Dim - 1)
             + eval(N1, N2, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const Number<1> &) const
    {
      return iterA(0, N1) * iterB(N2, 0);
    }

  public:
    Tensor2_carat_Tensor2(const Tensor2_Expr<A, T, Dim1, Dim, j, i> &a,
                          const Tensor2_Expr<B, U, Dim, Dim1, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim1>());
    }
  };

  template <class A, class B, class T, class U, int Dim0_0, int Dim1_0,
            int Dim0_1, int Dim1_1, char i0, char j0, char i1, char j1>
  auto operator^(const Tensor2_Expr<A, T, Dim0_0, Dim1_0, i0, j0> &a,
                 const Tensor2_Expr<B, U, Dim0_1, Dim1_1, i1, j1> &b)
  {
    using TensorExpr = Tensor2_carat_Tensor2<A, B, T, U, Dim0_0, Dim1_0,
                                             Dim0_1, Dim1_1, i0, j0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the ^ operator");

    // Definition of Helper constexpr variables
    constexpr int Dim = (i0 == i1 || i0 == j1) ? Dim1_0 : Dim0_0;
    constexpr char i = (i0 == i1 || i0 == j1) ? j0 : i0,
                   j = (i1 == i0 || i1 == j0) ? j1 : i1;

    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim,
                                  i, j>(TensorExpr(a, b));
  }

  /* I don't think that this product actually gives a Ddg. */

  //  /* A(i,j)^B(k,l) -> Ddg(i,k,j,l) */

  //  template<class A, class B, class T, class U, int Dim, int Dim1,
  //    char i, char j, char k>
  //  class Tensor2_carat_Tensor2_0213
  //  {
  //    const Tensor2_Expr<A,T,Dim01,Dim23,i,j> iterA;
  //    const Tensor2_Expr<B,U,Dim01,Dim23,k,l> iterB;
  //  public:
  //    Tensor2_carat_Tensor2_0213(const Tensor2_Expr<A,T,Dim01,Dim23,i,j> &a,
  //  			     const Tensor2_Expr<B,U,Dim01,Dim23,k,l> &b):
  //      iterA(a), iterB(b) {}
  //    typename promote<T,U>::V operator()(const int N1, const int N2, const
  //    int N3,
  //  			     const int N4) const
  //    {
  //      return iterA(N1,N3)*iterB(N2,N4);
  //    }
  //  };

  //  template<class A, class B, class T, class U, int Dim01, int Dim23,
  //    char i, char j, char k, char l>
  //  const Ddg_Expr<const
  //  Tensor2_carat_Tensor2_0213<A,B,T,U,Dim01,Dim23,i,j,k,l>,typename
  //  promote<T,U>::V,Dim01,Dim23,i,k,j,l> operator^(const
  //  Tensor2_Expr<A,T,Dim01,Dim23,i,j> &a, 	  const
  //  Tensor2_Expr<B,U,Dim01,Dim23,k,l> &b)
  //  {
  //    typedef Tensor2_carat_Tensor2_0213<A,B,T,U,Dim01,Dim23,i,j,k,l>
  //      TensorExpr;
  //    return Ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,k,j,l>
  //      (TensorExpr(a,b));
  //  }
}
