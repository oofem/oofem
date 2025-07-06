/* Creates a Tensor2_symmetric expression by contracting a
   Tensor2_symmetric and a Tensor2 together. There are different
   versions, depending on where the contracting indices are located
   (i.e. whether it is A(i,j)^B(j,k) or A(i,j)^B(k,j)).  The classes
   are numbered to differentiate between these.  Thus, A(i,j)^B(j,k)
   has 10 appended to the name because I count from 0. */

#pragma once

namespace FTensor
{
  /* Base Template */
  template <class A, class B, class T, class U, int Dim_0, int Dim0_1,
            int Dim1_1, char i0, char j0, char i1, char j1>
  class Tensor2_symmetric_carat_Tensor2
  {};

  /* A(i,j)*B(j,k) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Tensor2_symmetric_carat_Tensor2<A, B, T, U, Dim, Dim, Dim, i, j, j, k>
  {
    const Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    const Tensor2_Expr<B, U, Dim, Dim, j, k> iterB;

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
    Tensor2_symmetric_carat_Tensor2(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
      const Tensor2_Expr<B, U, Dim, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  /* A(i,j)*B(k,j) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Tensor2_symmetric_carat_Tensor2<A, B, T, U, Dim, Dim, Dim, i, j, k, j>
  {
    const Tensor2_symmetric_Expr<A, T, Dim, i, j> iterA;
    const Tensor2_Expr<B, U, Dim, Dim, k, j> iterB;

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
    Tensor2_symmetric_carat_Tensor2(
      const Tensor2_symmetric_Expr<A, T, Dim, i, j> &a,
      const Tensor2_Expr<B, U, Dim, Dim, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  /* A(j,i)*B(j,k) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Tensor2_symmetric_carat_Tensor2<A, B, T, U, Dim, Dim, Dim, j, i, j, k>
  {
    const Tensor2_symmetric_Expr<A, T, Dim, j, i> iterA;
    const Tensor2_Expr<B, U, Dim, Dim, j, k> iterB;

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
    Tensor2_symmetric_carat_Tensor2(
      const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a,
      const Tensor2_Expr<B, U, Dim, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  /* A(j,i)*B(k,j) */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Tensor2_symmetric_carat_Tensor2<A, B, T, U, Dim, Dim, Dim, j, i, k, j>
  {
    const Tensor2_symmetric_Expr<A, T, Dim, j, i> iterA;
    const Tensor2_Expr<B, U, Dim, Dim, k, j> iterB;

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
    Tensor2_symmetric_carat_Tensor2(
      const Tensor2_symmetric_Expr<A, T, Dim, j, i> &a,
      const Tensor2_Expr<B, U, Dim, Dim, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      return eval(N1, N2, Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim_0, int Dim0_1,
            int Dim1_1, char i0, char j0, char i1, char j1>
  auto operator^(const Tensor2_symmetric_Expr<A, T, Dim_0, i0, j0> &a,
                 const Tensor2_Expr<B, U, Dim0_1, Dim1_1, i1, j1> &b)
  {
    using TensorExpr
      = Tensor2_symmetric_carat_Tensor2<A, B, T, U, Dim_0, Dim0_1, Dim1_1, i0,
                                        j0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the ^ operator");

    // Definition of Helper constexpr variables
    constexpr char i = (i0 == i1 || i0 == j1) ? j0 : i0,
                   j = (i1 == i0 || i1 == j0) ? j1 : i1;

    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim_0,
                                  i, j>(TensorExpr(a, b));
  }

  /* B(k,j)*A(j,i) */

  template <class A, class B, class T, class U, int Dim_0, int Dim0_1,
            int Dim1_1, char i0, char j0, char i1, char j1>
  auto operator^(const Tensor2_Expr<B, U, Dim0_1, Dim1_1, i1, j1> &b,
                 const Tensor2_symmetric_Expr<A, T, Dim_0, i0, j0> &a)
  {
    using TensorExpr
      = Tensor2_symmetric_carat_Tensor2<A, B, T, U, Dim_0, Dim0_1, Dim1_1, i0,
                                        j0, i1, j1>;
    static_assert(
      !std::is_empty<TensorExpr>::value,
      "Indexes or Dimensions are not compatible with the ^ operator");

    // Definition of Helper constexpr variables
    constexpr char i = (i0 == i1 || i0 == j1) ? j0 : i0,
                   j = (i1 == i0 || i1 == j0) ? j1 : i1;

    return Tensor2_symmetric_Expr<TensorExpr, typename promote<T, U>::V, Dim_0,
                                  i, j>(TensorExpr(a, b));
  }
}
