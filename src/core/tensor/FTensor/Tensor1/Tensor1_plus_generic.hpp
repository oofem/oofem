/* Adds a Tensor1 to a generic (or vice versa), yielding a Tensor1.
   Usually used for doubles, but could be used for complex, etc.  All
   that it requires is that you can add an element of the Tensor1 to
   it.  */

/* A(i) + d0 -> Tensor1 */

#pragma once

namespace FTensor
{
  template <class A, class T, class U, int Dim, char i>
  class Tensor1_plus_generic
  {
    Tensor1_Expr<A, T, Dim, i> iterA;
    U d;

  public:
    typename promote<T, U>::V operator()(const int N) const
    {
      return iterA(N) + d;
    }

    Tensor1_plus_generic(const Tensor1_Expr<A, T, Dim, i> &a, const U &d0)
        : iterA(a), d(d0)
    {}
  };

  template <class A, class T, class U, int Dim, char i>
  Tensor1_Expr<Tensor1_plus_generic<A, T, U, Dim, i>,
               typename promote<T, U>::V, Dim, i>
  operator+(const Tensor1_Expr<A, T, Dim, i> &a, const U &d0)
  {
    using TensorExpr = Tensor1_plus_generic<A, T, U, Dim, i>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, d0));
  }

  /* d0 + A(i) -> Tensor1 */

  template <class A, class T, class U, int Dim, char i>
  Tensor1_Expr<Tensor1_plus_generic<A, T, U, Dim, i>,
               typename promote<T, U>::V, Dim, i>
  operator+(const U &d0, const Tensor1_Expr<A, T, Dim, i> &a)
  {
    using TensorExpr = Tensor1_plus_generic<A, T, U, Dim, i>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, d0));
  }
}
