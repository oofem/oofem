/* Various assignment operators.  I have to explicitly declare the
   second operator= because otherwise the compiler will generate its
   own and not use the template code. */

#pragma once

namespace FTensor
{
  /* T3as=T3as */

  template <class A, class B, class U, int Dim0, int Dim12, char i, char j,
            char k, int Current_Dim0, int Current_Dim1, int Current_Dim2>
  void T3as_equals_T3as(
    A &iter,
    const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, j, k> &result,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
    const Number<Current_Dim2> &)
  {
    iter.unsafe(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
      = result(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1);
    T3as_equals_T3as(iter, result, Number<Current_Dim0>(),
                     Number<Current_Dim1 - 1>(), Number<Current_Dim2>());
  }

  template <class A, class B, class U, int Dim0, int Dim12, char i, char j,
            char k, int Current_Dim0, int Current_Dim2>
  void T3as_equals_T3as(
    A &iter,
    const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, j, k> &result,
    const Number<Current_Dim0> &, const Number<1> &,
    const Number<Current_Dim2> &)
  {
    iter.unsafe(Current_Dim0 - 1, 0, Current_Dim2 - 1)
      = result(Current_Dim0 - 1, 0, Current_Dim2 - 1);
    T3as_equals_T3as(iter, result, Number<Current_Dim0>(),
                     Number<Current_Dim2 - 2>(), Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class U, int Dim0, int Dim12, char i, char j,
            char k, int Current_Dim0>
  void T3as_equals_T3as(
    A &iter,
    const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, j, k> &result,
    const Number<Current_Dim0> &, const Number<1> &, const Number<2> &)
  {
    iter.unsafe(Current_Dim0 - 1, 0, 1) = result(Current_Dim0 - 1, 0, 1);
    T3as_equals_T3as(iter, result, Number<Current_Dim0 - 1>(),
                     Number<Dim12 - 1>(), Number<Dim12>());
  }

  template <class A, class B, class U, int Dim0, int Dim12, char i, char j,
            char k>
  void T3as_equals_T3as(
    A &iter,
    const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, j, k> &result,
    const Number<1> &, const Number<1> &, const Number<2> &)
  {
    iter.unsafe(0, 0, 1) = result(0, 0, 1);
  }

  template <class A, class T, int Dim0, int Dim12, char i, char j, char k>
  template <class B, class U>
  Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T, Dim0,
                             Dim12, i, j, k> &
  Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T, Dim0,
                             Dim12, i, j, k>::
  operator=(
    const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, j, k> &result)
  {
    T3as_equals_T3as(iter, result, Number<Dim0>(), Number<Dim12 - 1>(),
                     Number<Dim12>());
    return *this;
  }

  /* T3as=T3as_Expr(T3as) */

  template <class A, class T, int Dim0, int Dim12, char i, char j, char k>
  Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T, Dim0,
                             Dim12, i, j, k> &
  Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T, Dim0,
                             Dim12, i, j, k>::
  operator=(
    const Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T,
                                     Dim0, Dim12, i, j, k> &result)
  {
    return operator=<Tensor3_antisymmetric<A, Dim0, Dim12>, T>(result);
  }

  /* This is for when the indices are switched (i,j,k) -> (i,k,j). */

  template <class A, class B, class U, int Dim0, int Dim12, char i, char j,
            char k, int Current_Dim0, int Current_Dim1, int Current_Dim2>
  void T3as_switched_equals_T3as(
    A &iter,
    const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, k, j> &result,
    const Number<Current_Dim0> &, const Number<Current_Dim1> &,
    const Number<Current_Dim2> &)
  {
    iter.unsafe(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1)
      = -result(Current_Dim0 - 1, Current_Dim1 - 1, Current_Dim2 - 1);
    T3as_switched_equals_T3as(iter, result, Number<Current_Dim0>(),
                              Number<Current_Dim1 - 1>(),
                              Number<Current_Dim2>());
  }

  template <class A, class B, class U, int Dim0, int Dim12, char i, char j,
            char k, int Current_Dim0, int Current_Dim2>
  void T3as_switched_equals_T3as(
    A &iter,
    const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, k, j> &result,
    const Number<Current_Dim0> &, const Number<1> &,
    const Number<Current_Dim2> &)
  {
    iter.unsafe(Current_Dim0 - 1, 0, Current_Dim2 - 1)
      = -result(Current_Dim0 - 1, 0, Current_Dim2 - 1);
    T3as_switched_equals_T3as(iter, result, Number<Current_Dim0>(),
                              Number<Current_Dim2 - 2>(),
                              Number<Current_Dim2 - 1>());
  }

  template <class A, class B, class U, int Dim0, int Dim12, char i, char j,
            char k, int Current_Dim0>
  void T3as_switched_equals_T3as(
    A &iter,
    const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, k, j> &result,
    const Number<Current_Dim0> &, const Number<1> &, const Number<2> &)
  {
    iter.unsafe(Current_Dim0 - 1, 0, 1) = -result(Current_Dim0 - 1, 0, 1);
    T3as_switched_equals_T3as(iter, result, Number<Current_Dim0 - 1>(),
                              Number<Dim12 - 1>(), Number<Dim12>());
  }

  template <class A, class B, class U, int Dim0, int Dim12, char i, char j,
            char k>
  void T3as_switched_equals_T3as(
    A &iter,
    const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, k, j> &result,
    const Number<1> &, const Number<1> &, const Number<2> &)
  {
    iter.unsafe(0, 0, 1) = -result(0, 0, 1);
  }

  template <class A, class T, int Dim0, int Dim12, char i, char j, char k>
  template <class B, class U>
  Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T, Dim0,
                             Dim12, i, j, k> &
  Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T, Dim0,
                             Dim12, i, j, k>::
  operator=(
    const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, k, j> &result)
  {
    T3as_switched_equals_T3as(iter, result, Number<Dim0>(),
                              Number<Dim12 - 1>(), Number<Dim12>());
    return *this;
  }
}
