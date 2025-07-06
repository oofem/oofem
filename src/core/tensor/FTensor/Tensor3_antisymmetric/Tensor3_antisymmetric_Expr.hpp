/* Declare a wrapper class for rank 3 Tensor expressions,
   antisymmetric in the last two indices.  I specialize it for when I
   wrap a simple Tensor3_antisymmetric(_ptr) so that it has a
   reference to the Tensor3_antisymmetric(_ptr) and not a copy.
   Otherwise assignment wouldn't work. */

#pragma once

#include "Tensor3_antisymmetric_or_Tensor3_antisymmetric.hpp"
#include "Tensor3_antisymmetric_plus_Tensor3_antisymmetric.hpp"
#include "Tensor3_antisymmetric_times_Tensor3.hpp"
//  #include "Tensor3_antisymmetric_mod_Tensor1.hpp"
#include "Tensor3_antisymmetric_times_generic.hpp"

namespace FTensor
{
  template <class A, class T, int Dim0, int Dim12, char i, char j, char k>
  class Tensor3_antisymmetric_Expr
  {
    A iter;

  public:
    Tensor3_antisymmetric_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3) const
    {
      return iter(N1, N2, N3);
    }
  };

  template <class A, class T, int Dim0, int Dim12, char i, char j, char k>
  class Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T,
                                   Dim0, Dim12, i, j, k>
  {
    Tensor3_antisymmetric<A, Dim0, Dim12> &iter;

  public:
    Tensor3_antisymmetric_Expr(Tensor3_antisymmetric<A, Dim0, Dim12> &a)
        : iter(a)
    {}
    T operator()(const int N1, const int N2, const int N3) const
    {
      return iter(N1, N2, N3);
    }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T, Dim0,
                               Dim12, i, j, k> &
    operator=(
      const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, j, k> &result);

    Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T, Dim0,
                               Dim12, i, j, k> &
    operator=(
      const Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>,
                                       T, Dim0, Dim12, i, j, k> &result);

    /* This is for when the indices are switched (i,j,k) -> (i,k,j). */

    template <class B, class U>
    Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<A, Dim0, Dim12>, T, Dim0,
                               Dim12, i, j, k> &
    operator=(
      const Tensor3_antisymmetric_Expr<B, U, Dim0, Dim12, i, k, j> &result);
  };
}

#include "Tensor3_antisymmetric_Expr_equals.hpp"
