/* Declares a wrapper class for rank 4 Tensor expressions symmetric on
   the first two and last two indices. */

#pragma once

#include "Ddg_and_Ddg.hpp"
#include "Ddg_and_Tensor2_symmetric.hpp"
#include "Ddg_carat_Ddg.hpp"
#include "Ddg_carat_Tensor2_symmetric.hpp"
#include "Ddg_minus_Ddg.hpp"
#include "Ddg_mod_Tensor2_symmetric.hpp"
#include "Ddg_or_Ddg.hpp"
#include "Ddg_plus_Ddg.hpp"
#include "Ddg_times_Ddg.hpp"
#include "Ddg_times_Tensor1.hpp"
#include "Ddg_times_Tensor2.hpp"
#include "Ddg_times_Tensor2_symmetric.hpp"
#include "Ddg_times_generic.hpp"
#include "minus_Ddg.hpp"
//  #include "Ddg_mod_Ddg.hpp"

namespace FTensor
{
  template <class A, class T, int Dim01, int Dim23, char i, char j, char k,
            char l>
  class Ddg_Expr
  {
    A iter;

  public:
    Ddg_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iter(N1, N2, N3, N4);
    }
  };

  template <class A, class T, int Tensor_Dim01, int Tensor_Dim23, int Dim01,
            int Dim23, char i, char j, char k, char l>
  class Ddg_Expr<Ddg<A, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i, j, k,
                 l>
  {
    Ddg<A, Tensor_Dim01, Tensor_Dim23> &iter;

  public:
    Ddg_Expr(Ddg<A, Tensor_Dim01, Tensor_Dim23> &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iter(N1, N2, N3, N4);
    }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    Ddg_Expr<Ddg<A, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i, j, k, l> &
    operator=(const Ddg_Expr<B, U, Dim01, Dim23, i, j, k, l> &result);

    Ddg_Expr<Ddg<A, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i, j, k, l> &
    operator=(const Ddg_Expr<Ddg<A, Tensor_Dim01, Tensor_Dim23>, T, Dim01,
                             Dim23, i, j, k, l> &result);
  };
}

#include "Ddg_Expr_equals.hpp"
