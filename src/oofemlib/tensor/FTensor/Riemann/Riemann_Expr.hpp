/* Declares a wrapper class for rank 4 Tensor expressions with Riemann
   symmetries.  */

#pragma once

#include "Riemann_minus_Riemann.hpp"
#include "Riemann_plus_Riemann.hpp"
#include "Riemann_times_Ddg.hpp"
#include "Riemann_times_Tensor1.hpp"
#include "Riemann_times_Tensor2_symmetric.hpp"

namespace FTensor
{
  template <class A, class T, int Dim, char i, char j, char k, char l>
  class Riemann_Expr
  {
    A iter;

  public:
    Riemann_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iter(N1, N2, N3, N4);
    }
  };

  template <class A, class T, int Dim, char i, char j, char k, char l>
  class Riemann_Expr<Riemann<A, Dim>, T, Dim, i, j, k, l>
  {
    Riemann<A, Dim> &iter;

  public:
    Riemann_Expr(Riemann<A, Dim> &a) : iter(a) {}
    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return iter.eval(N1, N2, N3, N4);
    }

    /* Various assignment operators.  I have to explicitly declare the
       second operator= because otherwise the compiler will generate its
       own and not use the template code. */

    template <class B, class U>
    const Riemann_Expr<Riemann<A, Dim>, T, Dim, i, j, k, l> &
    operator=(const Riemann_Expr<B, U, Dim, i, j, k, l> &result)
    {
      iter(0, 1, 0, 1) = result(0, 1, 0, 1);
      iter(0, 1, 0, 2) = result(0, 1, 0, 2);
      iter(0, 2, 0, 2) = result(0, 2, 0, 2);
      iter(0, 1, 1, 2) = result(0, 1, 1, 2);
      iter(0, 2, 1, 2) = result(0, 2, 1, 2);
      iter(1, 2, 1, 2) = result(1, 2, 1, 2);
      return *this;
    }

    const Riemann_Expr<Riemann<A, Dim>, T, Dim, i, j, k, l> &
    operator=(const Riemann_Expr<Riemann<A, Dim>, T, Dim, i, j, k, l> &result)
    {
      return operator=<Riemann<A, Dim>, T>(result);
    }
    template <class B, class U>
    const Riemann_Expr<Riemann<A, Dim>, T, Dim, i, j, k, l> &
    operator+=(const Riemann_Expr<B, U, Dim, i, j, k, l> &result)
    {
      iter(0, 1, 0, 1) += result(0, 1, 0, 1);
      iter(0, 1, 0, 2) += result(0, 1, 0, 2);
      iter(0, 2, 0, 2) += result(0, 2, 0, 2);
      iter(0, 1, 1, 2) += result(0, 1, 1, 2);
      iter(0, 2, 1, 2) += result(0, 2, 1, 2);
      iter(1, 2, 1, 2) += result(1, 2, 1, 2);
      return *this;
    }

    /* Add a Riemann with the indices switched, making it a
       subtraction. */

    template <class B, class U>
    const Riemann_Expr<Riemann<A, Dim>, T, Dim, i, j, k, l> &
    operator+=(const Riemann_Expr<B, U, Dim, j, i, k, l> &result)
    {
      iter(0, 1, 0, 1) -= result(0, 1, 0, 1);
      iter(0, 1, 0, 2) -= result(0, 1, 0, 2);
      iter(0, 2, 0, 2) -= result(0, 2, 0, 2);
      iter(0, 1, 1, 2) -= result(0, 1, 1, 2);
      iter(0, 2, 1, 2) -= result(0, 2, 1, 2);
      iter(1, 2, 1, 2) -= result(1, 2, 1, 2);
      return *this;
    }
  };
}
