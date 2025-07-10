/* Cross product in 3 dimensions */

#pragma once

#include "Levi_Civita.hpp"

namespace FTensor
{
  template <class A, class B, class T, class U, char i, char j, char k>
  auto cross(const Tensor1_Expr<A, T, 3, i> &a,
             const Tensor1_Expr<B, U, 3, j> &b, const Index<k, 3> &)
  {
    return a * b * levi_civita(Index<i, 3>(), Index<j, 3>(), Index<k, 3>());
  }
}
