/* Declaration for Riemann, which has the symmetries of the
   Riemann tensor.  Antisymmetric on the first two indices and the
   last two indices, and symmetric on a cyclic permutation of the last
   three indices. */

#pragma once

namespace FTensor
{
  template <class T, int Dim> class Riemann
  {};

  template <class T> class Riemann<T, 3>
  {
    /* The zero variable is because some of the components of a tensor
       with these symmetries are identically zero. */

    T zero;
    T data0101, data0102, data0112, data0202, data0212, data1212;

  public:
    Riemann() : zero(0.0) {}

    /* There are two operator(int,int,int,int)'s, one for non-consts
       that lets you change the value, and one for consts that doesn't.
       The non-const one will give you a wrong answer if you aren't
       careful.  The problem is that we only store the minimal set of
       components, but some have different signs.  We can't return the
       negative of a component, and assign something to it, because that
       would assign something to a temporary.  To get the correct answer
       if you don't want to change the value, use eval instead. */

    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return N1 == 0
               ? (N2 == 0
                    ? zero
                    : (N2 == 1
                         ? (N3 == 0
                              ? (N4 == 0 ? zero
                                         : (N4 == 1 ? data0101 : data0102))
                              : (N3 == 1
                                   ? (N4 == 0 ? -data0101
                                              : (N4 == 1 ? zero : data0112))
                                   : (N4 == 0 ? -data0102
                                              : (N4 == 1 ? -data0112 : zero))))
                         : (N3 == 0
                              ? (N4 == 0 ? zero
                                         : (N4 == 1 ? data0102 : data0202))
                              : (N3 == 1
                                   ? (N4 == 0 ? -data0102
                                              : (N4 == 1 ? zero : data0212))
                                   : (N4 == 0
                                        ? -data0202
                                        : (N4 == 1 ? -data0212 : zero))))))
               : (N1 == 1
                    ? (N2 == 0
                         ? (N3 == 0
                              ? (N4 == 0 ? zero
                                         : (N4 == 1 ? -data0101 : -data0102))
                              : (N3 == 1
                                   ? (N4 == 0 ? data0101
                                              : (N4 == 1 ? zero : -data0112))
                                   : (N4 == 0 ? data0102
                                              : (N4 == 1 ? data0112 : zero))))
                         : (N2 == 1
                              ? zero
                              : (N3 == 0
                                   ? (N4 == 0 ? zero
                                              : (N4 == 1 ? data0112 : data0212))
                                   : (N3 == 1
                                        ? (N4 == 0
                                             ? -data0112
                                             : (N4 == 1 ? zero : data1212))
                                        : (N4 == 0 ? -data0212
                                                   : (N4 == 1 ? -data1212
                                                              : zero))))))
                    : (N2 == 0
                         ? (N3 == 0
                              ? (N4 == 0 ? zero
                                         : (N4 == 1 ? -data0102 : -data0202))
                              : (N3 == 1
                                   ? (N4 == 0 ? data0102
                                              : (N4 == 1 ? zero : -data0212))
                                   : (N4 == 0 ? data0202
                                              : (N4 == 1 ? data0212 : zero))))
                         : (N2 == 1
                              ? (N3 == 0
                                   ? (N4 == 0
                                        ? zero
                                        : (N4 == 1 ? -data0112 : -data0212))
                                   : (N3 == 1
                                        ? (N4 == 0
                                             ? data0112
                                             : (N4 == 1 ? zero : -data1212))
                                        : (N4 == 0
                                             ? data0212
                                             : (N4 == 1 ? data1212 : zero))))
                              : zero)));
    }

    T eval(const int N1, const int N2, const int N3, const int N4) const
    {
      return N1 == 0
               ? (N2 == 0
                    ? zero
                    : (N2 == 1
                         ? (N3 == 0
                              ? (N4 == 0 ? zero
                                         : (N4 == 1 ? data0101 : data0102))
                              : (N3 == 1
                                   ? (N4 == 0 ? -data0101
                                              : (N4 == 1 ? zero : data0112))
                                   : (N4 == 0 ? -data0102
                                              : (N4 == 1 ? -data0112 : zero))))
                         : (N3 == 0
                              ? (N4 == 0 ? zero
                                         : (N4 == 1 ? data0102 : data0202))
                              : (N3 == 1
                                   ? (N4 == 0 ? -data0102
                                              : (N4 == 1 ? zero : data0212))
                                   : (N4 == 0
                                        ? -data0202
                                        : (N4 == 1 ? -data0212 : zero))))))
               : (N1 == 1
                    ? (N2 == 0
                         ? (N3 == 0
                              ? (N4 == 0 ? zero
                                         : (N4 == 1 ? -data0101 : -data0102))
                              : (N3 == 1
                                   ? (N4 == 0 ? data0101
                                              : (N4 == 1 ? zero : -data0112))
                                   : (N4 == 0 ? data0102
                                              : (N4 == 1 ? data0112 : zero))))
                         : (N2 == 1
                              ? zero
                              : (N3 == 0
                                   ? (N4 == 0 ? zero
                                              : (N4 == 1 ? data0112 : data0212))
                                   : (N3 == 1
                                        ? (N4 == 0
                                             ? -data0112
                                             : (N4 == 1 ? zero : data1212))
                                        : (N4 == 0 ? -data0212
                                                   : (N4 == 1 ? -data1212
                                                              : zero))))))
                    : (N2 == 0
                         ? (N3 == 0
                              ? (N4 == 0 ? zero
                                         : (N4 == 1 ? -data0102 : -data0202))
                              : (N3 == 1
                                   ? (N4 == 0 ? data0102
                                              : (N4 == 1 ? zero : -data0212))
                                   : (N4 == 0 ? data0202
                                              : (N4 == 1 ? data0212 : zero))))
                         : (N2 == 1
                              ? (N3 == 0
                                   ? (N4 == 0
                                        ? zero
                                        : (N4 == 1 ? -data0112 : -data0212))
                                   : (N3 == 1
                                        ? (N4 == 0
                                             ? data0112
                                             : (N4 == 1 ? zero : -data1212))
                                        : (N4 == 0
                                             ? data0212
                                             : (N4 == 1 ? data1212 : zero))))
                              : zero)));
    }

    T &operator()(const int N1, const int N2, const int N3, const int N4)
    {
      return (N1 == 0 && N2 == 1 && N3 == 0 && N4 == 1)
               ? data0101
               : ((N1 == 0 && N2 == 1 && N3 == 0 && N4 == 2)
                    ? data0102
                    : ((N1 == 0 && N2 == 1 && N3 == 1 && N4 == 2)
                         ? data0112
                         : ((N1 == 0 && N2 == 2 && N3 == 1 && N4 == 2)
                              ? data0212
                              : ((N1 == 0 && N2 == 2 && N3 == 0 && N4 == 2)
                                   ? data0202
                                   : ((N1 == 1 && N2 == 2 && N3 == 1 && N4 == 2)
                                        ? data1212
                                        : zero)))));
    }

    /* These operator()'s are the first part in constructing template
       expressions. */

    template <char i, char j, char k, char l>
    Riemann_Expr<Riemann<T, 3>, T, 3, i, j, k, l>
    operator()(const Index<i, 3>, const Index<j, 3>, const Index<k, 3>,
               const Index<l, 3>)
    {
      return Riemann_Expr<Riemann<T, 3>, T, 3, i, j, k, l>(*this);
    }
  };
}
#include "Riemann/Riemann_Expr.hpp"
