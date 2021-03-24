/* A Tensor0 class that really just ends up being an alias for a
   pointer to a templated type T. Note that the T * is mutable, so the
   pointer can change, allowing iterating over a array.  Also note
   that diffusion and interpolate are included at the end of this
   file, because it needs the full definition of Tensor0. */

#pragma once

#include "Tensor0/dTensor0.hpp"
#include "Tensor0/d_boundary_Tensor0.hpp"
#include "Tensor0/ddTensor0.hpp"
#include "Tensor0/dd_boundary_Tensor0.hpp"

namespace FTensor
{
  template <class T> class Tensor0
  {};

  template <class T> class Tensor0<T *>
  {
    mutable T *__restrict data;

  public:
    Tensor0(T *d) : data(d) {}

    const Tensor0 &operator=(const Tensor0 &a)
    {
      *data = *(a.data);
      return *this;
    }

    template <class U> const Tensor0<T *> &operator=(const U &d)
    {
      *data = d;
      return *this;
    }
    template <class U> const Tensor0<T *> &operator+=(const U &d)
    {
      *data += d;
      return *this;
    }
    template <class U> const Tensor0<T *> &operator-=(const U &d)
    {
      *data -= d;
      return *this;
    }
    template <class U> const Tensor0<T *> &operator*=(const U &d)
    {
      *data *= d;
      return *this;
    }
    template <class U> const Tensor0<T *> &operator/=(const U &d)
    {
      *data /= d;
      return *this;
    }

    /* Note that the conversion operator& to T * only works on
       consts, so it doesn't allow you to change the value of *data.
       You have to use the = operators to change that.  The idea is that
       operator& is only used for stencils and such.  */

    const T *operator&() const { return data; }
    operator T() const { return *data; }

    /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

    const Tensor0<T *> &operator++() const
    {
      ++data;
      return *this;
    }
  };
}

#include "Tensor0/d_one_sided_Tensor0.hpp"
#include "Tensor0/diffusion_Tensor0.hpp"
#include "Tensor0/interpolate_Tensor0.hpp"
