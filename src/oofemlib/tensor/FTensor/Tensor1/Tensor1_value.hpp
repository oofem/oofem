/* The general version, not for pointers. */

#include <iostream>
#pragma once

namespace FTensor
{
  template <class T, int Tensor_Dim> class Tensor1
  {
  protected:
    T data[Tensor_Dim];

  public:
    /* Initializations for varying numbers of elements. */
    template <class... U> constexpr Tensor1(U... d) : data{d...}
    {
      static_assert(sizeof...(d) == sizeof(data) / sizeof(T),
                    "Incorrect number of Arguments. Constructor should "
                    "initialize the entire Tensor");
    };

    constexpr Tensor1() {}

    /* There are two operator(int)'s, one for non-consts that lets you
       change the value, and one for consts that doesn't. */

    T &operator()(const int N)
    {
#ifdef FTENSOR_DEBUG
      if(N >= Tensor_Dim || N < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor1<T," << Tensor_Dim << ">.operator(" << N
            << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N];
    }
    T operator()(const int N) const
    {
#ifdef FTENSOR_DEBUG
      if(N >= Tensor_Dim || N < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor1<T," << Tensor_Dim << ">.operator(" << N
            << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N];
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. They are not entirely safe, since you can accidentaly use a
       higher dimension than what is really allowed (like Dim=5). */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim >= Dim), Tensor1_Expr<Tensor1<T, Tensor_Dim>, T, Dim, i>>::type
    operator()(const Index<i, Dim> &)
    {
      return Tensor1_Expr<Tensor1<T, Tensor_Dim>, T, Dim, i>(*this);
    }

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim >= Dim),
      Tensor1_Expr<const Tensor1<T, Tensor_Dim>, T, Dim, i>>::type
    operator()(const Index<i, Dim> &) const
    {
      return Tensor1_Expr<const Tensor1<T, Tensor_Dim>, T, Dim, i>(*this);
    }

    /* Convenience functions */

    Tensor1<T, Tensor_Dim> normalize()
    {
      const Index<'a', Tensor_Dim> a;
      (*this)(a) /= l2();
      return *this;
    }

    T l2() const { return sqrt(l2_squared(Number<Tensor_Dim>())); }

    template <int Current_Dim> T l2_squared(const Number<Current_Dim> &) const
    {
      return data[Current_Dim - 1] * data[Current_Dim - 1]
             + l2_squared(Number<Current_Dim - 1>());
    }
    T l2_squared(const Number<1> &) const { return data[0] * data[0]; }
  };
}
/// JSON compatible output

template <class T, int Tensor_Dim>
std::ostream &
operator<<(std::ostream &os, const FTensor::Tensor1<T, Tensor_Dim> &t)
{
  os << '[';
  for(int i = 0; i + 1 < Tensor_Dim; ++i)
    {
      os << t(i) << ',';
    }
  if(Tensor_Dim > 0)
    {
      os << t(Tensor_Dim - 1);
    }
  os << ']';
  return os;
}

template <class T, int Tensor_Dim>
std::istream &operator>>(std::istream &is, FTensor::Tensor1<T, Tensor_Dim> &t)
{
  char c;
  is >> c;
  for(int i = 0; i + 1 < Tensor_Dim; ++i)
    {
      is >> t(i) >> c;
    }
  if(Tensor_Dim > 0)
    {
      is >> t(Tensor_Dim - 1);
    }
  is >> c;
  return is;
}
