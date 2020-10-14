#pragma once

/* A general version, not for pointers. */

#include <ostream>

namespace FTensor
{
  template <class T, int Tensor_Dim> class Tensor2_antisymmetric
  {
    T data[(Tensor_Dim * (Tensor_Dim - 1)) / 2];

  public:
    template <class... U> Tensor2_antisymmetric(U... d) : data{d...}
    {
      static_assert(sizeof...(d) == sizeof(data) / sizeof(T),
                    "Incorrect number of Arguments. Constructor should "
                    "initialize the entire Tensor");
    }

    Tensor2_antisymmetric() {}

    /* There are two ways of accessing the values inside,
       unsafe(int,int) and operator(int,int).  unsafe(int,int) will give
       you a wrong answer if you aren't careful.  The problem is that we
       only store the minimal set of components, but some have different
       signs.  We can't return the negative of a component, and assign
       something to it, because that would assign something to a
       temporary.  To get the correct answer if you don't want to change
       the value, just use operator(int,int). */

    T &unsafe(const int N1, const int N2)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim || N1 < 0 || N2 >= Tensor_Dim || N2 < 0 || N1 >= N2)
        {
          std::stringstream s;
          s << "Bad index in Tensor2_antisymmetric<T," << Tensor_Dim
            << ">.operator(" << N1 << "," << N2 << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[(N2 - 1) + (N1 * (2 * (Tensor_Dim - 1) - N1 - 1)) / 2];
    }

    T operator()(const int N1, const int N2) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim || N1 < 0 || N2 >= Tensor_Dim || N2 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor2_antisymmetric<T," << Tensor_Dim
            << ">.operator(" << N1 << "," << N2 << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 == N2
               ? 0
               : (N1 < N2
                    ? data[(N2 - 1)
                           + (N1 * (2 * (Tensor_Dim - 1) - N1 - 1)) / 2]
                    : -data[(N1 - 1)
                            + (N2 * (2 * (Tensor_Dim - 1) - N2 - 1)) / 2]);
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. They are not entirely safe, since you can accidently use a
       higher dimension than what is really allowed (like Dim=5). */

    /* This returns a Tensor2_Expr, since the indices are not really
       antisymmetric anymore since they cover different dimensions. */

    template <char i, char j, int Dim0, int Dim1>
    typename std::enable_if<(Tensor_Dim >= Dim0 && Tensor_Dim >= Dim1),
                            Tensor2_Expr<Tensor2_antisymmetric<T, Tensor_Dim>,
                                         T, Dim0, Dim1, i, j>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>)
    {
      return Tensor2_Expr<Tensor2_antisymmetric<T, Tensor_Dim>, T, Dim0, Dim1,
                          i, j>(*this);
    }

    template <char i, char j, int Dim0, int Dim1>
    typename std::enable_if<
      (Tensor_Dim >= Dim0 && Tensor_Dim >= Dim1),
      Tensor2_Expr<const Tensor2_antisymmetric<T, Tensor_Dim>, T, Dim0, Dim1,
                   i, j>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>) const
    {
      return Tensor2_Expr<const Tensor2_antisymmetric<T, Tensor_Dim>, T, Dim0,
                          Dim1, i, j>(*this);
    }

    /* This returns a Tensor2_antisymmetric_Expr, since the indices are still
       antisymmetric on the lower dimensions. */

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim >= Dim),
      Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<T, Tensor_Dim>, T, Dim,
                                 i, j>>::type
    operator()(const Index<i, Dim> index1, const Index<j, Dim> index2)
    {
      return Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<T, Tensor_Dim>,
                                        T, Dim, i, j>(*this);
    }

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim >= Dim),
      Tensor2_antisymmetric_Expr<const Tensor2_antisymmetric<T, Tensor_Dim>, T,
                                 Dim, i, j>>::type
    operator()(const Index<i, Dim> index1, const Index<j, Dim> index2) const
    {
      return Tensor2_antisymmetric_Expr<
        const Tensor2_antisymmetric<T, Tensor_Dim>, T, Dim, i, j>(*this);
    }

    /* This is for expressions where a number is used for one slot, and
       an index for another, yielding a Tensor1_Expr.  The non-const
       versions don't actually create a Tensor2_number_rhs_[01] object.
       They create a Tensor1_Expr directly, which provides the
       appropriate indexing operators.  The const versions do create a
       Tensor2_number_[01]. */

    template <char i, int N, int Dim>
    typename std::enable_if<
      (Tensor_Dim >= Dim && Tensor_Dim > N),
      Tensor1_Expr<
        Tensor2_number_rhs_1<Tensor2_antisymmetric<T, Tensor_Dim>, T, N>, T,
        Dim, i>>::type
    operator()(const Index<i, Dim> index1, const Number<N>)
    {
      using TensorExpr
        = Tensor2_number_rhs_1<Tensor2_antisymmetric<T, Tensor_Dim>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N, int Dim>
    typename std::enable_if<
      (Tensor_Dim >= Dim && Tensor_Dim > N),
      Tensor1_Expr<const Tensor2_number_1<
                     const Tensor2_antisymmetric<T, Tensor_Dim>, T, N>,
                   T, Dim, i>>::type
    operator()(const Index<i, Dim> index1, const Number<N>) const
    {
      using TensorExpr
        = Tensor2_number_1<const Tensor2_antisymmetric<T, Tensor_Dim>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    template <char i, int N, int Dim>
    typename std::enable_if<
      (Tensor_Dim > N && Tensor_Dim >= Dim),
      Tensor1_Expr<
        Tensor2_number_rhs_0<Tensor2_antisymmetric<T, Tensor_Dim>, T, N>, T,
        Dim, i>>::type
    operator()(const Number<N>, const Index<i, Dim> index1)
    {
      using TensorExpr
        = Tensor2_number_rhs_0<Tensor2_antisymmetric<T, Tensor_Dim>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N, int Dim>
    typename std::enable_if<
      (Tensor_Dim > N && Tensor_Dim >= Dim),
      Tensor1_Expr<const Tensor2_number_0<
                     const Tensor2_antisymmetric<T, Tensor_Dim>, T, N>,
                   T, Dim, i>>::type
    operator()(const Number<N> &n1, const Index<i, Dim> index1) const
    {
      using TensorExpr
        = Tensor2_number_0<const Tensor2_antisymmetric<T, Tensor_Dim>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Specializations for using actual numbers instead of Number<> */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim >= Dim),
      Tensor1_Expr<
        const Tensor2_numeral_1<const Tensor2_antisymmetric<T, Tensor_Dim>, T>,
        T, Dim, i>>::type
    operator()(const Index<i, Dim> index1, const int N) const
    {
      using TensorExpr
        = Tensor2_numeral_1<const Tensor2_antisymmetric<T, Tensor_Dim>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N));
    }

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim >= Dim),
      Tensor1_Expr<
        const Tensor2_numeral_0<const Tensor2_antisymmetric<T, Tensor_Dim>, T>,
        T, Dim, i>>::type
    operator()(const int N, const Index<i, Dim> index1) const
    {
      using TensorExpr
        = Tensor2_numeral_0<const Tensor2_antisymmetric<T, Tensor_Dim>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N));
    }

    /* These two operator()'s return the Tensor2 with internal
       contractions, yielding a T.  I have to specify one for both
       const and non-const because otherwise the compiler will use the
       operator() which gives a Tensor2_Expr<>. */

    /* TODO Here is a good question. It wont create a problem if i put
       a higher dimension Index in here but it would be
       grammatically wrong. Im gonna add it for now, we can discuss it*/

    template <char i, int Dim>
    typename std::enable_if<(Tensor_Dim >= Dim), T>::type
    operator()(const Index<i, Dim> index1, const Index<i, Dim> index2)
    {
      return 0;
    }

    template <char i, int Dim>
    typename std::enable_if<(Tensor_Dim >= Dim), T>::type
    operator()(const Index<i, Dim> index1, const Index<i, Dim> index2) const
    {
      return 0;
    }
  };
}

/// JSON compatible output.  It only outputs unique, non-zero elments,
/// so a 3x3 antisymmetric matrix only outputs 3 elements.

namespace FTensor
{
  template <class T, int Tensor_Dim>
  std::ostream &Tensor2_antisymmetric_ostream_row(
    std::ostream &os, const FTensor::Tensor2_antisymmetric<T, Tensor_Dim> &t,
    const int &i)
  {
    os << '[';
    for(int j = i + 1; j + 1 < Tensor_Dim; ++j)
      {
        os << t(i, j) << ',';
      }
    if(Tensor_Dim > 0)
      {
        os << t(i, Tensor_Dim - 1);
      }
    os << ']';
    return os;
  }
}

template <class T, int Tensor_Dim>
std::ostream &
operator<<(std::ostream &os,
           const FTensor::Tensor2_antisymmetric<T, Tensor_Dim> &t)
{
  os << '[';
  for(int i = 0; i + 2 < Tensor_Dim; ++i)
    {
      FTensor::Tensor2_antisymmetric_ostream_row(os, t, i);
      os << ',';
    }
  if(Tensor_Dim > 1)
    {
      FTensor::Tensor2_antisymmetric_ostream_row(os, t, Tensor_Dim - 2);
    }
  os << ']';
  return os;
}

namespace FTensor
{
  template <class T, int Tensor_Dim>
  std::istream &Tensor2_antisymmetric_istream_row(
    std::istream &is, FTensor::Tensor2_antisymmetric<T, Tensor_Dim> &t,
    const int &i)
  {
    char c;
    is >> c;
    for(int j = i + 1; j + 1 < Tensor_Dim; ++j)
      {
        is >> t.unsafe(i, j) >> c;
      }
    if(Tensor_Dim > 0)
      {
        is >> t.unsafe(i, Tensor_Dim - 1);
      }
    is >> c;
    return is;
  }
}

template <class T, int Tensor_Dim>
std::istream &
operator>>(std::istream &is, FTensor::Tensor2_antisymmetric<T, Tensor_Dim> &t)
{
  char c;
  is >> c;
  for(int i = 0; i + 2 < Tensor_Dim; ++i)
    {
      FTensor::Tensor2_antisymmetric_istream_row(is, t, i);
      is >> c;
    }
  if(Tensor_Dim > 1)
    {
      FTensor::Tensor2_antisymmetric_istream_row(is, t, Tensor_Dim - 2);
    }
  is >> c;
  return is;
}
