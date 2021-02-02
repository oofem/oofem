/* A general version, not for pointers. */

#pragma once

#include "Tensor3_contracted.hpp"

#include <ostream>

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  class Tensor3
  {
    T data[Tensor_Dim0][Tensor_Dim1][Tensor_Dim2];

  public:
    template <class... U> Tensor3(U... d) : data{d...}
    {
      static_assert(sizeof...(d) == sizeof(data) / sizeof(T),
                    "Incorrect number of Arguments. Constructor should "
                    "initialize the entire Tensor");
    }

    Tensor3() {}

    /* There are two operator(int,int,int)'s, one for non-consts that lets you
       change the value, and one for consts that doesn't. */

    T &operator()(const int N1, const int N2, const int N3)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0
         || N3 >= Tensor_Dim2 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor3<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << "," << Tensor_Dim2 << ">.operator(" << N1 << "," << N2 << ","
            << N3 << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N1][N2][N3];
    }

    T operator()(const int N1, const int N2, const int N3) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0
         || N3 >= Tensor_Dim2 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor3<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << "," << Tensor_Dim2 << ">.operator(" << N1 << "," << N2 << ","
            << N3 << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N1][N2][N3];
    }

    /* These operator()'s are the first part in constructing template
       expressions. */

    template <char i, char j, char k, int Dim0, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor3_Expr<Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim0,
                   Dim1, Dim2, i, j, k>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>,
               const Index<k, Dim2>)
    {
      return Tensor3_Expr<Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T,
                          Dim0, Dim1, Dim2, i, j, k>(*this);
    }

    template <char i, char j, char k, int Dim0, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor3_Expr<const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T,
                   Dim0, Dim1, Dim2, i, j, k>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>,
               const Index<k, Dim2>) const
    {
      return Tensor3_Expr<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim0, Dim1,
        Dim2, i, j, k>(*this);
    }

    /* These operators are for internal contractions. We can not do
       something like A(i,j,j) where i and j have different dimensions,
       because it becomes ambiguous. */
    /*TODO I dont know how i and j having different dimensions can be
     * ambiguous?? Since you do Use that in A(j,j,i)*/
    /* A(i,j,j) */

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim && Tensor_Dim2 >= Dim),
      Tensor1_Expr<
        Tensor3_contracted_12<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim>,
        T, Dim, i>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>,
               const Index<j, Dim>) const
    {
      using TensorExpr = Tensor3_contracted_12<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* A(j,i,j) */

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim && Tensor_Dim2 >= Dim),
      Tensor1_Expr<
        Tensor3_contracted_02<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim>,
        T, Dim, i>>::type
    operator()(const Index<j, Dim>, const Index<i, Dim>,
               const Index<j, Dim>) const
    {
      using TensorExpr = Tensor3_contracted_02<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* A(j,j,i) */

    template <char i, char j, int Dim, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim && Tensor_Dim2 >= Dim2),
      Tensor1_Expr<
        Tensor3_contracted_01<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim>,
        T, Dim2, i>>::type
    operator()(const Index<j, Dim>, const Index<j, Dim>,
               const Index<i, Dim2>) const
    {
      using TensorExpr = Tensor3_contracted_01<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, Dim2>;
      return Tensor1_Expr<TensorExpr, T, Dim2, i>(TensorExpr(*this));
    }

    /* This is for expressions where a number is used for one slot, and
       indices for the others, yielding a Tensor2_Expr or
       Tensor2_symmetric_Expr.  The non-const versions don't actually
       create a Tensor3_number_rhs_* object, while the const versions
       do create a Tensor3_number_*. */

    /* First slot. */

    template <char i, char j, int N, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 > N && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Tensor3_number_rhs_0<
                     Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>,
                   T, Dim1, Dim2, i, j>>::type
    operator()(const Number<N>, const Index<i, Dim1>, const Index<j, Dim2>)
    {
      using TensorExpr = Tensor3_number_rhs_0<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim2, i, j>(*this);
    }

    template <char i, char j, int N, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 > N && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<
        Tensor3_number_0<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>,
        T, Dim1, Dim2, i, j>>::type
    operator()(const Number<N>, const Index<i, Dim1>,
               const Index<j, Dim2>) const
    {
      using TensorExpr = Tensor3_number_0<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim2, i, j>(TensorExpr(*this));
    }

    /* Second slot. */

    template <char i, char j, int N, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 > N && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Tensor3_number_rhs_0<
                     Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>,
                   T, Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const Number<N>, const Index<j, Dim2>)
    {
      using TensorExpr = Tensor3_number_rhs_0<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(*this);
    }

    template <char i, char j, int N, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 > N && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<
        Tensor3_number_0<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>,
        T, Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const Number<N>,
               const Index<j, Dim2>) const
    {
      using TensorExpr = Tensor3_number_0<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(TensorExpr(*this));
    }

    /* Third slot. */
    // TODO allow two different dimensions here
    template <char i, char j, int N, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim && Tensor_Dim2 > N),
      Tensor2_symmetric_Expr<
        Tensor3_number_rhs_2<Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>,
                             T, N>,
        T, Dim, i, j>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>, const Number<N>)
    {
      using TensorExpr = Tensor3_number_rhs_2<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(*this);
    }
    // TODO allow two different dimensions here
    template <char i, char j, int N, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim && Tensor_Dim2 > N),
      Tensor2_symmetric_Expr<
        Tensor3_number_2<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>,
        T, Dim, i, j>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>, const Number<N>) const
    {
      using TensorExpr = Tensor3_number_2<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this));
    }

    /* This is for expressions where a number is used for two slots, and
       an Index for the other, yielding a Tensor1_Expr.  The non-const
       versions don't actually create a Tensor3_number_rhs_* object,
       while the const versions do create a Tensor3_number_*. */

    /* Index in first slot. */

    template <char i, int N1, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 > N1 && Tensor_Dim2 > N2),
      Tensor1_Expr<
        Tensor3_number_rhs_12<
          Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N1, N2>,
        T, Dim, i>>::type
    operator()(const Index<i, Dim> index, const Number<N1>, const Number<N2>)
    {
      using TensorExpr = Tensor3_number_rhs_12<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N1, N2>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N1, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 > N1 && Tensor_Dim2 > N2),
      Tensor1_Expr<
        Tensor3_number_12<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N1, N2>,
        T, Dim, i>>::type
    operator()(const Index<i, Dim> index, const Number<N1>,
               const Number<N2>) const
    {
      using TensorExpr = Tensor3_number_12<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N1, N2>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Index in second slot.  I use the same structures as for the Index
       in the first slot since the tensor is symmetric on the first two
       indices. */

    template <char i, int N0, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 > N0 && Tensor_Dim1 >= Dim && Tensor_Dim2 > N2),
      Tensor1_Expr<
        Tensor3_number_rhs_12<
          Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N2>,
        T, Dim, i>>::type
    operator()(const Number<N0>, const Index<i, Dim> index, const Number<N2>)
    {
      using TensorExpr = Tensor3_number_rhs_12<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N2>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N0, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 > N0 && Tensor_Dim1 >= Dim && Tensor_Dim2 > N2),
      Tensor1_Expr<
        Tensor3_number_12<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N2>,
        T, Dim, i>>::type
    operator()(const Number<N0>, const Index<i, Dim> index,
               const Number<N2>) const
    {
      using TensorExpr = Tensor3_number_12<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N2>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Index in third slot. */

    template <char i, int N0, int N1, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 > N0 && Tensor_Dim1 > N1 && Tensor_Dim2 >= Dim),
      Tensor1_Expr<
        Tensor3_number_rhs_01<
          Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N1>,
        T, Dim, i>>::type
    operator()(const Number<N0>, const Number<N1>, const Index<i, Dim> index)
    {
      using TensorExpr = Tensor3_number_rhs_01<
        Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N1>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N0, int N1, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 > N0 && Tensor_Dim1 > N1 && Tensor_Dim2 >= Dim),
      Tensor1_Expr<
        Tensor3_number_01<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N1>,
        T, Dim, i>>::type
    operator()(const Number<N0>, const Number<N1>,
               const Index<i, Dim> index) const
    {
      using TensorExpr = Tensor3_number_01<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T, N0, N1>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Specializations for using actual numbers instead of Number<>.
       This is for expressions where an actual number is used for one
       slot, and indices for the others, yielding a Tensor2_Expr or
       Tensor2_symmetric_Expr. I only define the const versions. */

    /* First slot. */

    template <char i, char j, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Tensor3_numeral_0<
                     const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
                   T, Dim1, Dim2, i, j>>::type
    operator()(const int N, const Index<i, Dim1>, const Index<j, Dim2>) const
    {
      using TensorExpr = Tensor3_numeral_0<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim2, i, j>(
        TensorExpr(*this, N));
    }

    /* Second slot. */

    template <char i, char j, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Tensor3_numeral_0<
                     const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
                   T, Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const int N, const Index<j, Dim2>) const
    {
      using TensorExpr = Tensor3_numeral_0<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(
        TensorExpr(*this, N));
    }

    /* Third slot. */
    // TODO Allow multiple dimensions here
    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim1 >= Dim),
      Tensor2_symmetric_Expr<
        Tensor3_numeral_2<
          const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
        T, Dim, i, j>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>, const int N) const
    {
      using TensorExpr = Tensor3_numeral_2<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this, N));
    }

    /* This is for expressions where a numeral is used for two slots, and
       an Index for the other, yielding a Tensor1_Expr. */

    /* Index in first slot. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim),
      Tensor1_Expr<Tensor3_numeral_12<
                     const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
                   T, Dim, i>>::type
    operator()(const Index<i, Dim> index, const int N1, const int N2) const
    {
      using TensorExpr = Tensor3_numeral_12<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
    }

    /* Index in second slot.  I use the same structures as for the Index
       in the first slot since the tensor is symmetric on the first two
       indices. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim1 >= Dim),
      Tensor1_Expr<Tensor3_numeral_12<
                     const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
                   T, Dim, i>>::type
    operator()(const int N1, const Index<i, Dim> index, const int N2) const
    {
      using TensorExpr = Tensor3_numeral_12<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
    }

    /* Index in third slot. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim2 >= Dim),
      Tensor1_Expr<Tensor3_numeral_01<
                     const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>,
                   T, Dim, i>>::type
    operator()(const int N1, const int N2, const Index<i, Dim> index) const
    {
      using TensorExpr = Tensor3_numeral_01<
        const Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
    }
  };
}

/// JSON compatible output for numbers.  If T is something funky
/// (strings, arrays, etc.) then you are going to have a bad time.

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  std::ostream &Tensor3_ostream_row(
    std::ostream &os,
    const FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t,
    const int &i, const int &j)
  {
    os << '[';
    for(int k = 0; k + 1 < Tensor_Dim2; ++k)
      {
        os << t(i, j, k) << ',';
      }
    if(Tensor_Dim2 > 0)
      {
        os << t(i, j, Tensor_Dim2 - 1);
      }
    os << ']';
    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  std::ostream &Tensor3_ostream_block(
    std::ostream &os,
    const FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t,
    const int &i)
  {
    os << '[';
    for(int j = 0; j + 1 < Tensor_Dim1; ++j)
      {
        FTensor::Tensor3_ostream_row(os, t, i, j);
        os << ',';
      }
    if(Tensor_Dim1 > 0)
      {
        FTensor::Tensor3_ostream_row(os, t, i, Tensor_Dim1 - 1);
      }
    os << ']';
    return os;
  }
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
std::ostream &
operator<<(std::ostream &os,
           const FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t)
{
  os << '[';
  for(int i = 0; i + 1 < Tensor_Dim0; ++i)
    {
      FTensor::Tensor3_ostream_block(os, t, i);
      os << ',';
    }
  if(Tensor_Dim0 > 0)
    {
      FTensor::Tensor3_ostream_block(os, t, Tensor_Dim0 - 1);
    }
  os << ']';
  return os;
}

/// JSON compatible input.  It does not validate that separators are
/// actually braces '[' or commas ','.  It also ignores errors from
/// missing trailing characters.  So you could do something like
///
///   Tensor3<double,2,2,2> t3_1;
///   std::stringstream ss(":::3:4:::7:8:::::11:12:::13:14");
///   ss >> t3_1;

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  std::istream &Tensor3_istream_row(
    std::istream &is,
    FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t,
    const int &i, const int &j)
  {
    char c;
    is >> c;
    for(int k = 0; k + 1 < Tensor_Dim2; ++k)
      {
        is >> t(i, j, k) >> c;
      }
    if(Tensor_Dim2 > 0)
      {
        is >> t(i, j, Tensor_Dim2 - 1);
      }
    is >> c;
    return is;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  std::istream &Tensor3_istream_block(
    std::istream &is,
    FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t,
    const int &i)
  {
    char c;
    is >> c;
    for(int j = 0; j + 1 < Tensor_Dim1; ++j)
      {
        Tensor3_istream_row(is, t, i, j);
        is >> c;
      }
    if(Tensor_Dim1 > 0)
      {
        Tensor3_istream_row(is, t, i, Tensor_Dim1 - 1);
      }
    is >> c;
    return is;
  }
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
std::istream &
operator>>(std::istream &is,
           FTensor::Tensor3<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2> &t)
{
  char c;
  is >> c;
  for(int i = 0; i + 1 < Tensor_Dim0; ++i)
    {
      FTensor::Tensor3_istream_block(is, t, i);
      is >> c;
    }
  if(Tensor_Dim0 > 0)
    {
      FTensor::Tensor3_istream_block(is, t, Tensor_Dim0 - 1);
    }
  is >> c;
  return is;
}
