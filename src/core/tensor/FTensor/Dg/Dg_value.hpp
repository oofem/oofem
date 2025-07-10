/* A general version, not for pointers. */

#pragma once

#include "../Tensor3.hpp"

namespace FTensor
{
  template <class T, int Tensor_Dim01, int Tensor_Dim2> class Dg
  {
    T data[(Tensor_Dim01 * (Tensor_Dim01 + 1)) / 2][Tensor_Dim2];

  public:
    template <class... U> Dg(U... d) : data{d...}
    {
      static_assert(sizeof...(d) == sizeof(data) / sizeof(T),
                    "Incorrect number of Arguments. Constructor should "
                    "initialize the entire Tensor");
    }

    Dg() {}

    /* There are two operator(int,int,int)'s, one for non-consts that lets you
       change the value, and one for consts that doesn't. */

    T &operator()(const int N1, const int N2, const int N3)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim01 || N1 < 0 || N2 >= Tensor_Dim01 || N2 < 0
         || N3 >= Tensor_Dim2 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Dg<T," << Tensor_Dim01 << "," << Tensor_Dim2
            << ">.operator(" << N1 << "," << N2 << "," << N3 << ")"
            << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 > N2 ? data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2][N3]
                     : data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2][N3];
    }

    T operator()(const int N1, const int N2, const int N3) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim01 || N1 < 0 || N2 >= Tensor_Dim01 || N2 < 0
         || N3 >= Tensor_Dim2 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Dg<T," << Tensor_Dim01 << "," << Tensor_Dim2
            << ">.operator(" << N1 << "," << N2 << "," << N3 << ") const"
            << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 > N2 ? data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2][N3]
                     : data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2][N3];
    }

    /* These operator()'s are the first part in constructing template
       expressions. */

    template <char i, char j, char k, int Dim01, int Dim2>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim01 && Tensor_Dim2 >= Dim2),
      Dg_Expr<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k>>::type
    operator()(const Index<i, Dim01>, const Index<j, Dim01>,
               const Index<k, Dim2>)
    {
      return Dg_Expr<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k>(
        *this);
    }
    // TODO Add support for multiple dimensions ^^^^^
    template <char i, char j, char k, int Dim01, int Dim2>
    typename std::enable_if<(Tensor_Dim01 >= Dim01 && Tensor_Dim2 >= Dim2),
                            Dg_Expr<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T,
                                    Dim01, Dim2, i, j, k>>::type
    operator()(const Index<i, Dim01>, const Index<j, Dim01>,
               const Index<k, Dim2>) const
    {
      return Dg_Expr<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i,
                     j, k>(*this);
    }

    /* These operators are for internal contractions. We can not do
       something like A(i,j,j) where i and j have different dimensions,
       because it becomes ambiguous. */

    /* A(i,j,j) */

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim && Tensor_Dim2 >= Dim),
      Tensor1_Expr<
        Tensor3_contracted_12<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim>, T,
        Dim, i>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>,
               const Index<j, Dim>) const
    {
      using TensorExpr
        = Tensor3_contracted_12<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* A(j,i,j) */

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim && Tensor_Dim2 >= Dim),
      Tensor1_Expr<
        Tensor3_contracted_02<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim>, T,
        Dim, i>>::type
    operator()(const Index<j, Dim>, const Index<i, Dim>,
               const Index<j, Dim>) const
    {
      using TensorExpr
        = Tensor3_contracted_02<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* A(j,j,i) */

    template <char i, char j, int Dim01, int Dim2>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim01 && Tensor_Dim2 >= Dim2),
      Tensor1_Expr<
        Tensor3_contracted_01<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01>, T,
        Dim2, i>>::type
    operator()(const Index<j, Dim01>, const Index<j, Dim01>,
               const Index<i, Dim2>) const
    {
      using TensorExpr
        = Tensor3_contracted_01<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01>;
      return Tensor1_Expr<TensorExpr, T, Dim2, i>(TensorExpr(*this));
    }

    /* This is for expressions where a number is used for one slot, and
       indices for the others, yielding a Tensor2_Expr or
       Tensor2_symmetric_Expr.  The non-const versions don't actually
       create a Dg_number_rhs_* object, while the const versions
       do create a Dg_number_*. */

    /* First slot. */

    template <char i, char j, int N, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim01 > N && Tensor_Dim01 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Dg_number_rhs_0<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>, T,
                   Dim1, Dim2, i, j>>::type
    operator()(const Number<N>, const Index<i, Dim1>, const Index<j, Dim2>)
    {
      using TensorExpr
        = Dg_number_rhs_0<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim2, i, j>(*this);
    }

    template <char i, char j, int N, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim01 > N && Tensor_Dim01 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Dg_number_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>,
                   T, Dim1, Dim2, i, j>>::type
    operator()(const Number<N>, const Index<i, Dim1>,
               const Index<j, Dim2>) const
    {
      using TensorExpr
        = Dg_number_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim2, i, j>(TensorExpr(*this));
    }

    /* Second slot. */

    template <char i, char j, int N, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim0 && Tensor_Dim01 > N && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Dg_number_rhs_0<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>, T,
                   Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const Number<N>, const Index<j, Dim2>)
    {
      using TensorExpr
        = Dg_number_rhs_0<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(*this);
    }

    template <char i, char j, int N, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim0 && Tensor_Dim01 > N && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Dg_number_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>,
                   T, Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const Number<N>,
               const Index<j, Dim2>) const
    {
      using TensorExpr
        = Dg_number_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(TensorExpr(*this));
    }

    /* Third slot. */

    template <char i, char j, int N, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim && Tensor_Dim2 > N),
      Tensor2_symmetric_Expr<
        Dg_number_rhs_2<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>, T, Dim, i,
        j>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>, const Number<N>)
    {
      using TensorExpr
        = Dg_number_rhs_2<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(*this);
    }

    template <char i, char j, int N, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim && Tensor_Dim2 > N),
      Tensor2_symmetric_Expr<
        Dg_number_2<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>, T, Dim, i,
        j>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>, const Number<N>) const
    {
      using TensorExpr
        = Dg_number_2<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this));
    }
    // TODO add multiple dimensions up here ^^^^^
    /* This is for expressions where a number is used for two slots, and
       an Index for the other, yielding a Tensor1_Expr.  The non-const
       versions don't actually create a Dg_number_rhs_* object,
       while the const versions do create a Dg_number_*. */

    /* Index in first slot. */

    template <char i, int N1, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim && Tensor_Dim01 > N1 && Tensor_Dim2 > N2),
      Tensor1_Expr<Dg_number_rhs_12<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>,
                   T, Dim, i>>::type
    operator()(const Index<i, Dim> index, const Number<N1>, const Number<N2>)
    {
      using TensorExpr
        = Dg_number_rhs_12<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N1, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim && Tensor_Dim01 > N1 && Tensor_Dim2 > N2),
      Tensor1_Expr<
        Dg_number_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>, T,
        Dim, i>>::type
    operator()(const Index<i, Dim> index, const Number<N1>,
               const Number<N2>) const
    {
      using TensorExpr
        = Dg_number_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Index in second slot.  I use the same structures as for the Index
       in the first slot since the tensor is symmetric on the first two
       indices. */

    template <char i, int N0, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 > N0 && Tensor_Dim01 >= Dim && Tensor_Dim2 > N2),
      Tensor1_Expr<Dg_number_rhs_12<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N0, N2>,
                   T, Dim, i>>::type
    operator()(const Number<N0>, const Index<i, Dim>, const Number<N2>)
    {
      using TensorExpr
        = Dg_number_rhs_12<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N0, N2>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N0, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 > N0 && Tensor_Dim01 >= Dim && Tensor_Dim2 > N2),
      Tensor1_Expr<
        Dg_number_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N0, N2>, T,
        Dim, i>>::type
    operator()(const Number<N0>, const Index<i, Dim> index,
               const Number<N2>) const
    {
      using TensorExpr
        = Dg_number_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N0, N2>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Index in third slot. */

    template <char i, int N0, int N1, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 > N0 && Tensor_Dim01 > N1 && Tensor_Dim2 >= Dim),
      Tensor1_Expr<Dg_number_rhs_01<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N0, N1>,
                   T, Dim, i>>::type
    operator()(const Number<N0>, const Number<N1>, const Index<i, Dim> index)
    {
      using TensorExpr
        = Dg_number_rhs_01<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N0, N1>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N0, int N1, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 > N0 && Tensor_Dim01 > N1 && Tensor_Dim2 >= Dim),
      Tensor1_Expr<
        Dg_number_01<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N0, N1>, T,
        Dim, i>>::type
    operator()(const Number<N0>, const Number<N1>, const Index<i, Dim>) const
    {
      using TensorExpr
        = Dg_number_01<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N0, N1>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Specializations for using actual numbers instead of Number<>.
       This is for expressions where an actual number is used for one
       slot, and indices for the others, yielding a Tensor2_Expr or
       Tensor2_symmetric_Expr. I only define the const versions. */

    /* First slot. */

    template <char i, char j, int Dim1, int Dim2>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim1 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Dg_numeral_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T,
                   Dim1, Dim2, i, j>>::type
    operator()(const int N, const Index<i, Dim1>, const Index<j, Dim2>) const
    {
      using TensorExpr
        = Dg_numeral_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim2, i, j>(
        TensorExpr(*this, N));
    }

    /* Second slot. */

    template <char i, char j, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim0 && Tensor_Dim2 >= Dim2),
      Tensor2_Expr<Dg_numeral_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T,
                   Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const int N, const Index<j, Dim2>) const
    {
      using TensorExpr
        = Dg_numeral_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(
        TensorExpr(*this, N));
    }

    /* Third slot. */

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim),
      Tensor2_symmetric_Expr<
        Dg_numeral_2<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T, Dim, i,
        j>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>, const int N) const
    {
      using TensorExpr
        = Dg_numeral_2<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this, N));
    }
    // TODO Allow two different dimensions in 1 and 2
    /* This is for expressions where a numeral is used for two slots, and
       an Index for the other, yielding a Tensor1_Expr. */

    /* Index in first slot. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim),
      Tensor1_Expr<Dg_numeral_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T,
                   Dim, i>>::type
    operator()(const Index<i, Dim> index, const int N1, const int N2) const
    {
      using TensorExpr
        = Dg_numeral_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
    }

    /* Index in second slot.  I use the same structures as for the Index
       in the first slot since the tensor is symmetric on the first two
       indices. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim),
      Tensor1_Expr<Dg_numeral_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T,
                   Dim, i>>::type
    operator()(const int N1, const Index<i, Dim> index, const int N2) const
    {
      using TensorExpr
        = Dg_numeral_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
    }

    /* Index in third slot. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim2 >= Dim),
      Tensor1_Expr<Dg_numeral_01<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T,
                   Dim, i>>::type
    operator()(const int N1, const int N2, const Index<i, Dim>) const
    {
      using TensorExpr
        = Dg_numeral_01<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
    }
  };
}
