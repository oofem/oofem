/* A version for pointers */

#pragma once

#include "../Tensor3.hpp"

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim12>
  class Christof<T *, Tensor_Dim0, Tensor_Dim12>
  {
    mutable T *restrict data[Tensor_Dim0][(Tensor_Dim12 * (Tensor_Dim12 + 1)) / 2];

  public:
    template <class... U> Christof(U *... d) : data{d...}
    {
      static_assert(sizeof...(d) == sizeof(data) / sizeof(T),
                    "Incorrect number of Arguments. Constructor should "
                    "initialize the entire Tensor");
    }

    Christof() {}

    /* There are two operator(int,int,int)'s, one for non-consts that lets you
       change the value, and one for consts that doesn't. */

    T &operator()(const int N1, const int N2, const int N3)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim12 || N2 < 0
         || N3 >= Tensor_Dim12 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Christof<T*," << Tensor_Dim0 << ","
            << Tensor_Dim12 << ">.operator(" << N1 << "," << N2 << "," << N3
            << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N2 > N3 ? *data[N1][N2 + (N3 * (2 * Tensor_Dim12 - N3 - 1)) / 2]
                     : *data[N1][N3 + (N2 * (2 * Tensor_Dim12 - N2 - 1)) / 2];
    }

    T operator()(const int N1, const int N2, const int N3) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim12 || N2 < 0
         || N3 >= Tensor_Dim12 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Christof<T*," << Tensor_Dim0 << ","
            << Tensor_Dim12 << ">.operator(" << N1 << "," << N2 << "," << N3
            << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N2 > N3 ? *data[N1][N2 + (N3 * (2 * Tensor_Dim12 - N3 - 1)) / 2]
                     : *data[N1][N3 + (N2 * (2 * Tensor_Dim12 - N2 - 1)) / 2];
    }

    T *ptr(const int N1, const int N2, const int N3) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim12 || N2 < 0
         || N3 >= Tensor_Dim12 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Christof<T*," << Tensor_Dim0 << ","
            << Tensor_Dim12 << ">.ptr(" << N1 << "," << N2 << "," << N3 << ")"
            << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N2 > N3 ? data[N1][N2 + (N3 * (2 * Tensor_Dim12 - N3 - 1)) / 2]
                     : data[N1][N3 + (N2 * (2 * Tensor_Dim12 - N2 - 1)) / 2];
    }

    /* These operator()'s are the first part in constructing template
       expressions.  I mix up the indices here so that it behaves like a
       Dg.  That way I don't have to have a separate wrapper
       class Christof_Expr, which simplifies things. */

    template <char i, char j, char k, int Dim0, int Dim12>
    Dg_Expr<Christof<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i, j, k>
    operator()(const Index<k, Dim0> index1, const Index<i, Dim12> index2,
               const Index<j, Dim12> index3)
    {
      return Dg_Expr<Christof<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0,
                     i, j, k>(*this);
    }

    template <char i, char j, char k, int Dim0, int Dim12>
    Dg_Expr<const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i,
            j, k>
    operator()(const Index<k, Dim0> index1, const Index<i, Dim12> index2,
               const Index<j, Dim12> index3) const
    {
      return Dg_Expr<const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim12,
                     Dim0, i, j, k>(*this);
    }

    /* These operators are for internal contractions. */

    /* const versions */

    template <char i, char j, int Dim0, int Dim12>
    Tensor1_Expr<const Tensor3_contracted_12<
                   const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim12>,
                 T, Dim0, i>
    operator()(const Index<i, Dim0> index1, const Index<j, Dim12> index2,
               const Index<j, Dim12> index3) const
    {
      using TensorExpr
        = Tensor3_contracted_12<const Christof<T *, Tensor_Dim0, Tensor_Dim12>,
                                T, Dim12>;
      return Tensor1_Expr<TensorExpr, T, Dim0, i>(TensorExpr(*this));
    }

    template <char i, char j, int Dim02, int Dim1>
    Tensor1_Expr<const Tensor3_contracted_02<
                   const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim02>,
                 T, Dim1, i>
    operator()(const Index<j, Dim02> index1, const Index<i, Dim1> index2,
               const Index<j, Dim02> index3) const
    {
      using TensorExpr
        = Tensor3_contracted_02<const Christof<T *, Tensor_Dim0, Tensor_Dim12>,
                                T, Dim02>;
      return Tensor1_Expr<TensorExpr, T, Dim1, i>(TensorExpr(*this));
    }

    template <char i, char j, int Dim01, int Dim2>
    Tensor1_Expr<const Tensor3_contracted_01<
                   const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim01>,
                 T, Dim2, i>
    operator()(const Index<j, Dim01> index1, const Index<j, Dim01> index2,
               const Index<i, Dim2> index3) const
    {
      using TensorExpr
        = Tensor3_contracted_01<const Christof<T *, Tensor_Dim0, Tensor_Dim12>,
                                T, Dim01>;
      return Tensor1_Expr<TensorExpr, T, Dim2, i>(TensorExpr(*this));
    }

    /* non-const versions */

    template <char i, char j, int Dim0, int Dim12>
    Tensor1_Expr<const Tensor3_contracted_12<
                   const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim12>,
                 T, Dim0, i>
    operator()(const Index<i, Dim0> index1, const Index<j, Dim12> index2,
               const Index<j, Dim12> index3)
    {
      using TensorExpr
        = Tensor3_contracted_12<const Christof<T *, Tensor_Dim0, Tensor_Dim12>,
                                T, Dim12>;
      return Tensor1_Expr<TensorExpr, T, Dim0, i>(TensorExpr(*this));
    }

    template <char i, char j, int Dim02, int Dim1>
    Tensor1_Expr<const Tensor3_contracted_02<
                   const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim02>,
                 T, Dim1, i>
    operator()(const Index<j, Dim02> index1, const Index<i, Dim1> index2,
               const Index<j, Dim02> index3)
    {
      using TensorExpr
        = Tensor3_contracted_02<const Christof<T *, Tensor_Dim0, Tensor_Dim12>,
                                T, Dim02>;
      return Tensor1_Expr<TensorExpr, T, Dim1, i>(TensorExpr(*this));
    }

    template <char i, char j, int Dim01, int Dim2>
    Tensor1_Expr<const Tensor3_contracted_01<
                   const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim01>,
                 T, Dim2, i>
    operator()(const Index<j, Dim01> index1, const Index<j, Dim01> index2,
               const Index<i, Dim2> index3)
    {
      using TensorExpr
        = Tensor3_contracted_01<const Christof<T *, Tensor_Dim0, Tensor_Dim12>,
                                T, Dim01>;
      return Tensor1_Expr<TensorExpr, T, Dim2, i>(TensorExpr(*this));
    }

    /* This is for expressions where a number is used for one slot, and
       an index for the others, yielding a Tensor2_symmetric_Expr. */

    template <char i, char j, int N, int Dim12>
    Tensor2_symmetric_Expr<
      const Christof_number_0<const Christof<T *, Tensor_Dim0, Tensor_Dim12>,
                              T, N>,
      T, Dim12, i, j>
    operator()(const Number<N> n1, const Index<i, Dim12> index1,
               const Index<j, Dim12> index2) const
    {
      using TensorExpr
        = Christof_number_0<const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T,
                            N>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim12, i, j>(
        TensorExpr(*this));
    }

    /* An int in one spot, Index for the others, yielding a Tensor2.  I
       can use the same structure for both, since Christof is
       symmetric on the last two indices. */

    template <char i, char j, int Dim0, int Dim2>
    Tensor2_Expr<const Christof_numeral_1<
                   const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T>,
                 T, Dim0, Dim2, i, j>
    operator()(const Index<i, Dim0> index1, const int N,
               const Index<j, Dim2> index2) const
    {
      using TensorExpr
        = Christof_numeral_1<const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(
        TensorExpr(*this, N));
    }

    template <char i, char j, int Dim0, int Dim2>
    Tensor2_Expr<const Christof_numeral_1<
                   const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T>,
                 T, Dim0, Dim2, i, j>
    operator()(const Index<i, Dim0> index1, const Index<j, Dim2> index2,
               const int N) const
    {
      using TensorExpr
        = Christof_numeral_1<const Christof<T *, Tensor_Dim0, Tensor_Dim12>, T>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(
        TensorExpr(*this, N));
    }

    /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

    const Christof<T *, Tensor_Dim0, Tensor_Dim12> &operator++() const
    {
      for(int i = 0; i < Tensor_Dim0; ++i)
        for(int j = 0; j < (Tensor_Dim12 * (Tensor_Dim12 + 1)) / 2; ++j)
          ++data[i][j];
      return *this;
    }
  };
}
