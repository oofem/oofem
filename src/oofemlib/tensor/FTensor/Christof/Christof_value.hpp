/* A general version, not for pointers */

#pragma once

#include "../Dg.hpp"

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim12> class Christof
  {
    T data[Tensor_Dim0][(Tensor_Dim12 * (Tensor_Dim12 + 1)) / 2];

  public:
    template <class... U> Christof(U... d) : data{d...}
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
          s << "Bad index in Christof<T," << Tensor_Dim0 << "," << Tensor_Dim12
            << ">.operator(" << N1 << "," << N2 << "," << N3 << ")"
            << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N2 > N3 ? data[N1][N2 + (N3 * (2 * Tensor_Dim12 - N3 - 1)) / 2]
                     : data[N1][N3 + (N2 * (2 * Tensor_Dim12 - N2 - 1)) / 2];
    }

    T operator()(const int N1, const int N2, const int N3) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim12 || N2 < 0
         || N3 >= Tensor_Dim12 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Christof<T," << Tensor_Dim0 << "," << Tensor_Dim12
            << ">.operator(" << N1 << "," << N2 << "," << N3 << ") const"
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
    typename std::enable_if<(Tensor_Dim0 >= Dim0 && Tensor_Dim12 >= Dim12),
                            Dg_Expr<Christof<T, Tensor_Dim0, Tensor_Dim12>, T,
                                    Dim12, Dim0, i, j, k>>::type
    operator()(const Index<k, Dim0>, const Index<i, Dim12>,
               const Index<j, Dim12>)
    {
      return Dg_Expr<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim12, Dim0, i,
                     j, k>(*this);
    }

    template <char i, char j, char k, int Dim0, int Dim12>
    typename std::enable_if<(Tensor_Dim0 >= Dim0 && Tensor_Dim12 >= Dim12),
                            Dg_Expr<const Christof<T, Tensor_Dim0, Tensor_Dim12>,
                                    T, Dim12, Dim0, i, j, k>>::type
    operator()(const Index<k, Dim0>, const Index<i, Dim12>,
               const Index<j, Dim12>) const
    {
      return Dg_Expr<const Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim12,
                     Dim0, i, j, k>(*this);
    }

    /* This is for expressions where a number is used for two slots, and
       an Index for the other, yielding a Tensor1_Expr.  The non-const
       versions don't actually create a Dg_number_rhs_* object,
       while the const versions do create a Dg_number_*. */

    /* Index in first slot. */

    template <char i, int N1, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim12 > N1 && Tensor_Dim12 > N2),
      Tensor1_Expr<
        Dg_number_rhs_12<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, N1, N2>, T,
        Dim, i>>::type
    operator()(const Index<i, Dim>, const Number<N1>, const Number<N2>)
    {
      using TensorExpr
        = Dg_number_rhs_12<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, N1, N2>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N1, int N2, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim12 > N1 && Tensor_Dim12 > N2),
      Tensor1_Expr<
        Dg_number_12<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, N1, N2>, T,
        Dim, i>>::type
    operator()(const Index<i, Dim>, const Number<N1>, const Number<N2>) const
    {
      using TensorExpr
        = Dg_number_12<const Christof<T, Tensor_Dim0, Tensor_Dim12>, T, N1, N2>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }
    // TODO index on second and third position

    /* These operators are for internal contractions.  Some of them are
       less general, because, for example, A(j,i,j) with i and j having
       different dimensions is ambiguous. */

    /* const versions */

    template <char i, char j, int Dim0, int Dim12>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim12 >= Dim12),
      Tensor1_Expr<
        Tensor3_contracted_12<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim12>,
        T, Dim0, i>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim12>,
               const Index<j, Dim12>) const
    {
      using TensorExpr
        = Tensor3_contracted_12<Christof<T, Tensor_Dim0, Tensor_Dim12>, T,
                                Dim12>;
      return Tensor1_Expr<TensorExpr, T, Dim0, i>(TensorExpr(*this));
    }

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim12 >= Dim),
      Tensor1_Expr<
        Tensor3_contracted_02<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim>,
        T, Dim, i>>::type
    operator()(const Index<j, Dim>, const Index<i, Dim>,
               const Index<j, Dim>) const
    {
      using TensorExpr
        = Tensor3_contracted_02<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }
    // TODO allow different dimensions here ^^^^ with the one below
    //    template<char i, char j, int Dim02, int Dim1>
    //    Tensor1_Expr<const Tensor3_contracted_02
    //    <const Christof<T,Tensor_Dim0,Tensor_Dim12>,T,Dim02>,T,Dim1,i>
    //    operator()(const Index<j,Dim02> index1, const Index<i,Dim1> index2,
    //  	     const Index<j,Dim02> index3) const
    //    {
    //      typedef Tensor3_contracted_02
    //        <const Christof<T,Tensor_Dim0,Tensor_Dim12>,T,Dim02>
    //        TensorExpr;
    //      return Tensor1_Expr<TensorExpr,T,Dim1,i>(TensorExpr(*this));
    //    }

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim12 >= Dim),
      Tensor1_Expr<
        Tensor3_contracted_01<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim>,
        T, Dim, i>>::type
    operator()(const Index<j, Dim>, const Index<j, Dim>,
               const Index<i, Dim>) const
    {
      using TensorExpr
        = Tensor3_contracted_01<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }
    // TODO allow different dimensions here ^^^^ with the one below
    //    template<char i, char j, int Dim01, int Dim2>
    //    Tensor1_Expr<const Tensor3_contracted_01
    //    <const Christof<T,Tensor_Dim0,Tensor_Dim12>,T,Dim01>,T,Dim2,i>
    //    operator()(const Index<j,Dim01> index1, const Index<j,Dim01> index2,
    //  	     const Index<i,Dim2> index3) const
    //    {
    //      typedef Tensor3_contracted_01
    //        <const Christof<T,Tensor_Dim0,Tensor_Dim12>,T,Dim01>
    //        TensorExpr;
    //      return Tensor1_Expr<TensorExpr,T,Dim2,i>(TensorExpr(*this));
    //    }

    /* non-const versions, needed so that overload resolution doesn't
       pick the more general indexing operator. */

    template <char i, char j, int Dim0, int Dim12>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim12 >= Dim12),
      Tensor1_Expr<
        Tensor3_contracted_12<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim12>,
        T, Dim0, i>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim12>,
               const Index<j, Dim12>)
    {
      using TensorExpr
        = Tensor3_contracted_12<Christof<T, Tensor_Dim0, Tensor_Dim12>, T,
                                Dim12>;
      return Tensor1_Expr<TensorExpr, T, Dim0, i>(TensorExpr(*this));
    }

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim12 >= Dim),
      Tensor1_Expr<
        Tensor3_contracted_02<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim>,
        T, Dim, i>>::type
    operator()(const Index<j, Dim>, const Index<i, Dim>, const Index<j, Dim>)
    {
      using TensorExpr
        = Tensor3_contracted_02<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }
    // TODO allow different dimensions here ^^^^ with the one below
    //    template<char i, char j, int Dim02, int Dim1>
    //    Tensor1_Expr<const Tensor3_contracted_02
    //    <const Christof<T,Tensor_Dim0,Tensor_Dim12>,T,Dim02>,T,Dim1,i>
    //    operator()(const Index<j,Dim02> index1, const Index<i,Dim1> index2,
    //  	     const Index<j,Dim02> index3)
    //    {
    //      typedef Tensor3_contracted_02
    //        <const Christof<T,Tensor_Dim0,Tensor_Dim12>,T,Dim02>
    //        TensorExpr;
    //      return Tensor1_Expr<TensorExpr,T,Dim1,i>(TensorExpr(*this));
    //    }

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim && Tensor_Dim12 >= Dim),
      Tensor1_Expr<
        Tensor3_contracted_01<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim>,
        T, Dim, i>>::type
    operator()(const Index<j, Dim>, const Index<j, Dim>, const Index<i, Dim>)
    {
      using TensorExpr
        = Tensor3_contracted_01<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, Dim>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    //    template<char i, char j, int Dim01, int Dim2>
    //    Tensor1_Expr<const Tensor3_contracted_01
    //    <const Christof<T,Tensor_Dim0,Tensor_Dim12>,T,Dim01>,T,Dim2,i>
    //    operator()(const Index<j,Dim01> index1, const Index<j,Dim01> index2,
    //  	     const Index<i,Dim2> index3)
    //    {
    //      typedef Tensor3_contracted_01
    //        <const Christof<T,Tensor_Dim0,Tensor_Dim12>,T,Dim01>
    //        TensorExpr;
    //      return Tensor1_Expr<TensorExpr,T,Dim2,i>(TensorExpr(*this));
    //    }

    /* This is for expressions where a Number<> is used for one slot, and
       an index for the others, yielding a Tensor2_symmetric_Expr. */

    template <char i, char j, int N, int Dim12>
    typename std::enable_if<
      (Tensor_Dim0 > N && Tensor_Dim12 >= Dim12),
      Tensor2_symmetric_Expr<
        Christof_number_0<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, N>, T,
        Dim12, i, j>>::type
    operator()(const Number<N>, const Index<i, Dim12>,
               const Index<j, Dim12>) const
    {
      using TensorExpr
        = Christof_number_0<Christof<T, Tensor_Dim0, Tensor_Dim12>, T, N>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim12, i, j>(
        TensorExpr(*this));
    }
    // TODO Allow multiple dimensions here
    /* This is for expressions where an int is used for one slot, and
       an index for the others, yielding a Tensor2_symmetric_Expr. */

    template <char i, char j, int Dim12>
    typename std::enable_if<
      (Tensor_Dim12 >= Dim12),
      Tensor2_symmetric_Expr<
        Christof_numeral_0<Christof<T, Tensor_Dim0, Tensor_Dim12>, T>, T,
        Dim12, i, j>>::type
    operator()(const int N, const Index<i, Dim12>, const Index<j, Dim12>) const
    {
      using TensorExpr
        = Christof_numeral_0<Christof<T, Tensor_Dim0, Tensor_Dim12>, T>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim12, i, j>(
        TensorExpr(*this, N));
    }
    // TODO allow multiple dimensions here
    /* An int in one spot, Index for the others, yielding a Tensor2.  I
       can use the same structure for both, since Christof is
       symmetric on the last two indices. */

    template <char i, char j, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim12 >= Dim2),
      Tensor2_Expr<Christof_numeral_1<Christof<T, Tensor_Dim0, Tensor_Dim12>, T>,
                   T, Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const int N, const Index<j, Dim2>) const
    {
      using TensorExpr
        = Christof_numeral_1<Christof<T, Tensor_Dim0, Tensor_Dim12>, T>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(
        TensorExpr(*this, N));
    }

    template <char i, char j, int Dim0, int Dim2>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim12 >= Dim2),
      Tensor2_Expr<Christof_numeral_1<Christof<T, Tensor_Dim0, Tensor_Dim12>, T>,
                   T, Dim0, Dim2, i, j>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim2>, const int N) const
    {
      using TensorExpr
        = Christof_numeral_1<Christof<T, Tensor_Dim0, Tensor_Dim12>, T>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, j>(
        TensorExpr(*this, N));
    }
  };
}
