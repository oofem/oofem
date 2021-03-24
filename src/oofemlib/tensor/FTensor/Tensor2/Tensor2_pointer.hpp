/* A version for pointers. */

#pragma once

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1>
  class Tensor2<T *, Tensor_Dim0, Tensor_Dim1>
  {
    //mutable T *restrict data[Tensor_Dim0][Tensor_Dim1];
    mutable T *__restrict data[Tensor_Dim0][Tensor_Dim1];

  public:
    /* Initializations for varying numbers of elements. */
    template <class... U> Tensor2(U *... d) : data{d...} {}

    Tensor2() {}

    /* There are two operator(int,int)'s, one for non-consts that lets you
       change the value, and one for consts that doesn't. */

    T &operator()(const int N1, const int N2)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor2<T*," << Tensor_Dim0 << "," << Tensor_Dim1
            << ">.operator(" << N1 << "," << N2 << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return *data[N1][N2];
    }

    T operator()(const int N1, const int N2) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor2<T*," << Tensor_Dim0 << "," << Tensor_Dim1
            << ">.operator(" << N1 << "," << N2 << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return *data[N1][N2];
    }

    T *ptr(const int N1, const int N2) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor2<T*," << Tensor_Dim0 << "," << Tensor_Dim1
            << ">.ptr(" << N1 << "," << N2 << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N1][N2];
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. They are not entirely safe, since you can accidentaly use a
       higher dimension than what is really allowed (like Dim=5). */

    template <char i, char j, int Dim0, int Dim1>
    Tensor2_Expr<Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, Dim0, Dim1, i, j>
    operator()(const Index<i, Dim0>, const Index<j, Dim1> index2)
    {
      return Tensor2_Expr<Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, Dim0,
                          Dim1, i, j>(*this);
    }

    template <char i, char j, int Dim0, int Dim1>
    Tensor2_Expr<const Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, Dim0, Dim1,
                 i, j>
    operator()(const Index<i, Dim0>, const Index<j, Dim1> index2) const
    {
      return Tensor2_Expr<const Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T,
                          Dim0, Dim1, i, j>(*this);
    }

    /* This is for expressions where a number is used for one slot, and
       an index for another, yielding a Tensor1_Expr.  The non-const
       versions don't actually create a Tensor2_number_rhs_[01] object.
       They create a Tensor1_Expr directly, which provides the
       appropriate indexing operators.  The const versions do create a
       Tensor2_number_[01]. */

    template <char i, int Dim, int N>
    Tensor1_Expr<
      Tensor2_number_rhs_1<Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, N>, T,
      Dim, i>
    operator()(const Index<i, Dim>, const Number<N>)
    {
      using TensorExpr
        = Tensor2_number_rhs_1<Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int Dim, int N>
    Tensor1_Expr<
      Tensor2_number_rhs_0<Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, N>, T,
      Dim, i>
    operator()(const Number<N>, const Index<i, Dim>)
    {
      using TensorExpr
        = Tensor2_number_rhs_0<Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int Dim, int N>
    Tensor1_Expr<const Tensor2_number_1<
                   const Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, N>,
                 T, Dim, i>
    operator()(const Index<i, Dim>, const Number<N>) const
    {
      using TensorExpr
        = Tensor2_number_1<const Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    template <char i, int Dim, int N>
    Tensor1_Expr<const Tensor2_number_0<
                   const Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, N>,
                 T, Dim, i>
    operator()(const Number<N>, const Index<i, Dim>) const
    {
      using TensorExpr
        = Tensor2_number_0<const Tensor2<T *, Tensor_Dim0, Tensor_Dim1>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

    const Tensor2<T *, Tensor_Dim0, Tensor_Dim1> &operator++() const
    {
      for(int i = 0; i < Tensor_Dim0; ++i)
        for(int j = 0; j < Tensor_Dim1; ++j)
          ++data[i][j];
      return *this;
    }

    /* These two operator()'s return the Tensor2 with internal
       contractions, yielding a T.  I have to specify one for both
       const and non-const because otherwise they compiler will use the
       operator() which gives a Tensor2_Expr<>. */

    template <char i, int Dim>
    T operator()(const Index<i, Dim>, const Index<i, Dim> index2)
    {
      return internal_contract(Number<Dim>());
    }

    template <char i, int Dim>
    T operator()(const Index<i, Dim>, const Index<i, Dim> index2) const
    {
      return internal_contract(Number<Dim>());
    }

  private:
    template <int N> T internal_contract(Number<N>) const
    {
      return *data[N - 1][N - 1] + internal_contract(Number<N - 1>());
    }

    T internal_contract(Number<1>) const { return *data[0][0]; }
  };
}
