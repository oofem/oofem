/* A version for pointers. */

#pragma once

namespace FTensor
{
  template <class T, int Tensor_Dim> class Tensor2_symmetric<T *, Tensor_Dim>
  {
    //mutable T *restrict data[(Tensor_Dim * (Tensor_Dim + 1)) / 2];
    mutable T *__restrict data[(Tensor_Dim * (Tensor_Dim + 1)) / 2];

  public:
    template <class... U> Tensor2_symmetric(U *... d) : data{d...}
    {
      static_assert(sizeof...(d) == sizeof(data) / sizeof(T),
                    "Incorrect number of Arguments. Constructor should "
                    "initialize the entire Tensor");
    }

    Tensor2_symmetric() {}

    /* There are two operator(int,int)'s, one for non-consts that lets you
       change the value, and one for consts that doesn't. */

    T &operator()(const int N1, const int N2)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim || N1 < 0 || N2 >= Tensor_Dim || N2 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor2_symmetric<T*," << Tensor_Dim
            << ">.operator(" << N1 << "," << N2 << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 > N2 ? *data[N1 + (N2 * (2 * Tensor_Dim - N2 - 1)) / 2]
                     : *data[N2 + (N1 * (2 * Tensor_Dim - N1 - 1)) / 2];
    }

    T operator()(const int N1, const int N2) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim || N1 < 0 || N2 >= Tensor_Dim || N2 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor2_symmetric<T*," << Tensor_Dim
            << ">.operator(" << N1 << "," << N2 << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 > N2 ? *data[N1 + (N2 * (2 * Tensor_Dim - N2 - 1)) / 2]
                     : *data[N2 + (N1 * (2 * Tensor_Dim - N1 - 1)) / 2];
    }

    T *ptr(const int N1, const int N2) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim || N1 < 0 || N2 >= Tensor_Dim || N2 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor2_symmetric<T*," << Tensor_Dim << ">.ptr("
            << N1 << "," << N2 << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 > N2 ? data[N1 + (N2 * (2 * Tensor_Dim - N2 - 1)) / 2]
                     : data[N2 + (N1 * (2 * Tensor_Dim - N1 - 1)) / 2];
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. They are not entirely safe, since you can accidentaly use a
       higher dimension than what is really allowed (like Dim=5). */

    /* This returns a Tensor2_Expr, since the indices are not really
       symmetric anymore since they cover different dimensions. */

    template <char i, char j, int Dim0, int Dim1>
    Tensor2_Expr<Tensor2_symmetric<T *, Tensor_Dim>, T, Dim0, Dim1, i, j>
    operator()(const Index<i, Dim0> index1, const Index<j, Dim1> index2)
    {
      return Tensor2_Expr<Tensor2_symmetric<T *, Tensor_Dim>, T, Dim0, Dim1, i,
                          j>(*this);
    }

    template <char i, char j, int Dim0, int Dim1>
    Tensor2_Expr<const Tensor2_symmetric<T *, Tensor_Dim>, T, Dim0, Dim1, i, j>
    operator()(const Index<i, Dim0> index1, const Index<j, Dim1> index2) const
    {
      return Tensor2_Expr<const Tensor2_symmetric<T *, Tensor_Dim>, T, Dim0,
                          Dim1, i, j>(*this);
    }

    /* This returns a Tensor2_symmetric_Expr, since the indices are still
       symmetric on the lower dimensions. */

    template <char i, char j, int Dim>
    Tensor2_symmetric_Expr<Tensor2_symmetric<T *, Tensor_Dim>, T, Dim, i, j>
    operator()(const Index<i, Dim> index1, const Index<j, Dim> index2)
    {
      return Tensor2_symmetric_Expr<Tensor2_symmetric<T *, Tensor_Dim>, T, Dim,
                                    i, j>(*this);
    }

    template <char i, char j, int Dim>
    Tensor2_symmetric_Expr<const Tensor2_symmetric<T *, Tensor_Dim>, T, Dim, i,
                           j>
    operator()(const Index<i, Dim> index1, const Index<j, Dim> index2) const
    {
      return Tensor2_symmetric_Expr<const Tensor2_symmetric<T *, Tensor_Dim>,
                                    T, Dim, i, j>(*this);
    }

    /* This is for expressions where a number is used for one slot, and
       an index for another, yielding a Tensor1_Expr.  The non-const
       versions don't actually create a Tensor2_number_rhs_[01] object.
       They create a Tensor1_Expr directly, which provides the
       appropriate indexing operators.  The const versions do create a
       Tensor2_number_[01]. */

    template <char i, int N, int Dim>
    Tensor1_Expr<Tensor2_number_rhs_1<Tensor2_symmetric<T *, Tensor_Dim>, T, N>,
                 T, Dim, i>
    operator()(const Index<i, Dim> index1, const Number<N> &n1)
    {
      using TensorExpr
        = Tensor2_number_rhs_1<Tensor2_symmetric<T *, Tensor_Dim>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N, int Dim>
    Tensor1_Expr<
      const Tensor2_number_1<const Tensor2_symmetric<T *, Tensor_Dim>, T, N>,
      T, Dim, i>
    operator()(const Index<i, Dim> index1, const Number<N> &n1) const
    {
      using TensorExpr
        = Tensor2_number_1<const Tensor2_symmetric<T *, Tensor_Dim>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    template <char i, int N, int Dim>
    Tensor1_Expr<Tensor2_number_rhs_0<Tensor2_symmetric<T *, Tensor_Dim>, T, N>,
                 T, Dim, i>
    operator()(const Number<N> &n1, const Index<i, Dim> index1)
    {
      using TensorExpr
        = Tensor2_number_rhs_0<Tensor2_symmetric<T *, Tensor_Dim>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
    }

    template <char i, int N, int Dim>
    Tensor1_Expr<
      const Tensor2_number_0<const Tensor2_symmetric<T *, Tensor_Dim>, T, N>,
      T, Dim, i>
    operator()(const Number<N> &n1, const Index<i, Dim> index1) const
    {
      using TensorExpr
        = Tensor2_number_0<const Tensor2_symmetric<T *, Tensor_Dim>, T, N>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
    }

    /* Specializations for using actual numbers instead of Number<> */

    template <char i, int Dim>
    Tensor1_Expr<
      const Tensor2_numeral_1<const Tensor2_symmetric<T *, Tensor_Dim>, T>, T,
      Dim, i>
    operator()(const Index<i, Dim> index1, const int N) const
    {
      using TensorExpr
        = Tensor2_numeral_1<const Tensor2_symmetric<T *, Tensor_Dim>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N));
    }

    template <char i, int Dim>
    Tensor1_Expr<
      const Tensor2_numeral_0<const Tensor2_symmetric<T *, Tensor_Dim>, T>, T,
      Dim, i>
    operator()(const int N, const Index<i, Dim> index1) const
    {
      using TensorExpr
        = Tensor2_numeral_0<const Tensor2_symmetric<T *, Tensor_Dim>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N));
    }

    /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

    const Tensor2_symmetric<T *, Tensor_Dim> &operator++() const
    {
      for(int i = 0; i < (Tensor_Dim * (Tensor_Dim + 1)) / 2; ++i)
        ++data[i];
      return *this;
    }

    /* These two operator()'s return the Tensor2 with internal
       contractions, yielding a T.  I have to specify one for both
       const and non-const because otherwise the compiler will use the
       operator() which gives a Tensor2_Expr<>. */

    template <char i, int Dim>
    T operator()(const Index<i, Dim> index1, const Index<i, Dim> index2)
    {
      return internal_contract(Number<Dim>());
    }

    template <char i, int Dim>
    T operator()(const Index<i, Dim> index1, const Index<i, Dim> index2) const
    {
      return internal_contract(Number<Dim>());
    }

  private:
    template <int N> T internal_contract(Number<N>) const
    {
      return *data[N - 1 + ((N - 1) * (2 * Tensor_Dim - N)) / 2]
             + internal_contract(Number<N - 1>());
    }

    T internal_contract(Number<1>) const { return *data[0]; }
  };
}
