/* A general version, not for pointers. */

#pragma once

namespace FTensor
{
  template <class T, int Tensor_Dim01, int Tensor_Dim23> class Ddg
  {
    T data[(Tensor_Dim01 * (Tensor_Dim01 + 1)) / 2]
          [(Tensor_Dim23 * (Tensor_Dim23 + 1)) / 2];

  public:
    /* There are two operator(int,int,int,int)'s, one for non-consts
       that lets you change the value, and one for consts that
       doesn't. */

    T &operator()(const int N1, const int N2, const int N3, const int N4)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim01 || N1 < 0 || N2 >= Tensor_Dim01 || N2 < 0
         || N3 >= Tensor_Dim23 || N3 < 0 || N4 >= Tensor_Dim23 || N4 < 0)
        {
          std::stringstream s;
          s << "Bad index in Dg<T," << Tensor_Dim01 << "," << Tensor_Dim23
            << ">.operator(" << N1 << "," << N2 << "," << N3 << "," << N4
            << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 > N2
               ? (N3 > N4 ? data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2]
                                [N3 + (N4 * (2 * Tensor_Dim23 - N4 - 1)) / 2]
                          : data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2]
                                [N4 + (N3 * (2 * Tensor_Dim23 - N3 - 1)) / 2])
               : (N3 > N4 ? data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2]
                                [N3 + (N4 * (2 * Tensor_Dim23 - N4 - 1)) / 2]
                          : data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2]
                                [N4 + (N3 * (2 * Tensor_Dim23 - N3 - 1)) / 2]);
    }

    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim01 || N1 < 0 || N2 >= Tensor_Dim01 || N2 < 0
         || N3 >= Tensor_Dim23 || N3 < 0 || N4 >= Tensor_Dim23 || N4 < 0)
        {
          std::stringstream s;
          s << "Bad index in Dg<T," << Tensor_Dim01 << "," << Tensor_Dim23
            << ">.operator(" << N1 << "," << N2 << "," << N3 << "," << N4
            << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N1 > N2
               ? (N3 > N4 ? data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2]
                                [N3 + (N4 * (2 * Tensor_Dim23 - N4 - 1)) / 2]
                          : data[N1 + (N2 * (2 * Tensor_Dim01 - N2 - 1)) / 2]
                                [N4 + (N3 * (2 * Tensor_Dim23 - N3 - 1)) / 2])
               : (N3 > N4 ? data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2]
                                [N3 + (N4 * (2 * Tensor_Dim23 - N4 - 1)) / 2]
                          : data[N2 + (N1 * (2 * Tensor_Dim01 - N1 - 1)) / 2]
                                [N4 + (N3 * (2 * Tensor_Dim23 - N3 - 1)) / 2]);
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. They are not entirely safe, since you can accidently use a
       higher dimension than what is really allowed (like Dim=5). */

    template <char i, char j, char k, char l, int Dim01, int Dim23>
    typename std::enable_if<(Tensor_Dim01 >= Dim01 && Tensor_Dim23 >= Dim23),
                            Ddg_Expr<Ddg<T, Tensor_Dim01, Tensor_Dim23>, T,
                                     Dim01, Dim23, i, j, k, l>>::type
    operator()(const Index<i, Dim01>, const Index<j, Dim01>,
               const Index<k, Dim23>, const Index<l, Dim23>)
    {
      return Ddg_Expr<Ddg<T, Tensor_Dim01, Tensor_Dim23>, T, Dim01, Dim23, i,
                      j, k, l>(*this);
    }

    template <char i, char j, char k, char l, int Dim01, int Dim23>
    typename std::enable_if<(Tensor_Dim01 >= Dim01 && Tensor_Dim23 >= Dim23),
                            Ddg_Expr<const Ddg<T, Tensor_Dim01, Tensor_Dim23>,
                                     T, Dim01, Dim23, i, j, k, l>>::type
    operator()(const Index<i, Dim01>, const Index<j, Dim01>,
               const Index<k, Dim23>, const Index<l, Dim23>) const
    {
      return Ddg_Expr<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T, Dim01,
                      Dim23, i, j, k, l>(*this);
    }

    /* This is for expressions where a number is used for two slots, and
       an index for the other two, yielding a Tensor2_symmetric_Expr. */

    template <char i, char j, int N0, int N1, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 > N0 && Tensor_Dim01 > N1 && Tensor_Dim23 >= Dim),
      Tensor2_symmetric_Expr<
        Ddg_number_01<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T, N0, N1>, T,
        Dim, i, j>>::type
    operator()(const Number<N0>, const Number<N1>, const Index<i, Dim>,
               const Index<j, Dim>) const
    {
      using TensorExpr
        = Ddg_number_01<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T, N0, N1>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this));
    }

    template <char i, char j, int N0, int N1, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 > N0 && Tensor_Dim01 > N1 && Tensor_Dim23 >= Dim),
      Tensor2_symmetric_Expr<
        Ddg_number_rhs_01<Ddg<T, Tensor_Dim01, Tensor_Dim23>, T, N0, N1>, T,
        Dim, i, j>>::type
    operator()(const Number<N0>, const Number<N1>, const Index<i, Dim>,
               const Index<j, Dim>)
    {
      using TensorExpr
        = Ddg_number_rhs_01<Ddg<T, Tensor_Dim01, Tensor_Dim23>, T, N0, N1>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(*this);
    }

    /* This is for expressions where a number is used for one slot, and
       an index for the other three, yielding a Dg_Expr. */

    template <char i, char j, char k, int N0, int Dim1, int Dim23>
    typename std::enable_if<
      (Tensor_Dim01 > N0 && Tensor_Dim01 >= Dim1 && Tensor_Dim23 >= Dim23),
      Dg_Expr<Ddg_number_0<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T, N0>, T,
              Dim23, Dim1, i, j, k>>::type
    operator()(const Number<N0>, const Index<k, Dim1>, const Index<i, Dim23>,
               const Index<j, Dim23>) const
    {
      using TensorExpr
        = Ddg_number_0<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T, N0>;
      return Dg_Expr<TensorExpr, T, Dim23, Dim1, i, j, k>(TensorExpr(*this));
    }

    template <char i, char j, char k, int N0, int Dim1, int Dim23>
    typename std::enable_if<
      (Tensor_Dim01 > N0 && Tensor_Dim01 >= Dim1 && Tensor_Dim23 >= Dim23),
      Dg_Expr<Ddg_number_rhs_0<Ddg<T, Tensor_Dim01, Tensor_Dim23>, T, N0>, T,
              Dim23, Dim1, i, j, k>>::type
    operator()(const Number<N0>, const Index<k, Dim1>, const Index<i, Dim23>,
               const Index<j, Dim23>)
    {
      using TensorExpr
        = Ddg_number_rhs_0<Ddg<T, Tensor_Dim01, Tensor_Dim23>, T, N0>;
      return Dg_Expr<TensorExpr, T, Dim23, Dim1, i, j, k>(*this);
    }

    /* This is for expressions where an int (not a Number) is used for
       two slots, and an index for the other two, yielding a
       Tensor2_symmetric_Expr. */

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim23 >= Dim),
      Tensor2_symmetric_Expr<
        Ddg_numeral_01<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>, T, Dim, i,
        j>>::type
    operator()(const int N0, const int N1, const Index<i, Dim>,
               const Index<j, Dim>) const
    {
      using TensorExpr
        = Ddg_numeral_01<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this, N0, N1));
    }

    template <char i, char j, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim),
      Tensor2_symmetric_Expr<
        Ddg_numeral_23<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>, T, Dim, i,
        j>>::type
    operator()(const Index<i, Dim>, const Index<j, Dim>, const int N2,
               const int N3) const
    {
      using TensorExpr
        = Ddg_numeral_23<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>;
      return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
        TensorExpr(*this, N2, N3));
    }

    /* int in two slots but yielding a Tensor2 */

    template <char i, char j, int Dim1, int Dim3>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim1 && Tensor_Dim23 >= Dim3),
      Tensor2_Expr<Ddg_numeral_02<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>,
                   T, Dim1, Dim3, i, j>>::type
    operator()(const int N0, const Index<i, Dim1>, const int N2,
               const Index<j, Dim3>) const
    {
      using TensorExpr
        = Ddg_numeral_02<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim3, i, j>(
        TensorExpr(*this, N0, N2));
    }

    /* int in three slots, Index in the other yielding a Tensor1_Expr. */

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim),
      Tensor1_Expr<Ddg_numeral_123<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>,
                   T, Dim, i>>::type
    operator()(const Index<i, Dim>, const int N1, const int N2, const int N3)
    {
      using TensorExpr
        = Ddg_numeral_123<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(
        TensorExpr(*this, N1, N2, N3));
    }

    template <char i, int Dim>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim),
      Tensor1_Expr<Ddg_numeral_123<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>,
                   T, Dim, i>>::type
    operator()(const int N1, const Index<i, Dim>, const int N2, const int N3)
    {
      using TensorExpr
        = Ddg_numeral_123<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>;
      return Tensor1_Expr<TensorExpr, T, Dim, i>(
        TensorExpr(*this, N1, N2, N3));
    }

    /* This is for expressions where an int (not a Number) is used for
       one slot, and an index for the other three, yielding a
       Dg_Expr. */

    template <char i, char j, char k, int Dim1, int Dim23>
    typename std::enable_if<
      (Tensor_Dim01 >= Dim1 && Tensor_Dim23 >= Dim23),
      Dg_Expr<Ddg_numeral_0<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>, T,
              Dim23, Dim1, i, j, k>>::type
    operator()(const int N0, const Index<k, Dim1>, const Index<i, Dim23>,
               const Index<j, Dim23>) const
    {
      using TensorExpr
        = Ddg_numeral_0<const Ddg<T, Tensor_Dim01, Tensor_Dim23>, T>;
      return Dg_Expr<TensorExpr, T, Dim23, Dim1, i, j, k>(
        TensorExpr(*this, N0));
    }
  };
}
