/* A version for pointers. */

#pragma once

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim12>
  class Tensor3_antisymmetric<T *, Tensor_Dim0, Tensor_Dim12>
  {
    //mutable T *restrict
    //  data[Tensor_Dim0][(Tensor_Dim12 * (Tensor_Dim12 - 1)) / 2];
    mutable T *__restrict
      data[Tensor_Dim0][(Tensor_Dim12 * (Tensor_Dim12 - 1)) / 2];
    
    
  public:
    template <class... U> Tensor3_antisymmetric(U *... d) : data{d...}
    {
      static_assert(sizeof...(d) == sizeof(data) / sizeof(T),
                    "Incorrect number of Arguments. Constructor should "
                    "initialize the entire Tensor");
    }

    Tensor3_antisymmetric() {}

    /* There are two ways of accessing the values inside,
       unsafe(int,int,int) and operator(int,int,int).
       unsafe(int,int,int) will give you a wrong answer if you aren't
       careful.  The problem is that we only store the minimal set of
       components, but some have different signs.  We can't return the
       negative of a component, and assign something to it, because that
       would assign something to a temporary.  To get the correct answer
       if you don't want to change the value, just use
       operator(int,int,int). */

    T &unsafe(const int N1, const int N2, const int N3)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim12 || N2 < 0
         || N3 >= Tensor_Dim12 || N3 < 0 || N2 >= N3)
        {
          std::stringstream s;
          s << "Bad index in Tensor3_antisymmetric<T*," << Tensor_Dim0 << ","
            << Tensor_Dim12 << ">.unsafe(" << N1 << "," << N2 << "," << N3
            << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return *data[N1][N3 - 1 + (N2 * (2 * (Tensor_Dim12 - 1) - N2 - 1)) / 2];
    }

    T operator()(const int N1, const int N2, const int N3) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim12 || N2 < 0
         || N3 >= Tensor_Dim12 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor3_antisymmetric<T*," << Tensor_Dim0 << ","
            << Tensor_Dim12 << ">.operator(" << N1 << "," << N2 << "," << N3
            << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N2 < N3
               ? *data[N1]
                      [N3 - 1 + (N2 * (2 * (Tensor_Dim12 - 1) - N2 - 1)) / 2]
               : (N2 > N3
                    ? -*data[N1][N2 - 1
                                 + (N3 * (2 * (Tensor_Dim12 - 1) - N3 - 1)) / 2]
                    : 0.0);
    }

    T *ptr(const int N1, const int N2, const int N3) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim12 || N2 < 0
         || N3 >= Tensor_Dim12 || N3 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor3_antisymmetric<T*," << Tensor_Dim0 << ","
            << Tensor_Dim12 << ">.ptr(" << N1 << "," << N2 << "," << N3 << ")"
            << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return N2 < N3
               ? data[N1][N3 - 1 + (N2 * (2 * (Tensor_Dim12 - 1) - N2 - 1)) / 2]
               : (N2 > N3
                    ? -data[N1][N2 - 1
                                + (N3 * (2 * (Tensor_Dim12 - 1) - N3 - 1)) / 2]
                    : 0);
    }

    /* These operator()'s are the first part in constructing template
       expressions. */

    template <char i, char j, char k, int Dim0, int Dim12>
    Tensor3_antisymmetric_Expr<
      Tensor3_antisymmetric<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim0, Dim12, i,
      j, k>
    operator()(const Index<i, Dim0> index1, const Index<j, Dim12> index2,
               const Index<k, Dim12> index3)
    {
      return Tensor3_antisymmetric_Expr<
        Tensor3_antisymmetric<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim0, Dim12,
        i, j, k>(*this);
    }

    template <char i, char j, char k, int Dim0, int Dim12>
    Tensor3_antisymmetric_Expr<
      const Tensor3_antisymmetric<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim0,
      Dim12, i, j, k>
    operator()(const Index<i, Dim0> index1, const Index<j, Dim12> index2,
               const Index<k, Dim12> index3) const
    {
      return Tensor3_antisymmetric_Expr<
        const Tensor3_antisymmetric<T *, Tensor_Dim0, Tensor_Dim12>, T, Dim0,
        Dim12, i, j, k>(*this);
    }

    /* The ++ operator increments the pointer, not the number that the
       pointer points to.  This allows iterating over a grid. */

    const Tensor3_antisymmetric<T *, Tensor_Dim0, Tensor_Dim12> &
    operator++() const
    {
      for(int i = 0; i < Tensor_Dim0; ++i)
        for(int j = 0; j < (Tensor_Dim12 * (Tensor_Dim12 - 1)) / 2; ++j)
          ++data[i][j];
      return *this;
    }
  };
}
