/* Declares a wrapper class for antisymmetric rank 2 Tensor expressions.
   It is specialized for Tensor3_number_rhs_2.  */

#pragma once

namespace FTensor
{
  template <class A, class T, int Dim, char i, char j>
  class Tensor2_antisymmetric_Expr
  {
    A iter;

  public:
    Tensor2_antisymmetric_Expr(const A &a) : iter(a) {}
    T operator()(const int N1, const int N2) const { return iter(N1, N2); }
  };

  template <class A, class T, int Tensor_Dim, int Dim, char i, char j>
  class Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A, Tensor_Dim>, T,
                                   Dim, i, j>
  {
    Tensor2_antisymmetric<A, Tensor_Dim> &iter;

  public:
    Tensor2_antisymmetric_Expr(Tensor2_antisymmetric<A, Tensor_Dim> &a)
        : iter(a)
    {}
    T &operator()(const int N1, const int N2) { return iter(N1, N2); }
    T operator()(const int N1, const int N2) const { return iter(N1, N2); }

    // /* Various assignment operators.  I have to explicitly declare the
    //    second operator= because otherwise the compiler will generate its
    //    own and not use the template code. */

    // template<class B, class U>
    // const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator=(const Tensor2_antisymmetric_Expr<B,U,Dim,i,j> &result);

    // const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator=(const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // &result);

    // template<class B, class U>
    // const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator+=(const Tensor2_antisymmetric_Expr<B,U,Dim,i,j> &result);

    // template<class B, class U>
    // const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator-=(const Tensor2_antisymmetric_Expr<B,U,Dim,i,j> &result);

    // template<class B, class U>
    // const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator&=(const Tensor2_antisymmetric_Expr<B,U,Dim,i,j> &result);

    // /* This is for when the indices are switched (i,j) -> (j,i). */

    // template<class B, class U>
    // const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator=(const Tensor2_antisymmetric_Expr<B,U,Dim,j,i> &result);

    // template<class B, class U>
    // const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator+=(const Tensor2_antisymmetric_Expr<B,U,Dim,j,i> &result);

    // template<class B, class U>
    // const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator-=(const Tensor2_antisymmetric_Expr<B,U,Dim,j,i> &result);

    // /* Operations with just generics. */

    // template<class U>
    // const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator=(const U &d); template<class U>  const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator+=(const U &d); template<class U>  const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator-=(const U &d); template<class U>  const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator*=(const U &d); template<class U>  const
    // Tensor2_antisymmetric_Expr<Tensor2_antisymmetric<A,Tensor_Dim>,T,Dim,i,j>
    // & operator/=(const U &d);
  };
}
