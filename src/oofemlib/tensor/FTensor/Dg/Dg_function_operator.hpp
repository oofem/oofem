/* This contains the definitions of the almost all of the indexing
   operators for Dg. */

#pragma once

namespace FTensor
{
  /* These operator()'s are the first part in constructing template
     expressions. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, char k, int Dim01, int Dim2>
  Dg_Expr<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim01> index1, const Index<j, Dim01> index2,
             const Index<k, Dim2> index3)
  {
    return Dg_Expr<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k>(
      *this);
  }

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, char k, int Dim01, int Dim2>
  Dg_Expr<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i, j, k>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim01> index1, const Index<j, Dim01> index2,
             const Index<k, Dim2> index3) const
  {
    return Dg_Expr<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01, Dim2, i,
                   j, k>(*this);
  }

  /* These operators are for internal contractions. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int Dim, int Dim12>
  Tensor1_Expr<
    Tensor3_contracted_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim12, i>,
    T, Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim> index1, const Index<j, Dim12> index2,
             const Index<j, Dim12> index3) const
  {
    using TensorExpr
      = Tensor3_contracted_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim12,
                              i>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
  }

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int Dim, int Dim02>
  Tensor1_Expr<
    Tensor3_contracted_02<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim02, i>,
    T, Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<j, Dim02> index1, const Index<i, Dim> index2,
             const Index<j, Dim02> index3) const
  {
    using TensorExpr
      = Tensor3_contracted_02<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim02,
                              i>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
  }

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int Dim, int Dim01>
  Tensor1_Expr<
    Tensor3_contracted_01<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01, i>,
    T, Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<j, Dim01> index1, const Index<j, Dim01> index2,
             const Index<i, Dim> index3) const
  {
    using TensorExpr
      = Tensor3_contracted_01<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, Dim01,
                              i>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
  }

  /* This is for expressions where a number is used for one slot, and
     indices for the others, yielding a Tensor2_Expr or
     Tensor2_symmetric_Expr.  The non-const versions don't actually
     create a Dg_number_rhs_* object, while the const versions
     do create a Dg_number_*. */

  /* First slot. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int N, int Dim0, int Dim1>
  Tensor2_Expr<Dg_number_rhs_0<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>, T,
               Dim0, Dim1, i, j>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Number<N> n1, const Index<i, Dim0> index1,
             const Index<j, Dim1> index2)
  {
    using TensorExpr = Dg_number_rhs_0<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
    return Tensor2_Expr<TensorExpr, T, Dim0, Dim1, i, j>(*this);
  }

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int N, int Dim0, int Dim1>
  Tensor2_Expr<Dg_number_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>, T,
               Dim0, Dim1, i, j>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Number<N> n1, const Index<i, Dim0> index1,
             const Index<j, Dim1> index2) const
  {
    using TensorExpr
      = Dg_number_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
    return Tensor2_Expr<TensorExpr, T, Dim0, Dim1, i, j>(TensorExpr(*this));
  }

  /* Second slot. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int N, int Dim0, int Dim1>
  Tensor2_Expr<Dg_number_rhs_0<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>, T,
               Dim0, Dim1, i, j>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim0> index1, const Number<N> n1,
             const Index<j, Dim1> index2)
  {
    using TensorExpr = Dg_number_rhs_0<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
    return Tensor2_Expr<TensorExpr, T, Dim0, Dim1, i, j>(*this);
  }

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int N, int Dim0, int Dim1>
  Tensor2_Expr<Dg_number_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>, T,
               Dim0, Dim1, i, j>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim0> index1, const Number<N> n1,
             const Index<j, Dim1> index2) const
  {
    using TensorExpr
      = Dg_number_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
    return Tensor2_Expr<TensorExpr, T, Dim0, Dim1, i, j>(TensorExpr(*this));
  }

  /* Third slot. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int N, int Dim>
  Tensor2_symmetric_Expr<
    Dg_number_rhs_2<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>, T, Dim, i, j>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim> index1, const Index<j, Dim> index2,
             const Number<N> n1)
  {
    using TensorExpr = Dg_number_rhs_2<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
    return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(*this);
  }

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int N, int Dim>
  Tensor2_symmetric_Expr<
    Dg_number_2<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>, T, Dim, i, j>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim> index1, const Index<j, Dim> index2,
             const Number<N> n1) const
  {
    using TensorExpr
      = Dg_number_2<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N>;
    return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(TensorExpr(*this));
  }

  /* This is for expressions where a number is used for two slots, and
     an Index for the other, yielding a Tensor1_Expr.  The non-const
     versions don't actually create a Dg_number_rhs_* object,
     while the const versions do create a Dg_number_*. */

  /* Index in first slot. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, int N1, int N2, int Dim>
  Tensor1_Expr<Dg_number_rhs_12<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>,
               T, Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim> index, const Number<N1> n1,
             const Number<N2> n2)
  {
    using TensorExpr
      = Dg_number_rhs_12<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
  }

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, int N1, int N2, int Dim>
  Tensor1_Expr<Dg_number_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>,
               T, Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim> index, const Number<N1> n1,
             const Number<N2> n2) const
  {
    using TensorExpr
      = Dg_number_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
  }

  /* Index in second slot.  I use the same structures as for the Index
     in the first slot since the tensor is symmetric on the first two
     indices. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, int N1, int N2, int Dim>
  Tensor1_Expr<Dg_number_rhs_12<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>,
               T, Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Number<N1> n1, const Index<i, Dim> index,
             const Number<N2> n2)
  {
    using TensorExpr
      = Dg_number_rhs_12<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
  }

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, int N1, int N2, int Dim>
  Tensor1_Expr<Dg_number_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>,
               T, Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Number<N1> n1, const Index<i, Dim> index,
             const Number<N2> n2) const
  {
    using TensorExpr
      = Dg_number_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
  }

  /* Index in third slot. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, int N1, int N2, int Dim>
  Tensor1_Expr<Dg_number_rhs_01<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>,
               T, Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Number<N1> n1, const Number<N2> n2,
             const Index<i, Dim> index)
  {
    using TensorExpr
      = Dg_number_rhs_01<Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(*this);
  }

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, int N1, int N2, int Dim>
  Tensor1_Expr<Dg_number_01<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>,
               T, Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Number<N1> n1, const Number<N2> n2,
             const Index<i, Dim> index) const
  {
    using TensorExpr
      = Dg_number_01<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T, N1, N2>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this));
  }

  /* Specializations for using actual numbers instead of Number<>.
     This is for expressions where an actual number is used for one
     slot, and indices for the others, yielding a Tensor2_Expr or
     Tensor2_symmetric_Expr. I only define the const versions. */

  /* First slot. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int Dim0, int Dim1>
  Tensor2_Expr<Dg_numeral_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T,
               Dim0, Dim1, i, j>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const int N, const Index<i, Dim0> index1,
             const Index<j, Dim1> index2) const
  {
    using TensorExpr = Dg_numeral_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
    return Tensor2_Expr<TensorExpr, T, Dim0, Dim1, i, j>(TensorExpr(*this, N));
  }

  /* Second slot. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int Dim0, int Dim1>
  Tensor2_Expr<Dg_numeral_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T,
               Dim0, Dim1, i, j>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim0> index1, const int N,
             const Index<j, Dim1> index2) const
  {
    using TensorExpr = Dg_numeral_0<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
    return Tensor2_Expr<TensorExpr, T, Dim0, Dim1, i, j>(TensorExpr(*this, N));
  }

  /* Third slot. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, char j, int Dim>
  Tensor2_symmetric_Expr<
    Dg_numeral_2<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T, Dim, i, j>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim> index1, const Index<j, Dim> index2,
             const int N) const
  {
    using TensorExpr = Dg_numeral_2<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
    return Tensor2_symmetric_Expr<TensorExpr, T, Dim, i, j>(
      TensorExpr(*this, N));
  }

  /* This is for expressions where an actual number is used for two
     slots, and an Index for the other, yielding a Tensor1_Expr. I
     only define the const versions. */

  /* Index in first slot. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, int Dim>
  Tensor1_Expr<Dg_numeral_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T,
               Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const Index<i, Dim> index, const int N1, const int N2) const
  {
    using TensorExpr
      = Dg_numeral_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
  }

  /* Index in second slot.  I use the same structures as for the Index
     in the first slot since the tensor is symmetric on the first two
     indices. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, int Dim>
  Tensor1_Expr<Dg_numeral_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T,
               Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const int N1, const Index<i, Dim> index, const int N2) const
  {
    using TensorExpr
      = Dg_numeral_12<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
  }

  /* Index in third slot. */

  template <class T, int Tensor_Dim01, int Tensor_Dim2>
  template <char i, int Dim>
  Tensor1_Expr<Dg_numeral_01<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>, T,
               Dim, i>
  Dg<T, Tensor_Dim01, Tensor_Dim2>::
  operator()(const int N1, const int N2, const Index<i, Dim> index) const
  {
    using TensorExpr
      = Dg_numeral_01<const Dg<T, Tensor_Dim01, Tensor_Dim2>, T>;
    return Tensor1_Expr<TensorExpr, T, Dim, i>(TensorExpr(*this, N1, N2));
  }
}
