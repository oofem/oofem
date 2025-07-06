/* This file has all of the declarations for expressions like
   Dg*Tensor2_symmetric and Tensor2_symmetric*Dg,
   yielding a Dg, Tensor3, or Tensor1. */

#pragma once

namespace FTensor
{
  /* A(i,j,k)*B(k,l)->Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  class Dg_times_Tensor2_symmetric_0
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Tensor2_symmetric_Expr<B, U, Dim2, k, l> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, Current_Dim - 1) * iterB(Current_Dim - 1, N3)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(0, N3);
    }

  public:
    Dg_times_Tensor2_symmetric_0(
      const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
      const Tensor2_symmetric_Expr<B, U, Dim2, k, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Dg_Expr<Dg_times_Tensor2_symmetric_0<A, B, T, U, Dim01, Dim2, i, j, k, l>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, l>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Tensor2_symmetric_Expr<B, U, Dim2, k, l> &b)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_0<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, l>(
      TensorExpr(a, b));
  }

  /* B(k,l)*A(i,j,k)->Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Dg_Expr<Dg_times_Tensor2_symmetric_0<A, B, T, U, Dim01, Dim2, i, j, k, l>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, l>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim2, k, l> &b,
            const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_0<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, l>(
      TensorExpr(a, b));
  }

  /* A(i,j,k)*B(l,k)->Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  class Dg_times_Tensor2_symmetric_1
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Tensor2_symmetric_Expr<B, U, Dim2, l, k> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, N2, Current_Dim - 1) * iterB(N3, Current_Dim - 1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, N2, 0) * iterB(N3, 0);
    }

  public:
    Dg_times_Tensor2_symmetric_1(
      const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
      const Tensor2_symmetric_Expr<B, U, Dim2, l, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim2>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Dg_Expr<Dg_times_Tensor2_symmetric_1<A, B, T, U, Dim01, Dim2, i, j, k, l>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, l>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Tensor2_symmetric_Expr<B, U, Dim2, l, k> &b)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_1<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, l>(
      TensorExpr(a, b));
  }

  /* B(l,k)*A(i,j,k)->Dg */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Dg_Expr<Dg_times_Tensor2_symmetric_1<A, B, T, U, Dim01, Dim2, i, j, k, l>,
          typename promote<T, U>::V, Dim01, Dim2, i, j, l>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim2, l, k> &b,
            const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_1<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Dg_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2, i, j, l>(
      TensorExpr(a, b));
  }

  /* A(i,j,k)*B(j,l)->Tensor3 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  class Dg_times_Tensor2_symmetric_1_0
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Tensor2_symmetric_Expr<B, U, Dim01, j, l> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1, N2) * iterB(Current_Dim - 1, N3)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, 0, N2) * iterB(0, N3);
    }

  public:
    Dg_times_Tensor2_symmetric_1_0(
      const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
      const Tensor2_symmetric_Expr<B, U, Dim01, j, l> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Tensor3_Expr<
    Dg_times_Tensor2_symmetric_1_0<A, B, T, U, Dim01, Dim2, i, j, k, l>,
    typename promote<T, U>::V, Dim01, Dim2, Dim01, i, k, l>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Tensor2_symmetric_Expr<B, U, Dim01, j, l> &b)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_1_0<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2,
                        Dim01, i, k, l>(TensorExpr(a, b));
  }

  /* B(j,l)*A(i,j,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Tensor3_Expr<
    Dg_times_Tensor2_symmetric_1_0<A, B, T, U, Dim01, Dim2, i, j, k, l>,
    typename promote<T, U>::V, Dim01, Dim2, Dim01, i, k, l>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim01, j, l> &b,
            const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_1_0<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2,
                        Dim01, i, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k)*B(l,j)->Tensor3 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  class Dg_times_Tensor2_symmetric_1_1
  {
    Dg_Expr<A, T, Dim01, Dim2, i, j, k> iterA;
    Tensor2_symmetric_Expr<B, U, Dim01, l, j> iterB;

    template <int Current_Dim>
    typename promote<T, U>::V eval(const int N1, const int N2, const int N3,
                                   const Number<Current_Dim> &) const
    {
      return iterA(N1, Current_Dim - 1, N2) * iterB(N3, Current_Dim - 1)
             + eval(N1, N2, N3, Number<Current_Dim - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const int N2, const int N3, const Number<1> &) const
    {
      return iterA(N1, 0, N2) * iterB(N3, 0);
    }

  public:
    Dg_times_Tensor2_symmetric_1_1(
      const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
      const Tensor2_symmetric_Expr<B, U, Dim01, l, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return eval(N1, N2, N3, Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Tensor3_Expr<
    Dg_times_Tensor2_symmetric_1_1<A, B, T, U, Dim01, Dim2, i, j, k, l>,
    typename promote<T, U>::V, Dim01, Dim2, Dim01, i, k, l>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a,
            const Tensor2_symmetric_Expr<B, U, Dim01, l, j> &b)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_1_1<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2,
                        Dim01, i, k, l>(TensorExpr(a, b));
  }

  /* B(l,j)*A(i,j,k)->Tensor3 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k, char l>
  Tensor3_Expr<
    Dg_times_Tensor2_symmetric_1_1<A, B, T, U, Dim01, Dim2, i, j, k, l>,
    typename promote<T, U>::V, Dim01, Dim2, Dim01, i, k, l>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim01, l, j> &b,
            const Dg_Expr<A, T, Dim01, Dim2, i, j, k> &a)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_1_1<A, B, T, U, Dim01, Dim2, i, j, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim01, Dim2,
                        Dim01, i, k, l>(TensorExpr(a, b));
  }

  /* A(i,j,k)*B(j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Dg_times_Tensor2_symmetric_12
  {
    Dg_Expr<A, T, Dim, Dim, i, j, k> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, j, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(N1, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(N1, 0, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, 0, 0) * iterB(0, 0);
    }

  public:
    Dg_times_Tensor2_symmetric_12(
      const Dg_Expr<A, T, Dim, Dim, i, j, k> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim>(), Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_12<A, B, T, U, Dim, i, j, k>,
               typename promote<T, U>::V, Dim, i>
  operator*(const Dg_Expr<A, T, Dim, Dim, i, j, k> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
  {
    using TensorExpr = Dg_times_Tensor2_symmetric_12<A, B, T, U, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(i,j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_12<A, B, T, U, Dim, i, j, k>,
               typename promote<T, U>::V, Dim, i>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b,
            const Dg_Expr<A, T, Dim, Dim, i, j, k> &a)
  {
    using TensorExpr = Dg_times_Tensor2_symmetric_12<A, B, T, U, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, b));
  }

  /* A(i,j,k)*B(k,j)->Tensor1 */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Dg_times_Tensor2_symmetric_21
  {
    Dg_Expr<A, T, Dim, Dim, i, j, k> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, k, j> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(N1, Current_Dim0 - 1, Current_Dim1 - 1)
               * iterB(Current_Dim1 - 1, Current_Dim0 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(N1, 0, Current_Dim1 - 1) * iterB(Current_Dim1 - 1, 0)
             + eval(N1, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(N1, 0, 0) * iterB(0, 0);
    }

  public:
    Dg_times_Tensor2_symmetric_21(
      const Dg_Expr<A, T, Dim, Dim, i, j, k> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim>(), Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_21<A, B, T, U, Dim, i, j, k>,
               typename promote<T, U>::V, Dim, i>
  operator*(const Dg_Expr<A, T, Dim, Dim, i, j, k> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, k, j> &b)
  {
    using TensorExpr = Dg_times_Tensor2_symmetric_21<A, B, T, U, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, b));
  }

  /* B(k,j)*A(i,j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_21<A, B, T, U, Dim, i, j, k>,
               typename promote<T, U>::V, Dim, i>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, k, j> &b,
            const Dg_Expr<A, T, Dim, Dim, i, j, k> &a)
  {
    using TensorExpr = Dg_times_Tensor2_symmetric_21<A, B, T, U, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, b));
  }

  /* A(j,i,k)*B(j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Dg_times_Tensor2_symmetric_02
  {
    Dg_Expr<A, T, Dim, Dim, j, i, k> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, j, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, N1, Current_Dim1 - 1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(0, N1, Current_Dim1 - 1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(0, N1, 0) * iterB(0, 0);
    }

  public:
    Dg_times_Tensor2_symmetric_02(
      const Dg_Expr<A, T, Dim, Dim, j, i, k> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim>(), Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_02<A, B, T, U, Dim, i, j, k>,
               typename promote<T, U>::V, Dim, i>
  operator*(const Dg_Expr<A, T, Dim, Dim, j, i, k> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
  {
    using TensorExpr = Dg_times_Tensor2_symmetric_02<A, B, T, U, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(j,i,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_02<A, B, T, U, Dim, i, j, k>,
               typename promote<T, U>::V, Dim, i>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b,
            const Dg_Expr<A, T, Dim, Dim, j, i, k> &a)
  {
    using TensorExpr = Dg_times_Tensor2_symmetric_02<A, B, T, U, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, b));
  }

  /* A(k,i,j)*B(j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  class Dg_times_Tensor2_symmetric_20
  {
    Dg_Expr<A, T, Dim, Dim, k, i, j> iterA;
    Tensor2_symmetric_Expr<B, U, Dim, j, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, N1, Current_Dim1 - 1)
               * iterB(Current_Dim1 - 1, Current_Dim0 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(0, N1, Current_Dim1 - 1) * iterB(Current_Dim1 - 1, 0)
             + eval(N1, Number<Dim>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(0, N1, 0) * iterB(0, 0);
    }

  public:
    Dg_times_Tensor2_symmetric_20(
      const Dg_Expr<A, T, Dim, Dim, k, i, j> &a,
      const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim>(), Number<Dim>());
    }
  };

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_20<A, B, T, U, Dim, i, j, k>,
               typename promote<T, U>::V, Dim, i>
  operator*(const Dg_Expr<A, T, Dim, Dim, k, i, j> &a,
            const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b)
  {
    using TensorExpr = Dg_times_Tensor2_symmetric_20<A, B, T, U, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(k,i,j)->Tensor1 */

  template <class A, class B, class T, class U, int Dim, char i, char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_20<A, B, T, U, Dim, i, j, k>,
               typename promote<T, U>::V, Dim, i>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim, j, k> &b,
            const Dg_Expr<A, T, Dim, Dim, k, i, j> &a)
  {
    using TensorExpr = Dg_times_Tensor2_symmetric_20<A, B, T, U, Dim, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim, i>(
      TensorExpr(a, b));
  }

  /* A(j,k,i)*B(j,k)->Tensor1 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  class Dg_times_Tensor2_symmetric_01
  {
    Dg_Expr<A, T, Dim01, Dim2, j, k, i> iterA;
    Tensor2_symmetric_Expr<B, U, Dim01, j, k> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, Current_Dim1 - 1, N1)
               * iterB(Current_Dim0 - 1, Current_Dim1 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(0, Current_Dim1 - 1, N1) * iterB(0, Current_Dim1 - 1)
             + eval(N1, Number<Dim01>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(0, 0, N1) * iterB(0, 0);
    }

  public:
    Dg_times_Tensor2_symmetric_01(
      const Dg_Expr<A, T, Dim01, Dim2, j, k, i> &a,
      const Tensor2_symmetric_Expr<B, U, Dim01, j, k> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim01>(), Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_01<A, B, T, U, Dim01, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim2, i>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, j, k, i> &a,
            const Tensor2_symmetric_Expr<B, U, Dim01, j, k> &b)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_01<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, i>(
      TensorExpr(a, b));
  }

  /* B(j,k)*A(j,k,i)->Tensor1 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_01<A, B, T, U, Dim01, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim2, i>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim01, j, k> &b,
            const Dg_Expr<A, T, Dim01, Dim2, j, k, i> &a)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_01<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, i>(
      TensorExpr(a, b));
  }

  /* A(j,k,i)*B(k,j)->Tensor1 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  class Dg_times_Tensor2_symmetric_10
  {
    Dg_Expr<A, T, Dim01, Dim2, j, k, i> iterA;
    Tensor2_symmetric_Expr<B, U, Dim01, k, j> iterB;

    template <int Current_Dim0, int Current_Dim1>
    typename promote<T, U>::V eval(const int N1, const Number<Current_Dim0> &,
                                   const Number<Current_Dim1> &) const
    {
      return iterA(Current_Dim0 - 1, Current_Dim1 - 1, N1)
               * iterB(Current_Dim1 - 1, Current_Dim0 - 1)
             + eval(N1, Number<Current_Dim0 - 1>(), Number<Current_Dim1>());
    }
    template <int Current_Dim1>
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<Current_Dim1> &) const
    {
      return iterA(0, Current_Dim1 - 1, N1) * iterB(Current_Dim1 - 1, 0)
             + eval(N1, Number<Dim01>(), Number<Current_Dim1 - 1>());
    }
    typename promote<T, U>::V
    eval(const int N1, const Number<1> &, const Number<1> &) const
    {
      return iterA(0, 0, N1) * iterB(0, 0);
    }

  public:
    Dg_times_Tensor2_symmetric_10(
      const Dg_Expr<A, T, Dim01, Dim2, j, k, i> &a,
      const Tensor2_symmetric_Expr<B, U, Dim01, k, j> &b)
        : iterA(a), iterB(b)
    {}
    typename promote<T, U>::V operator()(const int N1) const
    {
      return eval(N1, Number<Dim01>(), Number<Dim01>());
    }
  };

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_10<A, B, T, U, Dim01, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim2, i>
  operator*(const Dg_Expr<A, T, Dim01, Dim2, j, k, i> &a,
            const Tensor2_symmetric_Expr<B, U, Dim01, k, j> &b)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_10<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, i>(
      TensorExpr(a, b));
  }

  /* B(k,j)*A(j,k,i)->Tensor1 */

  template <class A, class B, class T, class U, int Dim01, int Dim2, char i,
            char j, char k>
  Tensor1_Expr<Dg_times_Tensor2_symmetric_10<A, B, T, U, Dim01, Dim2, i, j, k>,
               typename promote<T, U>::V, Dim2, i>
  operator*(const Tensor2_symmetric_Expr<B, U, Dim01, k, j> &b,
            const Dg_Expr<A, T, Dim01, Dim2, j, k, i> &a)
  {
    using TensorExpr
      = Dg_times_Tensor2_symmetric_10<A, B, T, U, Dim01, Dim2, i, j, k>;
    return Tensor1_Expr<TensorExpr, typename promote<T, U>::V, Dim2, i>(
      TensorExpr(a, b));
  }
}
