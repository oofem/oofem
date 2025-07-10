/* Takes a derivative of a Tensor1, yielding a Tensor2. */

//  template<class T, int Dim, char i, char j>
//  class dTensor1
//  {
//    const Tensor1<T*,Dim> &a;
//    const int di,dj,dk;
//    const double dx,dy,dz;
//  public:
//    typename promote<T,double>::V operator()(const int N1, const int N2) const
//    {
//      return N2==0 ? (*(a.ptr(N1)+di)-*(a.ptr(N1)-di))*dx*0.5
//        : (N2==1 ? (*(a.ptr(N1)+dj)-*(a.ptr(N1)-dj))*dy*0.5
//  	 : (*(a.ptr(N1)+dk)-*(a.ptr(N1)-dk))*dz*0.5);
//    }
//    dTensor1(const Tensor1<T*,Dim> &A,
//  	   const int Di, const int Dj, const int Dk,
//  	   const double Dx, const double Dy, const double Dz):
//      a(A), di(Di), dj(Dj), dk(Dk), dx(Dx), dy(Dy), dz(Dz) {}
//  };

//  template<class T, int Dim, char i, char j>
//  const Tensor2_Expr<const dTensor1<T,Dim,i,j>,
//    typename promote<T,double>::V,Dim,3,i,j>
//  d(const Tensor1<T*,Dim> &a, const Index<i,Dim> index1, const Index<j,3>
//  index2,
//    const int &di, const int &dj, const int &dk,
//    const double &dx, const double &dy, const double &dz)
//  {
//    typedef dTensor1<T,Dim,i,j> TensorExpr;
//    return Tensor2_Expr<TensorExpr,typename promote<T,double>::V,Dim,3,i,j>
//      (TensorExpr(a,di,dj,dk,dx,dy,dz));
//  }

#pragma once

namespace FTensor
{
  template <class T, int Dim0, int Dim1, char i, char j> class dTensor1
  {
    const Tensor1<T *, Dim0> &a;
    const Tensor1<int, Dim1> &d_ijk;
    const Tensor1<double, Dim1> &d_xyz;

  public:
    typename promote<T, double>::V operator()(const int N1, const int N2) const
    {
      return (*(a.ptr(N1) + d_ijk(N2)) - *(a.ptr(N1) - d_ijk(N2))) * d_xyz(N2)
             * 0.5;
    }
    dTensor1(const Tensor1<T *, Dim0> &A, const Tensor1<int, Dim1> &D_ijk,
             const Tensor1<double, Dim1> &D_xyz)
        : a(A), d_ijk(D_ijk), d_xyz(D_xyz)
    {}
  };

  template <class T, int Dim0, int Dim1, char i, char j>
  const Tensor2_Expr<const dTensor1<T, Dim0, Dim1, i, j>,
                     typename promote<T, double>::V, Dim0, Dim1, i, j>
  d(const Tensor1<T *, Dim0> &a, const Index<i, Dim0> index1,
    const Index<j, Dim1> index2, const Tensor1<int, Dim1> &d_ijk,
    const Tensor1<double, Dim1> &d_xyz)
  {
    using TensorExpr = dTensor1<T, Dim0, Dim1, i, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, double>::V, Dim0, Dim1,
                        i, j>(TensorExpr(a, d_ijk, d_xyz));
  }
}
