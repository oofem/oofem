#pragma once
#include <nanobind/ndarray.h>
#include <iostream>
#include "floatmatrix.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "intarray.h"

NAMESPACE_BEGIN(NB_NAMESPACE)
NAMESPACE_BEGIN(detail)

namespace nb=::nanobind;

template<typename T, typename T2=void> struct _ArrayInfo{};


struct _ArrayInfo1d{
    typedef nb::shape<-1> Shape;
    typedef nb::f_contig Contig;
    template<typename T>
    static void resize(T& v, size_t shape0, size_t shape1){ v.resize(shape0); }
    template<typename T>
    static void shape_strides_size(const T& v, size_t shape[2], int64_t strides[2], size_t& size){
        shape[0]=v.size(); shape[1]=0; strides[0]=1; strides[1]=0; size=v.size();
    };

};
struct _ArrayInfo2d{
    typedef nb::shape<-1,-1> Shape;
    typedef nb::f_contig Contig;
    template<typename T>
    static void resize(T& v, size_t shape0, size_t shape1){ v.resize(shape0,shape1); }
    template<typename T>
    static void shape_strides_size(const T& v, size_t shape[2], int64_t strides[2], size_t& size){
        shape[0]=v.giveNumberOfRows(); shape[1]=v.giveNumberOfColumns(); strides[0]=1; strides[1]=v.giveNumberOfRows(); size=v.giveNumberOfRows()*v.giveNumberOfColumns();
    };
};
template<> struct _ArrayInfo<::oofem::FloatMatrix>: public _ArrayInfo2d{};
template<> struct _ArrayInfo<::oofem::FloatArray>: public _ArrayInfo1d{};
template<> struct _ArrayInfo<::oofem::IntArray>: public _ArrayInfo1d{};

template<typename T>
struct type_caster<
    T,std::enable_if_t<
        std::is_same<T,::oofem::FloatMatrix>::value || std::is_same<T,::oofem::FloatArray>::value || std::is_same<T,::oofem::IntArray>::value
       ,int
    >
>
{
    using Scalar = typename T::Scalar;
    using Info = _ArrayInfo<T>;
    using NDArray = nanobind::ndarray<Scalar,nanobind::numpy,typename Info::Shape,typename Info::Contig>;
    using NDArrayConst = nanobind::ndarray<Scalar,nanobind::numpy,typename Info::Shape,typename Info::Contig>;
    using NDArrayCaster = make_caster<NDArray>;
    using NDArrayCasterConst = make_caster<NDArrayConst>;

    NB_TYPE_CASTER(T, NDArrayCaster::Name)

    bool from_python(handle src, uint8_t flags, cleanup_list *cleanup) noexcept {
        // We're in any case making a copy, so non-writable inputs area also okay
        // std::cerr<<"@@ from_python["<<<<std::endl;
        // std::cerr<<"@@    ["<<nb::str(nb::steal(src).attr("__class__")).c_str()<<"]"<<std::endl;
        NDArrayCasterConst caster;
        if (!caster.from_python(src, flags & ~(uint8_t)cast_flags::accepts_none, cleanup)){
            // std::cerr<<"@@   → ERROR"<<std::endl;
            // std::cerr<<"from_python["<<typeid(T).name()<<"]: failure; the object passed must be a numpy.array."<<std::endl;
            return false;
        }
        const NDArrayConst &array = caster.value;
        Info::resize(value,array.shape(0),array.shape(1));
        // The layout is contiguous & compatible thanks to array_for_eigen_t<T>
        memcpy(value.givePointer(), array.data(), array.size() * sizeof(Scalar));
        return true;
    }

    template <typename T2>
    static handle from_cpp(T2 &&v, rv_policy policy, cleanup_list *cleanup) noexcept {
        // std::cerr<<"@@ from_cpp["<<typeid(v).name()<<"]"<<std::endl;
        policy = infer_policy<T2>(policy);
        if constexpr (std::is_pointer_v<T2>)
            return from_cpp_internal((const T &) *v, policy, cleanup);
        else
            return from_cpp_internal((const T &) v, policy, cleanup);
    }

    static handle from_cpp_internal(const T &v, rv_policy policy, cleanup_list *cleanup) noexcept {
        size_t shape[2];
        int64_t strides[2];
        size_t size;
        Info::template shape_strides_size<T>(v,shape,strides,size);
        // std::cerr<<"    from_cpp_internal[shape=["<<shape[0]<<","<<shape[1]<<"],strides=["<<strides[0]<<","<<strides[1]<<"],size="<<size<<"]"<<std::endl;
        void *ptr = (void *) v.givePointer();
        // for(int i=0; i<size; i++) std::cerr<<"["<<i<<"] → "<<*(double*)(((void*)v.givePointer())+i*sizeof(double))<<std::endl;
        if (policy == rv_policy::move) {
            // Don't bother moving when the data is static or occupies <1KB
            if (size < (1024 / sizeof(Scalar))) policy = rv_policy::copy;
        }

        object owner;
        if (policy == rv_policy::move) {
            T *temp = new T(std::move(v));
            owner = capsule(temp, [](void *p) noexcept { delete (T *) p; });
            ptr = temp->givePointer();
            policy = rv_policy::reference;
        } else if (policy == rv_policy::reference_internal && cleanup->self()) {
            owner = borrow(cleanup->self());
            policy = rv_policy::reference;
        }

        object o = steal(NDArrayCaster::from_cpp(
            NDArray(ptr, T::Dim, shape, owner, strides), policy, cleanup)
        );

        return o.release();
    }
};


NAMESPACE_END(detail)
NAMESPACE_END(NB_NAMESPACE)
