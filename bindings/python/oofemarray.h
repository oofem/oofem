NAMESPACE_BEGIN(NB_NAMESPACE)
NAMESPACE_BEGIN(detail)

namespace nb=::nanobind;

template<typename T, typename T2=void> struct _ArrayInfo{};


struct _ArrayInfo1d{
    static constexpr int Dim=1;
    typedef nb::shape<-1> Shape;
    typedef nb::c_contig Contig;
    template<typename T>
    static void resize(T& v, size_t shape0, size_t shape1){ v.resize(shape0); }
    template<typename T, typename Scalar>
    static void shape_strides_size(const T& v, size_t shape[2], int64_t strides[2], size_t& size){
        shape[0]=v.size(); shape[1]=0; strides[0]=sizeof(Scalar); strides[1]=0; size=v.size();
    };

};
struct _ArrayInfo2d{
    static constexpr int Dim=2;
    typedef nb::shape<-1,-1> Shape;
    typedef nb::f_contig Contig;
    template<typename T>
    static void resize(T& v, size_t shape0, size_t shape1){ v.resize(shape0,shape1); }
    template<typename T, typename Scalar>
    static void shape_strides_size(const T& v, size_t shape[2], int64_t strides[2], size_t& size){
        shape[0]=v.giveRowSize(); shape[1]=v.giveColSize(); strides[0]=sizeof(Scalar); strides[1]=sizeof(Scalar)*v.giveRowSize(); size=v.giveRowSize()*v.giveColSize();
    };
};
template<> struct _ArrayInfo<FloatMatrix>: public _ArrayInfo2d{
    typedef double Scalar;
    typedef FloatMatrix T;
};

template<> struct _ArrayInfo<FloatArray>: public _ArrayInfo1d{
    typedef double Scalar;
    typedef FloatArray T;
};
template<> struct _ArrayInfo<IntArray>: public _ArrayInfo1d{
    typedef int Scalar;
    typedef IntArray T;
};

template<>
struct type_caster<FloatMatrix> {
    using T = FloatMatrix;
    using Info = _ArrayInfo<T>;
    using NDArray = nanobind::ndarray<Info::Scalar,nanobind::numpy,Info::Shape,Info::Contig>;
    using NDArrayConst = nanobind::ndarray<const Info::Scalar,nanobind::numpy,Info::Shape,Info::Contig>;
    using NDArrayCaster = make_caster<NDArray>;
    using NDArrayCasterConst = make_caster<NDArrayConst>;

    NB_TYPE_CASTER(T, NDArrayCaster::Name)

    bool from_python(handle src, uint8_t flags, cleanup_list *cleanup) noexcept {
        // We're in any case making a copy, so non-writable inputs area also okay
        //using NDArrayConst = array_for_eigen_t<T, const typename T::Scalar>;
        //make_caster<NDArrayConst> caster;
        NDArrayCasterConst caster;
        if (!caster.from_python(src, flags & ~(uint8_t)cast_flags::accepts_none, cleanup))
            return false;
        const NDArrayConst &array = caster.value;
        Info::resize(value,array.shape(0),array.shape(1));
        // The layout is contiguous & compatible thanks to array_for_eigen_t<T>
        memcpy(value.givePointer(), array.data(), array.size() * sizeof(Info::Scalar));
        return true;
    }

    template <typename T2>
    static handle from_cpp(T2 &&v, rv_policy policy, cleanup_list *cleanup) noexcept {
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
        Info::shape_strides_size<T,Info::Scalar>(v,shape,strides,size);
        void *ptr = (void *) v.givePointer();
        if (policy == rv_policy::move) {
            // Don't bother moving when the data is static or occupies <1KB
            if (size < (1024 / sizeof(Info::Scalar))) policy = rv_policy::copy;
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
            NDArray(ptr, Info::Dim, shape, owner, strides), policy, cleanup)
        );

        return o.release();
    }
};


NAMESPACE_END(detail)
NAMESPACE_END(NB_NAMESPACE)
