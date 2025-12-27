/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <boost/python.hpp>
//#include <boost/python/return_internal_reference.hpp>
//#include <boost/python/manage_new_object.hpp>
//#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <boost/python/raw_function.hpp>
using namespace std;
using namespace boost::python;
namespace bp = boost::python;
#include <memory> // for std::auto_ptr<>



/*****************************************************
*
* O O F E M L I B   M O D U L E
*
*****************************************************/
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "engngm.h"
#include "unknownnumberingscheme.h"
#include "domain.h"
#include "sparsemtrx.h"
#include "load.h"
#include "initialcondition.h"
#include "function.h"
#include "material.h"
#include "crosssection.h"
#include "node.h"
#include "connectivitytable.h"
#include "spatiallocalizer.h"
#include "errorestimator.h"
#include "nodalrecoverymodel.h"
#include "oofemtxtinputrecord.h"
#include "oofemtxtdatareader.h"
#include "logger.h"
#include "util.h"
#include "problemmode.h"
#include "dofmanager.h"
#include "elementgeometrytype.h"
#include "field.h"
#include "uniformgridfield.h"
#include "unstructuredgridfield.h"
#include "datastream.h"
#include "dofmanvalfield.h"
#include "dofmantransftype.h"
#include "fieldmanager.h"
#include "generalboundarycondition.h"
#include "boundarycondition.h"
#include "Loads/foreigntempfieldload.h"
#include "integrationrule.h"
#include "gausspoint.h"
#include "internalstatetype.h"
#include "matresponsemode.h"
#include "sm/Materials/structuralmaterial.h"
#include "matstatus.h"
#include "sm/Materials/structuralms.h"
#include "exportmodulemanager.h"
#include "outputmanager.h"
#include "classfactory.h"

#include "sm/Materials/structmatsettable.h"

#include<iostream>

namespace oofem {

/*****************************************************
*
* O O F E M L I B   C L A S S E S
*
*****************************************************/


// not using auto_ptr avoids many warnings, but potentially causes crashes?
// #define OOFEM_USE_DEPRECATED_AUTOPTR


/**********************************
*********** CONVERTERS ************
**********************************/

/*

Converters serve to automatically convert python types to corresponding c++ types when calling
c++ functions from python with arguments. Many converters are provided by boost::python (such as
number conversions, str to std::string and such. The ones we define below add extra converters
specific to this project.

The FloatArray converters (from-python only) allow one to pass python sequence (list or tuple) where
FloatArray is expected, which results in a more natural syntax in scripts.

std::vector converters allow to pass sequence where std::vector<T> is required (to-python; provided
that all sequence's items are convertible to T themselves) and, on the other hand, allow to return
std::vector<T> from c++ functions which will be converted to list automatically, without having to
write a wrapper function which would just do this conversion.

The templates below don't define the converters yet, they are instantiated in the module initialization
routine at the end. If converters are missing for some type (such as std::vector<bool>), they can be
added at the end by simply instantiating these templates with the respective argument.
*/

// copied from https://github.com/eudoxos/minieigen/blob/master/src/converters.hpp (written by myself)
template<typename T>
bool pySeqItemCheck(PyObject* o, int i){ return bp::extract<T>(bp::object(bp::handle<>(PySequence_GetItem(o,i)))).check(); }
template<typename T>
T pySeqItemExtract(PyObject* o, int i){ return bp::extract<T>(bp::object(bp::handle<>(PySequence_GetItem(o,i))))(); }

// AnyArray will be FloatArray or IntArray, Scalar is double or int
template<typename AnyArray, typename Scalar>
struct custom_AnyArray_from_sequence
{
    custom_AnyArray_from_sequence(){ bp::converter::registry::push_back(&convertible,&construct,bp::type_id<AnyArray>()); }
    static void* convertible(PyObject* obj_ptr){
        // dont check size, since FloatArray/IntArray are always dynamic-sized
        if(!PySequence_Check(obj_ptr)) return 0;
        // check that sequence items are convertible to scalars (should be done in other converters as well?!); otherwise Matrix3 is convertible to Vector3, but then we fail in *construct* very unclearly (TypeError: No registered converter was able to produce a C++ rvalue of type double from this Python object of type Vector3)
        size_t len=PySequence_Size(obj_ptr);
        for(size_t i=0; i<len; i++) if(!pySeqItemCheck<Scalar>(obj_ptr,i)) return 0;
        return obj_ptr;
    }
    static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data){
        void* storage=((bp::converter::rvalue_from_python_storage<AnyArray>*)(data))->storage.bytes;
        new (storage) AnyArray;
        size_t len;
        len=PySequence_Size(obj_ptr);
        ((AnyArray*)storage)->resize(len);
        for(size_t i=0; i<len; i++) (*((AnyArray*)storage))[i]=pySeqItemExtract<Scalar>(obj_ptr,i);
        data->convertible=storage;
    }
};

// copied from https://github.com/woodem/woo/blob/master/lib/pyutil/converters.hpp (written by myself)
/*** python sequence to std::vector<T> ***/
template<typename containedType>
struct custom_vector_from_seq
{
    custom_vector_from_seq(){ bp::converter::registry::push_back(&convertible,&construct,bp::type_id<vector<containedType> >()); }
    static void* convertible(PyObject* obj_ptr){
        // the second condition is important, for some reason otherwise there were attempted conversions of Body to list which failed afterwards.
        if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
        return obj_ptr;
    }
    static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data){
        void* storage=((bp::converter::rvalue_from_python_storage<vector<containedType> >*)(data))->storage.bytes;
        new (storage) vector<containedType>();
        vector<containedType>* v=(vector<containedType>*)(storage);
        int l=PySequence_Size(obj_ptr); if(l<0) abort(); /*std::cerr<<"l="<<l<<"; "<<typeid(containedType).name()<<std::endl;*/ v->reserve(l); for(int i=0; i<l; i++) { v->push_back(bp::extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
        data->convertible=storage;
    }
};
/*** std::vector<T> to python list ***/
template<typename containedType>
struct custom_vector_to_list
{
    static PyObject* convert(const vector<containedType>& v) {
        bp::list ret; for(const containedType& e: v) ret.append(e);
        return bp::incref(ret.ptr());
    }
};

/** shorthand for registering both converter above **/
template<typename T>
void converters_cxxVector_pyList_2way()
{
    custom_vector_from_seq<T>(); bp::to_python_converter<vector<T>,custom_vector_to_list<T>>();
};



/*****************************************************
* FloatArray
*****************************************************/


void (FloatArray::*floatarray_add_1)(const FloatArray&) = &FloatArray::add;
void (FloatArray::*floatarray_add_2)(double, const FloatArray&) = &FloatArray::add;
void (FloatArray::*floatarray_add_3)(double) = &FloatArray::add;
double (FloatArray::*dotProduct_1)(const FloatArray&) const = &FloatArray::dotProduct;
double (FloatArray::*dotProduct_2)(const FloatArray&, int) const = &FloatArray::dotProduct;
void (FloatArray::*beDifferenceOf_1)(const FloatArray&, const FloatArray&) = &FloatArray::beDifferenceOf;
void (FloatArray::*beDifferenceOf_2)(const FloatArray&, const FloatArray&, int) = &FloatArray::beDifferenceOf;
void (FloatArray::*floatarray_printYourself_1)() const  = &FloatArray::printYourself;

std::string FloatArray_str(const FloatArray& a){
    std::string ret("[");
    for(int i=0; i<a.giveSize(); i++){ ret+=(i>0?", ":"")+std::to_string(a[i]); }
    return ret+"]";
}

void pyclass_FloatArray()
{
    class_<FloatArray>("FloatArray") // should not have boost::noncopyable so that it can be returned from c++ funcs
        .def(init< optional<int> >())
        .def(init< FloatArray& >())
        .def("resize", &FloatArray::resize, "Checks size of receiver towards requested bounds. If dimension mismatch, size is adjusted accordingly")
        .def("containsOnlyZeroes", &FloatArray::containsOnlyZeroes, "Returns nonzero if all coefficients of the receiver are 0, else returns zero")
        .def("giveSize", &FloatArray::giveSize, "Returns the size of receiver")
        .def("isNotEmpty", &FloatArray::isNotEmpty, "Returns true if receiver is not empty")
        .def("isEmpty", &FloatArray::isEmpty, "Returns true if receiver is empty")
        .def("negated", &FloatArray::negated, "Switches the sign of every coefficient of receiver")
        .def("printYourself", floatarray_printYourself_1, "Print receiver on stdout")
        .def("zero", &FloatArray::zero, "Zeroes all coefficients of receiver")
        .def("beProductOf", &FloatArray::beProductOf, "Receiver becomes the result of the product of aMatrix and anArray. Adjusts the size of receiver if necessary")
        .def("beTProductOf", &FloatArray::beTProductOf, "Receiver becomes the result of the product of aMatrix^T and anArray. Adjusts the size of receiver if necessary")
        .def("beVectorProductOf", &FloatArray::beVectorProductOf, "Computes vector product (or cross product) of vectors given as parameters")
        .def("add", floatarray_add_1, "Adds array src to receiver. If the receiver's size is zero, it adjusts its size to size of src array. If recever's size is nonzero and different from src array size an error is generated")
        .def("add", floatarray_add_2, "Adds array times factor to receiver. If the receiver's size is zero, it adjusts its size to size of b array. If receiver's size is nonzero and different from b")
        .def("add", floatarray_add_3, "Adds scalar to receiver")
        .def("subtract", &FloatArray::subtract, "Subtracts array src to receiver. If the receiver's size is zero, it adjusts its size to size of src array. If recever's size is nonzero and different from src")
        .def("times", &FloatArray::times, "Multiplies receiver with scalar")
        .def("beMaxOf", &FloatArray::beMaxOf, "Sets receiver to maximum of a or b's respective elements")
        .def("beMinOf", &FloatArray::beMinOf, "Sets receiver to be minimum of a or b's respective elements")
        .def("beDifferenceOf", beDifferenceOf_1, "Sets receiver to be a - b")
        .def("beDifferenceOf", beDifferenceOf_2, "Sets receiver to be a - b, using only the first n entries")
        .def("assemble", &FloatArray::assemble, "Assembles the array fe (typically, the load vector of a finite element) into the receiver, using loc as location array")
        .def("dotProduct", dotProduct_1, "Computes the dot product (or inner product) of receiver and argument")
        .def("dotProduct", dotProduct_2, "Computes the dot product (or inner product) of receiver and argument using first 'size' elements")
        .def("normalize", &FloatArray::normalize, "Normalizes receiver. Euclidean norm is used, after operation receiver will have this norm equal to 1.0")
        .def("computeNorm", &FloatArray::computeNorm, "Computes the norm (or length) of the vector")
        .def("computeSquaredNorm", &FloatArray::computeSquaredNorm, "Computes the square of the norm")
        .def("sum", &FloatArray::sum, "Computes the sum of receiver values")

        .def("__str__", &FloatArray_str,"Return printable representation.")

        .def("__len__", &FloatArray::giveSize, "Returns the size of receiver")
        .def("__getitem__", &FloatArray::__getitem__, "Coefficient access function. Provides 0-based indexing access")
        .def("__setitem__", &FloatArray::__setitem__, "Coefficient access function. Provides 0-based indexing access")
        .def("beCopyOf", &FloatArray::beCopyOf, "Modifies receiver to become copy of given parameter")
        ;
}


/*****************************************************
* FloatMatrix
*****************************************************/
bool (FloatMatrix::*solveForRhs_1)(const FloatArray &b, FloatArray &answer, bool) = &FloatMatrix::solveForRhs;
void (FloatMatrix::*solveForRhs_2)(const FloatMatrix &b, FloatMatrix &answer, bool) = &FloatMatrix::solveForRhs;
void (FloatMatrix::*assemble_1)(const FloatMatrix&, const IntArray&) = &FloatMatrix::assemble;
void (FloatMatrix::*assemble_2)(const FloatMatrix&, const IntArray&, const IntArray&) = &FloatMatrix::assemble;
void (FloatMatrix::*assemble_3)(const FloatMatrix&, const int*, const int*) = &FloatMatrix::assemble;
void (FloatMatrix::*floatmatrix_add_1)(const FloatMatrix&) = &FloatMatrix::add;
void (FloatMatrix::*floatmatrix_add_2)(double, const FloatMatrix&) = &FloatMatrix::add;
void (FloatMatrix::*floatmatrix_printYourself_1)() const = &FloatMatrix::printYourself;


void pyclass_FloatMatrix()
{
    class_<FloatMatrix, boost::noncopyable>("FloatMatrix")
        .def(init<int, int>())
        .def(init< FloatMatrix& >())
        .def("assemble", assemble_1, "Assembles the contribution using localization array into receiver. The receiver must have dimensions large enough to localize contribution")
        .def("assemble", assemble_2, "Assembles the contribution using localization array into receiver. The receiver must have dimensions large enough to localize contribution")
        .def("assemble", assemble_3, "Assembles the contribution using localization array into receiver. The receiver must have dimensions large enough to localize contribution")
        .def("assembleu", assemble_2, "Assembles the contribution using localization array into receiver. The receiver must have dimensions large enough to localize contribution")
        .def("giveDeterminant", &FloatMatrix::giveDeterminant, "Returns determinant of the receiver. Receiver should be square matrix")
        .def("zero", &FloatMatrix::zero, "Zeroes all coefficient of receiver")
        .def("beUnitMatrix", &FloatMatrix::beUnitMatrix, "Sets receiver to unity matrix")
        .def("beTranspositionOf", &FloatMatrix::beTranspositionOf, "Assigns to the receiver the transposition of parameter")
        .def("beProductOf", &FloatMatrix::beProductOf, "Assigns to the receiver product of a * b")
        .def("beTProductOf", &FloatMatrix::beTProductOf, "Assigns to the receiver product of a^T * b$")
        .def("beProductTOf", &FloatMatrix::beProductTOf, "Assigns to the receiver product of a * b^T$")
        .def("beDyadicProductOf", &FloatMatrix::beDyadicProductOf, "Assigns to the receiver the dyadic product v_1 * v_2^T$")
        .def("beInverseOf", &FloatMatrix::beInverseOf, "Modifies receiver to become inverse of given parameter. Size of receiver will be adjusted")
        .def("solveForRhs", solveForRhs_1, "Solves the  system of linear equations K * a = b. Uses Gaussian elimination with pivoting directly on receiver")
        .def("solveForRhs", solveForRhs_2, "Solves the  system of linear equations K * A = B. Uses Gaussian elimination with pivoting directly on receiver")
        .def("plusProductSymmUpper", &FloatMatrix::plusProductSymmUpper, "Adds to the receiver the product a^T * b dV. If the receiver has zero size, it is expanded. Assumes that receiver and product  a^T * b dV are symmetric matrices. Computes only the upper half of receiver")
        .def("plusProductUnsym", &FloatMatrix::plusProductUnsym, "Adds to the receiver the product a^T * b dV. If the receiver has zero size, it is expanded")
        .def("add", floatmatrix_add_1, "Adds matrix to the receiver. If receiver has zero size, size is accordingly adjusted")
        .def("add", floatmatrix_add_2, "Adds matrix to the receiver. If receiver has zero size, size is accordingly adjusted")
        .def("subtract", &FloatMatrix::subtract, "Subtracts matrix from the receiver")
        .def("times", &FloatMatrix::times, "Multiplies receiver by factor f")
        .def("negated", &FloatMatrix::negated, "Changes sign of receiver values")
        .def("symmetrized", &FloatMatrix::symmetrized, "Initializes the lower half of the receiver according to the upper half")
        .def("rotatedWith", &FloatMatrix::rotatedWith, "Returns the receiver 'a' transformed using give transformation matrix r. The method performs the operation  a = r^T * a * r")
        .def("resize", &FloatMatrix::resize, "Checks size of receiver towards requested bounds. If dimension mismatch, size is adjusted accordingly")
        .def("printYourself", floatmatrix_printYourself_1, "Prints matrix to stdout")

        .def("__setitem__", &FloatMatrix::__setitem__, "Coefficient access function. Implements 0-based indexing")
        .def("__getitem__", &FloatMatrix::__getitem__, "Coefficient access function. Implements 0-based indexing")
        .def("beCopyOf", &FloatMatrix::beCopyOf, "Modifies receiver to become copy of given parameter")
        ;
}


/*****************************************************
* IntArray
*****************************************************/
void (IntArray::*intarray_printYourself_1)() const = &IntArray::printYourself;
void pyclass_IntArray()
{
    class_<IntArray, boost::noncopyable>("IntArray")
        .def(init< optional<int> >())
        .def(init< IntArray& >())
        .def("resize", &IntArray::resize, "Checks size of receiver towards requested bounds. If dimension mismatch, size is adjusted accordingly")
        .def("giveSize", &IntArray::giveSize, "Returns size of receiver")
        .def("isEmpty", &IntArray::isEmpty, "Checks if receiver is empty (i.e., zero sized)")
        .def("containsOnlyZeroes", &IntArray::containsOnlyZeroes, "Checks if receiver is all zero")
        .def("minimum", &IntArray::minimum, "Finds the minimum component in the array")
        .def("contains", &IntArray::contains, "Retuns true if receiver contains given value")
        .def("add", &IntArray::add, "Adds given scalar to all values of receiver")
        .def("zero", &IntArray::zero, "Sets all component to zero")
        .def("printYourself", intarray_printYourself_1, "Prints receiver on stdout")

        .def("__len__", &IntArray::giveSize, "Returns size of receiver")
        .def("__getitem__", &IntArray::__getitem__, "Coefficient access function. Provides 0-based indexing access")
        .def("__setitem__", &IntArray::__setitem__, "Coefficient access function. Provides 0-based indexing access")
        .def("beCopyOf", &IntArray::beCopyOf, "Modifies receiver to become copy of given parameter")
        ;
}


/*****************************************************
* SparseMtrx
*****************************************************/
struct PySparseMtrx : SparseMtrx, wrapper<SparseMtrx>
{
    PySparseMtrx(int n=0, int m=0): SparseMtrx(n,m) {}

    void times(const FloatArray &x, FloatArray &answer) const override {
        this->get_override("times")();
    }

    void times(double x) override {
        this->get_override("times")();
    }

    int buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme& s) override {
        return this->get_override("buildInternalStructure")(eModel, di, s);
    }

    int buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme& r, const UnknownNumberingScheme &s) override {
        if (override f = this->get_override("buildInternalStructure")) {return f(eModel, di, r, s);}
        return this->get_override("buildInternalStructure")(eModel, di, r, s);
    }

    int assemble(const IntArray &loc, const FloatMatrix &mat) override {
        return this->get_override("assemble")();
    }

    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat) override {
        return this->get_override("assemble")();
    }

    bool canBeFactorized() const override {
        return this->get_override("canBeFactorized")();
    }

    void zero() override {
        this->get_override("zero")();
    }

    double &at(int i, int j) override {
        return this->get_override("at")();
    }

    double at(int i, int j) const override {
        return this->get_override("at")();
    }

    void toFloatMatrix(FloatMatrix &answer) const override {
        this->get_override("toFloatMatrix")();
    }

    void printYourself() const override {
        this->get_override("printYourself")();
    }

    SparseMtrxType giveType() const override {
        return this->get_override("giveType")();
    }

    bool isAsymmetric() const override {
        return this->get_override("isAsymmetric")();
    }

    FloatArray trans_mult(const FloatArray &x) const {
        return this->get_override("trans_mult")();
    }

    const char *giveClassName() const override {
        return this->get_override("giveClassName")();
    }
};

void (PySparseMtrx::*sm_times_1)(const FloatArray &x, FloatArray &ans) const  = &PySparseMtrx::times;
void (PySparseMtrx::*sm_times_2)(double x) = &PySparseMtrx::times;
int (PySparseMtrx::*sm_assemble_1)(const IntArray &loc, const FloatMatrix &mat) = &PySparseMtrx::assemble;
int (PySparseMtrx::*sm_assemble_2)(const IntArray &r, const IntArray& c, const FloatMatrix &mat) = &PySparseMtrx::assemble;

int (PySparseMtrx::*sm_buildInternalStructure_1)(EngngModel *eModel, int di, const UnknownNumberingScheme& s) = &PySparseMtrx::buildInternalStructure;
int (PySparseMtrx::*sm_buildInternalStructure_2)(EngngModel *eModel, int di, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) = &PySparseMtrx::buildInternalStructure;

void pyclass_SparseMtrx()
{
    class_<PySparseMtrx, boost::noncopyable>("SparseMtrx")
        //.def(init<>())
        .def("times", pure_virtual(sm_times_1))
        .def("timesd", pure_virtual(sm_times_2))
        .def("assemble", sm_assemble_1)
        .def("assembleu", sm_assemble_2)
        .def("buildInternalStructure", pure_virtual(sm_buildInternalStructure_1))
        .def("buildInternalStructure2", sm_buildInternalStructure_2)

        //.def("zero", zero)
        ;
}


/*****************************************************
* SpatialLocalizer
*****************************************************/
Element* (SpatialLocalizer::*giveElementContainingPoint_1)(const FloatArray &coords, const IntArray* regionList) = &SpatialLocalizer::giveElementContainingPoint;
void pyclass_SpatialLocalizer()
{
    class_<SpatialLocalizer, boost::noncopyable>("SpatialLocalizer", no_init)
      .def("giveElementContainingPoint", pure_virtual(giveElementContainingPoint_1), return_internal_reference<>())
        // XXX .def("giveElementCloseToPoint", &SpatialLocalizer::giveElementCloseToPoint, return_internal_reference<>())
        .def("giveElementClosestToPoint", &SpatialLocalizer::giveElementClosestToPoint, return_internal_reference<>())
        .def("init", &SpatialLocalizer::init)
        .def("giveClassName", &SpatialLocalizer::giveClassName)
        ;
}


void pyclass_UnknownNumberingScheme()
{
    class_<UnknownNumberingScheme,boost::noncopyable>("UnknownNumberingScheme",no_init)
        .def("giveDofEquationNumber",pure_virtual(&UnknownNumberingScheme::giveDofEquationNumber))
        .def("giveRequiredNumberOfDomainEquation",&UnknownNumberingScheme::giveRequiredNumberOfDomainEquation)
        .def("isDefault",&UnknownNumberingScheme::isDefault)
        .def("init",&UnknownNumberingScheme::init)
        ;
}

/*****************************************************
* EngngModelContext
*****************************************************/
void pyclass_EngngModelContext()
{
    class_<EngngModelContext, boost::noncopyable>("EngngModelContext", no_init)
        .def("giveFieldManager", &EngngModelContext::giveFieldManager, return_internal_reference<>())
        ;
}


/*****************************************************
* EngngModel
*****************************************************/
class PyEngngModel : public EngngModel, public wrapper<EngngModel>
{
public:
    //PyEngngModel(int i) : EngngModel(i, nullptr) {}
    PyEngngModel(int i, EngngModel *_master = nullptr) : EngngModel(i, _master) {}

    void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime) override {
        this->get_override("printDofOutputAt")();
    }

    void solveYourself() override {
        if (override f = this->get_override("solveYourself")) {f();return;}
        EngngModel::solveYourself();
    }
    void solveYourselfAt(TimeStep *t) override {
        if (override f = this->get_override("solveYourselfAt")) {f(t); return;}
        EngngModel::solveYourselfAt(t);
    }
    void terminate(TimeStep *t) override {
        if (override f = this->get_override("terminate")) {f(t); return;}
        EngngModel::terminate(t);
    }

    void updateYourself(TimeStep *t) override {
        if (override f = this->get_override("updateYourself")) {f(t); return;}
        EngngModel::updateYourself(t);
    }

    TimeStep* giveNextStep() override {
        if (override f = this->get_override("giveNextStep")) {return f();}
        return this->get_override("giveNextStep")();
    }

    void default_solveYourself() { return this->EngngModel::solveYourself(); }
    void default_solveYourselfAt(TimeStep* t) { return this->EngngModel::solveYourselfAt(t); }
    void default_terminate(TimeStep* t) { return this->EngngModel::terminate(t); }
    void default_updateYourself(TimeStep* t) { return this->EngngModel::updateYourself(t); }
    TimeStep* default_giveNextStep() { return this->EngngModel::giveNextStep(); }
    MetaStep* default_giveMetaStep(int i) { return this->EngngModel::giveMetaStep(i); }
};

void pyclass_EngngModel()
{
    class_<PyEngngModel, boost::noncopyable>("EngngModel", no_init)
        //.def(init<int, optional<EngngModel*> > ())
        //.def(init<int, char*, optional<EngngModel*> >())
        .def("giveDomain", &PyEngngModel::giveDomain, return_internal_reference<>())
        .def("setDomain", &EngngModel::setDomain)
        .def("giveNumberOfDomains", &PyEngngModel::giveNumberOfDomains)
        .add_property("numberOfDomains",&PyEngngModel::giveNumberOfDomains)
        .def("giveOutputBaseFileName", &PyEngngModel::giveOutputBaseFileName)
        .add_property("outputBaseFileName", &PyEngngModel::giveOutputBaseFileName)
        .def("giveDomainErrorEstimator", &PyEngngModel::giveDomainErrorEstimator, return_internal_reference<>())
        .def("terminateAnalysis", &PyEngngModel::terminateAnalysis)
        .def("solveYourself", &EngngModel::solveYourself, &PyEngngModel::default_solveYourself)
        .def("solveYourselfAt", &EngngModel::solveYourselfAt, &PyEngngModel::default_solveYourselfAt)
        .def("terminate", &EngngModel::terminate, &PyEngngModel::default_terminate)
        .def("updateYourself", &EngngModel::updateYourself, &PyEngngModel::default_updateYourself)
        .def("initializeYourself", &PyEngngModel::initializeYourself)
        .def("printDofOutputAt", pure_virtual(&EngngModel::printDofOutputAt))
        .def("printYourself", pure_virtual(&EngngModel::printYourself))
        .def("checkConsistency", &PyEngngModel::checkConsistency)
        .def("init", &PyEngngModel::init)
        .def("postInitialize", &PyEngngModel::postInitialize)
        .def("checkProblemConsistency", &PyEngngModel::checkProblemConsistency)
        .def("setRenumberFlag", &PyEngngModel::setRenumberFlag)
        .def("giveNumberOfSteps", &EngngModel::giveNumberOfSteps)
        .add_property("numberOfSteps",&PyEngngModel::giveNumberOfSteps)
        .def("giveCurrentStep", &EngngModel::giveCurrentStep,(boost::python::arg("force")=false), return_internal_reference<>())
        .def("givePreviousStep", &EngngModel::givePreviousStep,(boost::python::arg("force")=false), return_internal_reference<>())
        // .add_property("currentStep",make_function(&PyEngngModel::giveCurrentStep, return_internal_reference<>()))
        .def("giveNextStep",&EngngModel::giveNextStep, &PyEngngModel::default_giveNextStep, return_internal_reference<>())
        .def("generateNextStep", &EngngModel::generateNextStep, return_internal_reference<>()) // bp prototype
        .def("giveExportModuleManager",&EngngModel::giveExportModuleManager, return_internal_reference<>())
        .add_property("exportModuleManager",make_function(&PyEngngModel::giveExportModuleManager, return_internal_reference<>()))
        .def("giveContext", &PyEngngModel::giveContext, return_internal_reference<>())
        .def("giveField",&EngngModel::giveField)
        .add_property("context", make_function(&PyEngngModel::giveContext, return_internal_reference<>()))
        .def("giveMetaStep", &EngngModel::giveMetaStep, &PyEngngModel::default_giveMetaStep, return_internal_reference<>())
        .def("initMetaStepAttributes", &EngngModel::initMetaStepAttributes)
        .def("preInitializeNextStep", &EngngModel::preInitializeNextStep)
        .def("giveClassName", &EngngModel::giveClassName)
        ;
}

EngngModel *InstanciateProblem_1 (DataReader &dr, problemMode mode, int contextFlag)
{
    return InstanciateProblem(dr, mode, contextFlag).release();
}


/*****************************************************
* ExportModuleManager
*****************************************************/
void pyclass_ExportModuleManager()
{
    class_<ExportModuleManager, boost::noncopyable >("ExportModuleManager", no_init)
        .def("doOutput", &ExportModuleManager::doOutput)
        .def("giveModule", &ExportModuleManager::giveModule, return_internal_reference<>())
        ;
}


/*****************************************************
* ExportModule
*****************************************************/
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(exportmodule_overloads_dooutput, resize, 1, 2);
class PyExportModule : public ExportModule, public wrapper<ExportModule>
{
public:
    PyExportModule(int i, EngngModel * e) : ExportModule(i, e) {}

    void doOutput(TimeStep *tStep, bool forcedOutput=false) override {
        this->get_override("doOutput")();
    }

    void default_doOutput(TimeStep *tStep, bool forcedOutput=false) {this->ExportModule::doOutput(tStep,forcedOutput);}
};

void pyclass_ExportModule()
{
    class_<PyExportModule, boost::noncopyable>("ExportModule", no_init)
        .def("doOutput", pure_virtual(&ExportModule::doOutput), (bp::arg("tStep"), bp::arg("forcedOutput")=false ) )
        .def("doForcedOutput", &ExportModule::doForcedOutput)
        //.def("doOutput", pure_virtual(&ExportModule::doOutput), exportmodule_overloads_dooutput() )
        ;
}


/*****************************************************
* DataReader
*****************************************************/
struct PyDataReader : DataReader, wrapper<DataReader>
{
    InputRecord &giveInputRecord(InputRecordType irType, int recordId) override
    {
        return this->get_override("giveInputRecord")();
    }

    void finish() override {
        this->get_override("finish")();
    }

    std :: string giveReferenceName() const override {
        return this->get_override("giveReferenceName")();
    }
};

void pyclass_DataReader()
{
    class_<PyDataReader, boost::noncopyable>("DataReader")
        ;
}


/*****************************************************
* OOFEMTXTDataReader
*****************************************************/
void pyclass_OOFEMTXTDataReader()
{
    class_<OOFEMTXTDataReader, bases<DataReader> >("OOFEMTXTDataReader", init<char* >())
        ;
}


/*****************************************************
* Domain
*****************************************************/
void pyclass_Domain()
{
    class_<Domain, boost::noncopyable>("Domain", init<int, int, EngngModel* >())
        .def("giveNumber", &Domain::giveNumber)
        .def("setNumber", &Domain::setNumber)
        .add_property("number",&Domain::giveNumber,&Domain::setNumber)
        .def("giveElement", &Domain::giveElement, return_internal_reference<>())
        .def("giveEngngModel", &Domain::giveEngngModel, return_internal_reference<>())
        .def("giveLoad", &Domain::giveLoad, return_internal_reference<>())
        .def("giveBc", &Domain::giveBc, return_internal_reference<>())
        .def("giveIc", &Domain::giveIc, return_internal_reference<>())
        .def("giveFunction", &Domain::giveFunction, return_internal_reference<>())
        .def("giveMaterial", &Domain::giveMaterial, return_internal_reference<>())
        .def("giveCrossSection", &Domain::giveCrossSection, return_internal_reference<>())
        .def("giveNode", &Domain::giveNode, return_internal_reference<>())
        .def("giveDofManager", &Domain::giveDofManager, return_internal_reference<>())
        .def("giveMateriatal", &Domain::giveMaterial, return_internal_reference<>())
        .def("giveOutputManager", &Domain::giveOutputManager, return_internal_reference<>())
        .add_property("outputManager",make_function(&Domain::giveOutputManager, return_internal_reference<>()))

        .def("giveNumberOfDofManagers", &Domain::giveNumberOfDofManagers)
        .add_property("numberOfDofManagers",&Domain::giveNumberOfDofManagers)
        .def("giveNumberOfElements", &Domain::giveNumberOfElements)
        .add_property("numberOfElements",&Domain::giveNumberOfElements)
        .def("giveNumberOfMaterialModels", &Domain::giveNumberOfMaterialModels)
        .add_property("numberOfMaterialModels",&Domain::giveNumberOfMaterialModels)
        .def("giveNumberOfCrossSectionModels", &Domain::giveNumberOfCrossSectionModels)
        .add_property("numberOfCrossSectionModels",&Domain::giveNumberOfCrossSectionModels)
        .def("giveNumberOfBoundaryConditions", &Domain::giveNumberOfBoundaryConditions)
        .add_property("numberOfBoundaryConditions",&Domain::giveNumberOfBoundaryConditions)
        .def("giveNumberOfInitialConditions", &Domain::giveNumberOfInitialConditions)
        .add_property("numberOfInitialConditions",&Domain::giveNumberOfInitialConditions)
        .def("giveNumberOfFunctions", &Domain::giveNumberOfFunctions)
        .add_property("numberOfFunctions",&Domain::giveNumberOfFunctions)
        .def("giveNumberOfRegions", &Domain::giveNumberOfRegions)
        .add_property("numberOfRegions",&Domain::giveNumberOfRegions)

        .def("resizeDofManagers", &Domain::resizeDofManagers)
        .def("resizeElements", &Domain::resizeElements)
        .def("resizeCrossSectionModels", &Domain::resizeCrossSectionModels)
        .def("resizeMaterials", &Domain::resizeMaterials)
        .def("resizeBoundaryConditions", &Domain::resizeBoundaryConditions)
        .def("resizeInitialConditions", &Domain::resizeInitialConditions)
        .def("resizeFunctions", &Domain::resizeFunctions)

        .def("setDofManager", &Domain::py_setDofManager)
        .def("setElement", &Domain::py_setElement)
        .def("setCrossSection", &Domain::py_setCrossSection)
        .def("setMaterial", &Domain::py_setMaterial)
        .def("setBoundaryCondition", &Domain::py_setBoundaryCondition)
        .def("setInitialCondition", &Domain::py_setInitialCondition)
        .def("setFunction", &Domain::py_setFunction)

        .def("checkConsistency", &Domain::checkConsistency)
        .def("giveConnectivityTable", &Domain::giveConnectivityTable, return_internal_reference<>())
        .def("giveSpatialLocalizer", &Domain::giveSpatialLocalizer, return_internal_reference<>())
        .def("giveErrorEstimator", &Domain::giveErrorEstimator, return_internal_reference<>())
        .def("giveSmoother", &Domain::giveSmoother, return_internal_reference<>())

        .def("giveNumberOfSpatialDimensions", &Domain::giveNumberOfSpatialDimensions)
        .add_property("numberOfSpatialDimensions", &Domain::giveNumberOfSpatialDimensions)

        .def("checkConsistency", &Domain::checkConsistency)
        .def("giveArea", &Domain::giveArea)
        .def("giveVolume", &Domain::giveVolume)
        ;
}


/*****************************************************
* FEMComponent
*****************************************************/
void pyclass_FEMComponent()
{
    class_<FEMComponent, boost::noncopyable>("FEMComponent", no_init)
        .def("giveNumber", &FEMComponent::giveNumber)
        .def("setNumber", &FEMComponent::setNumber)
        .add_property("number", &FEMComponent::giveNumber, &FEMComponent::setNumber)
        .def("giveDomain", &FEMComponent::giveDomain, return_internal_reference<>())
        .def("setDomain", &FEMComponent::setDomain)
        .add_property("domain", make_function(&FEMComponent::giveDomain, return_internal_reference<>()), &FEMComponent::setDomain)
        .def("checkConsistency", &FEMComponent::checkConsistency)
        .def("printYourself", &FEMComponent::printYourself)
        ;
}


/*****************************************************
* DofManager
*****************************************************/
class PyDofManager : public DofManager, public wrapper<DofManager>
{
public:
    PyDofManager (int n, Domain *d) : DofManager (n,d) {}

    void giveUnknownVector(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *stepN, bool padding) {
        if ( override f = this->get_override("giveUnknownVector") ) {
            f (answer, dofMask, mode, stepN, padding);
            return;
        }
        DofManager::giveUnknownVector(answer, dofMask, mode, stepN);
    }
    bool computeL2GTransformation(FloatMatrix &answer, const IntArray &dofIDArry) override {
        if (override f = this->get_override("computeL2GTransformation")) { return f(answer, dofIDArry); }
        return DofManager::computeL2GTransformation(answer, dofIDArry);
    }

    void default_giveUnknownVector(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *stepN, bool padding)
    { return this->DofManager::giveUnknownVector(answer, dofMask, mode, stepN, padding); }
};

void (DofManager::*giveUnknownVector_1)(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *stepN, bool padding) = &DofManager::giveUnknownVector;

void pyclass_DofManager()
{
    class_<DofManager, bases<FEMComponent>, boost::noncopyable>("DofManager", no_init)
        //.def("hasCoordinates", &DofManager::hasCoordinates)
        // TODO return type (copy rather than pointer?)
        .def("giveCoordinates", &DofManager::giveCoordinates, return_internal_reference<>())
        .add_property("coordinates", make_function(&DofManager::giveCoordinates, return_internal_reference<>()))
        .def("giveCoordinate", &DofManager::giveCoordinate)
        // TODO return type (copy rather than pointer?)
        .def("giveLoadArray", &DofManager::giveLoadArray, return_internal_reference<>())
        // TODO return type (copy rather than pointer?)
        .add_property("loadArray", make_function(&DofManager::giveLoadArray, return_internal_reference<>()))
        .def("giveLabel", &DofManager::giveLabel)
        .add_property("label",&DofManager::giveLabel)
        .def("giveUnknownVector", giveUnknownVector_1)
        .def("computeL2GTransformation", &DofManager::computeL2GTransformation)
        ;
}


/*****************************************************
* Element
*****************************************************/
class PyElement : public Element, public wrapper<Element>
{
public:
    PyElement (int n, Domain *d) : Element (n,d) {}

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) override {
        this->get_override("giveCharacteristicMatrix")();
    }

    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override {
        this->get_override("giveCharacteristicVector")();
    }

    double giveCharacteristicValue(CharType, TimeStep *) override {
        return this->get_override("giveCharacteristicValue")();
    }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override {
        this->get_override("giveDofManDofIDMask")();
    }

    int evaluateAt(FloatArray &answer, FloatArray &coords,
            ValueModeType mode, TimeStep *atTime) {
        return this->get_override("evaluateAt")(answer, coords,mode,atTime);
    }
};


void (PyElement::*giveLocationArray_1)(IntArray &locationArray, const UnknownNumberingScheme &s, IntArray *dofIds) const  = &Element::giveLocationArray;
void (PyElement::*giveLocationArray_2)(IntArray &locationArray, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIds) const = &Element::giveLocationArray;

void pyclass_Element()
{
    class_<PyElement, bases<FEMComponent>, boost::noncopyable>("Element", no_init)
        .def("giveNumberOfDofs", &PyElement::giveNumberOfDofs)
        .add_property("numberOfDofs", &PyElement::giveNumberOfDofs)
        .def("giveNumberOfInternalDofManagers", &PyElement::giveNumberOfInternalDofManagers)
        .add_property("numberOfInternalDofManagers", &PyElement::giveNumberOfInternalDofManagers)
        .def("giveInternalDofManager", &PyElement::giveInternalDofManager, return_internal_reference<>())
        .def("giveCharacteristicMatrix", &PyElement::giveCharacteristicMatrix)
        .def("giveCharacteristicVector", &PyElement::giveCharacteristicVector)
        .def("giveCharacteristicValue", &PyElement::giveCharacteristicValue)
        .def("computeMeanSize", &PyElement::computeMeanSize)
        .add_property("meanSize", &PyElement::computeMeanSize)
        .def("computeVolume", &PyElement::computeVolume)
        .add_property("volume", &PyElement::computeVolume)
        .def("computeArea", &PyElement::computeArea)
        .add_property("area", &PyElement::computeArea)
        .def("computeLength", &PyElement::computeLength)
        .add_property("length", &PyElement::computeLength)
        .def("giveLabel", &Element::giveLabel)
        .add_property("label",&Element::giveLabel)
        .def("giveLocationArray", giveLocationArray_1)
        .def("giveLocationArray2", giveLocationArray_2)
        .def("giveNumberOfDofManagers", &PyElement::giveNumberOfDofManagers)
        .add_property("numberOfDofManagers", &PyElement::giveNumberOfDofManagers)
        .def("computeNumberOfDofs", &PyElement::computeNumberOfDofs)
        .def("giveNumberOfNodes", &PyElement::giveNumberOfNodes)
        .def("giveGeometryType", &PyElement::giveGeometryType)
        .add_property("geometryType", &PyElement::giveGeometryType)
        .def("giveDofManagerNumber", &Element::giveDofManagerNumber)
        .def("giveDofManager", &Element::giveDofManager, return_internal_reference<>())
        .def("giveNode", &Element::giveNode, return_internal_reference<>())
          .def("giveClassName",&Element::giveClassName)
        // TODO return type (copy rather than pointer?)
        .def("giveDofManArray", &Element::giveDofManArray, return_value_policy<reference_existing_object>())
        .def("giveNumberOfIntegrationRules", &Element::giveNumberOfIntegrationRules)
        .add_property("numberOfIntegrationRules", &PyElement::giveNumberOfIntegrationRules)
        .def("giveDefaultIntegrationRulePtr", &Element::giveDefaultIntegrationRulePtr, return_internal_reference<>())
        .add_property("defaultIntegrationRule",make_function(&PyElement::giveDefaultIntegrationRulePtr, return_internal_reference<>()))
        .def("giveMaterial", &Element::giveMaterial, return_internal_reference<>())
        .def("giveCrossSection", &Element::giveCrossSection, return_internal_reference<>())
        .def("setCrossSection", &Element::setCrossSection)
        .add_property("crossSection",make_function(&PyElement::giveCrossSection, return_internal_reference<>()), &PyElement::setCrossSection)
        .def("giveRegionNumber", &Element::giveRegionNumber)
        .add_property("regionNumber", &PyElement::giveRegionNumber)
        ;
}


/*****************************************************
* Material
*****************************************************/
void pyclass_Material()
{
    class_<Material, bases<FEMComponent>, boost::noncopyable >("Material", no_init)
        .def("giveIPValue", &Material::giveIPValue)
        .def("setIPValue", &Material::setIPValue)
        .def("giveStatus", &Material::giveStatus, return_internal_reference<>())
        ;
}


/*****************************************************
* StructuralMaterial
*****************************************************/
struct PyStructuralMaterial : StructuralMaterial , wrapper<StructuralMaterial>
{
    PyStructuralMaterial(int i, Domain *d) : StructuralMaterial(i,d) {}
    void giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) const override {
        if (override f = this->get_override("giveRealStressVector")) { f(answer,gp,reducedStrain,tStep); }
        this->get_override("giveRealStressVector")(answer,gp,reducedStrain,tStep);
    }
};
void pyclass_StructuralMaterial()
{
    class_<PyStructuralMaterial, bases<Material>, boost::noncopyable >("StructuralMaterial", no_init)
        .def("giveRealStressVector", pure_virtual( &StructuralMaterial::giveRealStressVector))
        .def("giveStiffnessMatrix", &StructuralMaterial::giveStiffnessMatrix)
        ;
}

/*****************************************************
* StructuralMaterialSettable //
*****************************************************/
void pyclass_StructuralMaterialSettable()
{
    class_<StructuralMaterialSettable, bases<StructuralMaterial>, boost::noncopyable >("StructuralMaterialSettable", no_init)
        .def("giveIPValue", &StructuralMaterialSettable::giveIPValue)
        .def("setIPValue", &StructuralMaterialSettable::setIPValue)
        ;
}


StructuralMaterial* material2structuralMaterial(Material *mat) { return (StructuralMaterial*)mat; }
StructuralMaterialSettable* material2structMatSettable(Material *mat) { return (StructuralMaterialSettable*)mat; }




/*****************************************************
* CrossSection
*****************************************************/
void pyclass_CrossSection()
{
    class_<CrossSection, bases<FEMComponent>, boost::noncopyable >("CrossSection", no_init)
        ;
}


/*****************************************************
* GeneralBoundaryCondition
*****************************************************/
void pyclass_GeneralBoundaryCondition()
{
    class_<GeneralBoundaryCondition, bases<FEMComponent>, boost::noncopyable >("GeneralBoundaryCondition", no_init)
        ;
}


/*****************************************************
* BoundaryCondition
*****************************************************/
void pyclass_BoundaryCondition()
{
    class_<BoundaryCondition, bases<GeneralBoundaryCondition>, boost::noncopyable >("BoundaryCondition", init<int, Domain*>())
        .def("setPrescribedValue", &BoundaryCondition::setPrescribedValue)
        // XXX .def("give", &BoundaryCondition::give)
        .def("isImposed", &BoundaryCondition::isImposed)
        ;
}


/*****************************************************
* Load
*****************************************************/
struct PyLoad : Load , wrapper<Load>
{
    PyLoad(int i, Domain *d) : Load(i,d) {}
    void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode) override {
        if (override f = this->get_override("computeValueAt")) { f(answer,tStep,coords,mode);}
        this->get_override("computeValueAt")(answer,tStep,coords,mode);
    }
    FloatArray* GiveComponentArray_copy() { return new FloatArray(giveComponentArray()); }
};


void pyclass_Load()
{
    class_<PyLoad, bases<GeneralBoundaryCondition>, boost::noncopyable >("Load", no_init)
        .def("setComponentArray", &Load::setComponentArray)
        .def("giveComponentArray", &PyLoad::GiveComponentArray_copy,return_value_policy<manage_new_object>())
        .def("computeValueAt", pure_virtual( &Load::computeValueAt))
        ;

    class_<ForeignTemperatureFieldLoad,bases<PyLoad>>("ForeignTemperatureFieldLoad",no_init)
        .def_readwrite("foreignField",&ForeignTemperatureFieldLoad::foreignField)
        ;
}


/*****************************************************
* LoadTimeFunction
*****************************************************/
void pyclass_Function()
{
    class_<Function, bases<FEMComponent>, boost::noncopyable >("Function", no_init)
        ;
}


/*****************************************************
* MaterialStatus
*****************************************************/
void pyclass_MaterialStatus()
{
    class_<MaterialStatus, boost::noncopyable >("MaterialStatus", no_init)
        .def("updateYourself", &MaterialStatus::updateYourself)
        ;
}


/*****************************************************
* StructuralMaterialStatus
*****************************************************/
void pyclass_StructuralMaterialStatus()
{
    class_<StructuralMaterialStatus, bases<MaterialStatus>, boost::noncopyable >("StructuralMaterialStatus", no_init)
        ;
}


/*****************************************************
* TimeStep
*****************************************************/
void pyclass_TimeStep()
{
    class_<TimeStep, boost::noncopyable>("TimeStep",no_init)
        .def("giveTargetTime", &TimeStep::giveTargetTime)
        .def("setTargetTime", &TimeStep::setTargetTime)
        .add_property("targetTime",&TimeStep::giveTargetTime, &TimeStep::setTargetTime)
        .def("giveIntrinsicTime", &TimeStep::giveIntrinsicTime)
        .def("setIntrinsicTime", &TimeStep::setIntrinsicTime)
        .add_property("intrinsicTime",&TimeStep::giveIntrinsicTime, &TimeStep::setIntrinsicTime)
        .def("giveTimeIncrement", &TimeStep::giveTimeIncrement)
        .def("setTimeIncrement", &TimeStep::setTimeIncrement)
        .add_property("timeIncrement",&TimeStep::giveTimeIncrement, &TimeStep::setTimeIncrement)
        .def("giveNumber", &TimeStep::giveNumber)
        .def("setNumber", &TimeStep::setNumber)
        .add_property("number",&TimeStep::giveNumber, &TimeStep::setNumber)
        ;
}

/*****************************************************
* MetaStep
*****************************************************/
void pyclass_MetaStep()
{
    class_<MetaStep, boost::noncopyable>("MetaStep",no_init)
        .def("giveNumber", &MetaStep::giveNumber)
        .add_property("number", &MetaStep::giveNumber)
        ;
}


/*****************************************************
* DataStream
*****************************************************/
void pyclass_DataStream()
{
    class_<DataStream, boost::noncopyable>("DataStream", no_init)
        ;
}

/*****************************************************
* Field
*****************************************************/
#if 0
class PyField : public Field, public wrapper<Field>
{
public:
    PyField (FieldType b): Field(b) {}
    int evaluateAt(FloatArray &answer, FloatArray &coords,
            ValueModeType mode, TimeStep *atTime) override {
        return this->get_override("evaluateAt")(answer, coords,mode,atTime);
    }
    int evaluateAt(FloatArray &answer, DofManager* dman,
            ValueModeType mode, TimeStep *atTime) override {
        return this->get_override("evaluateAt")(answer,dman,mode,atTime);
    }
    void saveContext(DataStream &stream, ContextMode mode) override {
        return this->get_override("saveContext")(stream, mode);
    }
    void restoreContext(DataStream &stream, ContextMode mode) override {
        return this->get_override("restoreContext")(stream, mode);
    }
};

// these two need to copy the FloatArray object, which is noncopyable and no by-value converter is found
FloatArray Field_evaluateAtDman(Field* field, DofManager* dman, ValueModeType mode, TimeStep* atTime){
    FloatArray ret;
    int err=field->evaluateAt(ret,dman,mode,atTime);
    if(err) throw std::runtime_error("Error evaluating field at given DofManager.");
    return ret;
}
FloatArray Field_evaluateAtPos(Field* field, FloatArray& coords, ValueModeType mode, TimeStep* atTime){
    FloatArray ret;
    int err=field->evaluateAt(ret,coords,mode,atTime);
    if(err) throw std::runtime_error("Error evaluating field at given position.");
    return ret;
}
    // these force Python user to check error code, which is very unpythonic
    int (Field::*Field_evaluateAtPos)(FloatArray &answer, FloatArray &coords, ValueModeType mode, TimeStep *atTime) = &Field::evaluateAt;
    int (Field::*Field_evaluateAtDman)(FloatArray &answer, DofManager* dman, ValueModeType mode, TimeStep *atTime) = &Field::evaluateAt;
}
#endif

// raise exception if there is a problem
FloatArray Field_evaluateAtDman(Field* field, DofManager* dman, ValueModeType mode, TimeStep* atTime)
{
    FloatArray answer;
    int err=field->evaluateAt(answer,dman,mode,atTime);
    if(err) throw std::runtime_error("Error evaluating field at given DofManager.");
    return answer;
}
FloatArray Field_evaluateAtPos(Field* field, FloatArray coords, ValueModeType mode, TimeStep* atTime)
{
    FloatArray answer;
    int err=field->evaluateAt(answer,coords,mode,atTime);
    if(err) throw std::runtime_error("Error evaluating field at given position.");
    return answer;
}


#ifndef OOFEM_USE_DEPRECATED_AUTOPTR
    std::shared_ptr<Field> DofManValueField_create (FieldType t, Domain* d) { return std::make_shared<DofManValueField>(t, d); }
#else
    std::auto_ptr<Field> DofManValueField_create (FieldType t, Domain* d) { return std::auto_ptr<Field>( new DofManValueField(t, d)); }
#endif

void pyclass_Field()
{
    //this class is held by auto_ptr:  to allow for taking ownership of a raw pointer
    //see http://www.boost.org/doc/libs/1_47_0/libs/python/doc/v2/faq.html#ownership for detail
    // also see http://stackoverflow.com/questions/4112561/boost-python-ownership-of-pointer-variables
    class_<Field, FieldPtr, boost::noncopyable>("Field", no_init)
        .def("evaluateAtPos",Field_evaluateAtPos,(bp::arg("coords"),bp::arg("mode")=VM_Unknown,bp::arg("atTime")=bp::object()))
        /*
        making atTime=py::object() (i.e. nullptr) by default leads to crash, look at that better;
        the default should use the current timestep or anything that is meaningful; mandatory for now
        */
        .def("evaluateAtDman",Field_evaluateAtDman,(bp::arg("dman"),bp::arg("mode"),bp::arg("atTime")))
        .def("giveType", &Field::giveType)
        ;

    class_<UniformGridField,std::shared_ptr<UniformGridField>,bases<Field>>("UniformGridField")
        .def("setGeometry",&UniformGridField::setGeometry,(bp::arg("lo"),bp::arg("hi"),bp::arg("div")))
        .def("setValues",&UniformGridField::setValues,(bp::arg("values")))
        ;

    class_<UnstructuredGridField, std::shared_ptr<UnstructuredGridField>, bases<Field>>("UnstructuredGridField", init<int,int,optional<double>>())
      .def("addVertex",&UnstructuredGridField::addVertex,(bp::arg("num"), bp::arg("coords")))
      .def("setVertexValue", &UnstructuredGridField::setVertexValue, (bp::arg("num"), bp::arg("vv")))
      .def("addCell", &UnstructuredGridField::addCell, (bp::arg("num"), bp::arg("type"), bp::arg("vertices")))
      ;
}

/*****************************************************
* FieldManager
*****************************************************/
//
// boost python DOES support std::shared_ptr
// references are released automatically if refcount drops to zero
//
void pyclass_FieldManager()
{
    class_<FieldManager, boost::noncopyable>("FieldManager", no_init)
        .def("registerField", &FieldManager::registerField)
        .def("isFieldRegistered", &FieldManager::isFieldRegistered)
        .def("giveField", &FieldManager::giveField)
        .def("unregisterField", &FieldManager::unregisterField)
        .def("giveRegisteredKeys",&FieldManager::giveRegisteredKeys) // returning std::vector<FieldPtr> as list is handled by custom converter
        ;
}


/*****************************************************
* DofManValueField
*****************************************************/
void pyclass_DofManValueField()
{
    class_<DofManValueField, bases<Field>, boost::noncopyable >("DofManValueField", no_init)
        .def("evaluateAt", Field_evaluateAtPos) //??
        .def("setDofManValue",&DofManValueField::setDofManValue)
        ;
}


/*****************************************************
* IntegrationRule
*****************************************************/
void pyclass_IntegrationRule()
{
    class_<IntegrationRule, boost::noncopyable >("IntegrationRule", init<int, Element*>())
        .def("getIntegrationPoint", &IntegrationRule::getIntegrationPoint, return_internal_reference<>())
        .def("giveNumberOfIntegrationPoints", &IntegrationRule::giveNumberOfIntegrationPoints)
        ;
}


/*****************************************************
* GaussPoint
*****************************************************/
void pyclass_GaussPoint()
{
    class_<GaussPoint, boost::noncopyable >("GaussPoint", no_init)
        .def("giveNumber", &GaussPoint::giveNumber)
        .add_property("number", &GaussPoint::giveNumber)
        .def("giveElement", &GaussPoint::giveElement, return_internal_reference<>())
        .add_property("element",make_function(&GaussPoint::giveElement, return_internal_reference<>()))
        .def("giveMaterial", &GaussPoint::giveMaterial, return_internal_reference<>())
        .add_property("material",make_function(&GaussPoint::giveMaterial, return_internal_reference<>()))
        ;
}


/*****************************************************
* OutputManager
*****************************************************/
void pyclass_OutputManager()
{
    class_<OutputManager, boost::noncopyable >("OutputManager", no_init)
        ;
}


/*****************************************************
* ClassFactory
*****************************************************/
void pyclass_ClassFactory()
{
    class_<ClassFactory, boost::noncopyable >("ClassFactory", no_init)
        ;
}







/*****************************************************
*
* O O F E M L I B   E N U M E R A T I O N S
*
*****************************************************/


/*****************************************************
* problemMode
*****************************************************/
void pyenum_problemMode()
{
    enum_<problemMode>("problemMode")
        .value("_processor", _processor)
        .value("_postProcessor", _postProcessor)
        ;
}


/*****************************************************
* Element_Geometry_Type
*****************************************************/
void pyenum_Element_Geometry_Type()
{
    enum_<Element_Geometry_Type>("Element_Geometry_Type")
        .value("EGT_line_1", EGT_line_1)
        .value("EGT_line_2", EGT_line_2) /* line element with three nodes 1---2---3 */
        .value("EGT_triangle_1",EGT_triangle_1) /* triangle element with three nodes */
        .value("EGT_triangle_2", EGT_triangle_2) /* triangle element with 6 nodes */
        .value("EGT_quad_1",EGT_quad_1)   /* quadrialateral with 4 nodes */
        .value("EGT_quad_2", EGT_quad_2)   /* quadratic quadrialateral with 8 nodes */
        .value("EGT_tetra_1", EGT_tetra_1)  /* tetrahedron with 4 nodes */
        .value("EGT_hexa_1", EGT_hexa_1)   /* hexahedron with 8 nodes */
        .value("EGT_hexa_2", EGT_hexa_2)   /* hexahedron with 20 nodes */
        .value("EGT_Composite", EGT_Composite)/* Composite" geometry, vtk export supported by individual elements */
        .value("EGT_unknown", EGT_unknown)  /* unknown" element geometry type */
        ;
}


/*****************************************************
* ValueModeType
*****************************************************/
void pyenum_ValueModeType()
{
    enum_<ValueModeType>("ValueModeType")
        .value("VM_Unknown", VM_Unknown)
        .value("VM_Total", VM_Total)
        .value("VM_Velocity", VM_Velocity)
        .value("VM_Acceleration", VM_Acceleration)
        .value("VM_Incremental", VM_Incremental)
        ;
}


/*****************************************************
* DofManTransfType
*****************************************************/
void pyenum_DofManTransfType()
{
    enum_<DofManTransfType>("DofManTransfType")
        .value("_toGlobalCS", _toGlobalCS)
        .value("_toNodalCS", _toNodalCS)
        ;
}


/*****************************************************
* DofIDItem
*****************************************************/
void pyenum_DofIDItem()
{
    enum_<DofIDItem>("DofIDItem")
        .value("Undef", Undef)
        .value("D_u", D_u)
        .value("D_v", D_v)
        .value("D_w", D_w)
        .value("R_u", R_u)
        .value("R_v", R_v)
        .value("R_w", R_w)
        .value("V_u", V_u)
        .value("V_v", V_v)
        .value("V_w", V_w)
        .value("T_f", T_f)
        .value("P_f", P_f)
        ;
}


/*****************************************************
* FieldType
*****************************************************/
void pyenum_FieldType()
{
    enum_<FieldType>("FieldType")
        .value("FT_Unknown", FT_Unknown)
        .value("FT_Velocity", FT_Velocity)
        .value("FT_Displacements", FT_Displacements)
        .value("FT_VelocityPressure", FT_VelocityPressure)
        .value("FT_Pressure", FT_Pressure)
        .value("FT_Temperature", FT_Temperature)
        .value("FT_HumidityConcentration", FT_HumidityConcentration)
        .value("FT_TransportProblemUnknowns", FT_TransportProblemUnknowns)
        .value("FT_TemperatureAmbient", FT_TemperatureAmbient)
        ;
}


/*****************************************************
* InternalStateType
*****************************************************/
void pyenum_InternalStateType()
{
    enum_<InternalStateType>("InternalStateType")
        .value("IST_StressTensor", IST_StressTensor)
        .value("IST_StrainTensor", IST_StrainTensor)
        .value("IST_PlasticStrainTensor", IST_PlasticStrainTensor)
        ;
}


/*****************************************************
* MatResponseMode
*****************************************************/
void pyenum_MatResponseMode()
{
    enum_<MatResponseMode>("MatResponseMode")
        .value("TangentStiffness", TangentStiffness)
        .value("SecantStiffness", SecantStiffness)
        .value("ElasticStiffness", ElasticStiffness)
        ;
}


/*****************************************************
* domainType
*****************************************************/
void pyenum_domainType()
{
    enum_<domainType>("domainType")
        .value("_unknownMode", _unknownMode)
        .value("_2dPlaneStressMode", _2dPlaneStressMode)
        .value("_PlaneStrainMode", _PlaneStrainMode)
        .value("_2dPlaneStressRotMode", _2dPlaneStressRotMode)
        .value("_3dMode", _3dMode)
        .value("_3dAxisymmMode", _3dAxisymmMode)
        .value("_2dMindlinPlateMode", _2dMindlinPlateMode)
        .value("_3dShellMode", _3dShellMode)
        .value("_2dTrussMode", _2dTrussMode)
        .value("_1dTrussMode", _1dTrussMode)
        .value("_2dBeamMode", _2dBeamMode)
        .value("_HeatTransferMode", _HeatTransferMode)
        .value("_HeatMass1Mode", _HeatMass1Mode)
        .value("_2dIncompressibleFlow", _2dIncompressibleFlow)
        .value("_3dIncompressibleFlow", _3dIncompressibleFlow)
        .value("_2dLatticeMode", _2dLatticeMode)
        ;
}








/*****************************************************
*
* O O F E M L I B   " C O N S T R U C T O R S "
*
*****************************************************/


/*****************************************************
* Auxiliary functions
*****************************************************/

/*
make oofem input line from **kw arguments. This function is used by "constructor" methods
e.g. in function:
  engngModel("nonLinearStatic",nSteps=10) # kw = {'nSteps':10}
*/
OOFEMTXTInputRecord makeOOFEMTXTInputRecordFrom(bp::dict &kw)
{
    //bp::dict temp(import("__main__").attr("__dict__"));
    //Python3
    bp::object temp(import("__main__").attr("__dict__"));
    temp["kw"] = kw;
    str command =
        "ret = ''\n"
        "for key,val in kw.items():\n" // iterate over key,val pairs
        "  if key.lower()=='number' or key.lower()=='domain': continue\n" // do not include "number" and "domain" kws
        "  if key=='f_t': key='f(t)'\n" // handle f(t) loadTimeFunction field name
        "  ret += ' %s'%key\n" // add key
        "  if val is True: pass\n" // e.g. tstep_all=True will produce only 'step_all'
        "  elif isinstance(val,(int,float,str)): ret += ' %s'%val\n" // add val if it is int, float or str
        "  elif isinstance(val,(list,tuple)):\n" // if val is tuple or list
        "    ret += ' %d'%len(val)\n" // add its length first
        "    for v in val:\n" // and then its values
        "      if isinstance(v,(int,float,str)): ret += ' %s'%v\n" // int, float, str
        "      else: ret += ' %d'%(v.giveNumber())\n" // or any other object, that have giveNumber method (e.g. nodes can be passed to elemnt function as they are, without the need of extracting their numbers first
        "  else: ret += ' %d'%(val.giveNumber())\n\n" // add arbitrary object with giveNumber method
        "ret = ret.lower()\n" // finally make it lower case
        "print(ret)\n"
        ;
    
    //Execute Python commands using boost::python::exec()
    exec(command,temp,temp);
    // extract string from globals["ret"], convert it to char* and return OOFEMTXTInputRecord from it
    //return OOFEMTXTInputRecord( ( extract<string>(temp["ret"])() ).c_str() );
    return OOFEMTXTInputRecord( 0, ( extract<string>(temp["ret"])() ) );
}

OOFEMTXTInputRecord makeOutputManagerOOFEMTXTInputRecordFrom(bp::dict &kw)
{
    bp::list keys = kw.keys();
    bp::list vals = kw.values();
    kw.clear();
    for (int i=0; i<len(keys); i++) {
        str s = extract<str>(keys[i])().lower();
        if (s=="tstep_all" || s=="tstep_step" || s=="tsteps_out" || s=="dofman_all" || s=="dofman_output" || s=="dofman_except" || s=="element_all" || s=="element_output" || s=="element_except") {
            kw[s] = vals[i];
        }
        kw[""] = "outputmanager";
    }
    return makeOOFEMTXTInputRecordFrom(kw);
}

// transform all keys of given dictionary to lower case
void makeDictKeysLowerCase(bp::dict &kw)
{
    bp::list keys = kw.keys();
    bp::list vals = kw.values();
    kw.clear();
    for (int i=0; i<len(keys); i++) {
        kw[extract<str>(keys[i])().lower()] = vals[i];
    }
}

/*
The process is almost same for all classes, therefore only Element part is documented line by line
*/

/*****************************************************
* EngngModel
*****************************************************/
// engngModel(aClass,number=0,master=None,**kw)
object engngModel(bp::tuple args, bp::dict kw)
{
    //args
    string aClass = extract<string>(args[0])();
    int number = len(args)>1? extract<int>(args[1])() : 0;
    EngngModel* master = len(args)>2? extract<EngngModel*>(args[2])() : nullptr;
    std::unique_ptr<EngngModel> engngm = classFactory.createEngngModel(aClass.c_str(),number,master);
    if ( !engngm ) { OOFEM_LOG_ERROR("engngModel: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    // instanciateYourself
    ///@todo Output filename isn't stored like this (and has never been!)!?
    string outFile;
    if ( ir.hasField("outfile") ) {
       ir.giveField(outFile, "outfile");
    } else {
       outFile = "oofem.out.XXXXXX";
    }
    //engngm->Instanciate_init(outFile.c_str(), engngm->giveNumberOfDomains());
    engngm->letOutputBaseFileNameBe(outFile);
    engngm->initializeFrom(ir);

    if ( ir.hasField(_IFT_EngngModel_nmsteps) ) {
      OOFEM_LOG_ERROR("engngModel: simulation with metasteps is not (yet) supported in Python");
    } else {
      engngm->instanciateDefaultMetaStep(ir);
    }

    engngm->Instanciate_init();
    //
    object ret = object(ptr(engngm.release()));
    /* ????????????????????
    // sets the last created engngModel as default one for further script
    temp_global["defaultEngngModel"] = ret;
    ???????????????????? */
    return ret;
}

object createEngngModelOfType(const char* type, bp::tuple args, bp::dict kw)
{
    args = len(args)>1? bp::make_tuple(type,args[0],args[1]) : len(args)>0? bp::make_tuple(type,args[0]) : bp::make_tuple(type);
    return engngModel(args,kw);
}
object linearStatic(bp::tuple args, bp::dict kw) { return createEngngModelOfType("linearstatic", args, kw); }


/*****************************************************
* Domain
*****************************************************/
// domain(number=1,serialNumber=1,engngModel*=None,dType=_unknownMode,**kw)
object domain(bp::tuple args, bp::dict kw)
{
    // args
    int number =             len(args)>0? extract<int>(args[0])() : 0;
    int serialNumber =       len(args)>1? extract<int>(args[1])() : 0;
    EngngModel *engngModel = len(args)>2? extract<EngngModel*>(args[2])() : nullptr;
    domainType dType =       len(args)>3? extract<domainType>(args[3])() : _unknownMode;
    auto d = std::make_unique<Domain>(number,serialNumber,engngModel);
    d->setDomainType(dType);
    // output manager record
    OOFEMTXTInputRecord omir = makeOutputManagerOOFEMTXTInputRecordFrom(kw);
    d->giveOutputManager()->initializeFrom(omir);
    object ret = object(ptr(d.release()));
    /* ????????????????????
    // sets the last created domain as default one for furtherscript
    temp_global["defaultDomain"] = ret;
    ???????????????????? */
    return ret;
}


/*****************************************************
* Element
*****************************************************/
// element(aClass,domain=defaultDomain,**kw)
object element(bp::tuple args, bp::dict kw)
{
    // extracts first python argument (string element type)
    string aClass = extract<string>(args[0])();
    // extracts values from args if they are specified
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : nullptr;
    /* ????????????????
    // if no material is specified, set it to 1
    if (!kw.has_key("mat")) { kw["mat"] = 1; }
    // if no cross section is specified, set it to 1
    if (!kw.has_key("crosssect")) { kw["crosssect"] = 1; }
    // if no domain is specified and one already exists in the script, use that one
    if (domain == nullptr && temp_global.has_key("defaultDomain")) { domain = extract<Domain*>(temp_global["defaultDomain"])(); }
    if (domain==nullptr) { LOG_ERROR(oofem_errLogger,"wrong Domain"); }
    ???????????????????? */
    // create Element (convert aClass to char* - expected by classFactory.createElement)
    auto elem = classFactory.createElement(aClass.c_str(),number,domain);
    // if elem==nullptr, something was wrong
    if (!elem) { OOFEM_LOG_ERROR("element: wrong input data"); }
    // sets globalNumber == number befor initializeFrom
    elem->setGlobalNumber(number);
    // construct OOFEMTXTInputRecord from bp::dict **kw
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    // pass input record to elem
    elem->initializeFrom(ir);
    // convert element to PyObject (expected by raw_function, which enables fun(*args,**kw) syntax in python)
    return object(ptr(elem.release()));
}

// auxiliary constructor for specific element type (name is passed as first argument)
object createElementOfType(const char* type, bp::tuple args, bp::dict kw)
{
    args = len(args)>1? bp::make_tuple(type,args[0],args[1]) : len(args)>0? bp::make_tuple(type,args[0]) : bp::make_tuple(type);
    return element(args, kw);
}
// specific elements
object beam2d(bp::tuple args, bp::dict kw) { return createElementOfType("beam2d",args,kw); }



/*****************************************************
* DofManager
*****************************************************/
// dofManager(aClass,domain=defaultDomain,**kw)
object dofManager(bp::tuple args, bp::dict kw)
{
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : nullptr;
    auto dofMan = classFactory.createDofManager(aClass.c_str(),number,domain);
    if (!dofMan) { OOFEM_LOG_ERROR("dofManager: wrong input data"); }
    dofMan->setGlobalNumber(number);
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    dofMan->initializeFrom(ir);
    return object(ptr(dofMan.release()));
}

object createDofManagerOfType(const char* type, bp::tuple args, bp::dict kw)
{
    args = len(args)>1? bp::make_tuple(type,args[0],args[1]) : len(args)>0? bp::make_tuple(type,args[0]) : bp::make_tuple(type);
    return dofManager(args, kw);
}
object node(bp::tuple args, bp::dict kw) { return createDofManagerOfType("node",args,kw); }


/*****************************************************
* GeneralBoundaryCondition
*****************************************************/
// generalBoundaryCondition(aClass,domain=defaultDomain,**kw)
object generalBoundaryCondition(bp::tuple args, bp::dict kw)
{
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : nullptr;
    auto bc = classFactory.createBoundaryCondition(aClass.c_str(),number,domain);
    if (!bc) { OOFEM_LOG_ERROR("generalBoundaryCondition: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    bc->initializeFrom(ir);
    return object(ptr(bc.release()));
}

object CreateBCOfType(const char* type, bp::tuple args, bp::dict kw)
{
    args = len(args)>1? bp::make_tuple(type,args[0],args[1]) : len(args)>0? bp::make_tuple(type,args[0]) : bp::make_tuple(type);
    return generalBoundaryCondition(args, kw);
}
object boundaryCondition(bp::tuple args, bp::dict kw) { return CreateBCOfType("boundarycondition",args,kw); }
object constantEdgeLoad(bp::tuple args, bp::dict kw) { return CreateBCOfType("constantedgeload",args,kw); }
object nodalLoad(bp::tuple args, bp::dict kw) { return CreateBCOfType("nodalload",args,kw); }
object structTemperatureLoad(bp::tuple args, bp::dict kw) { return CreateBCOfType("structtemperatureload",args,kw); }


/*****************************************************
* Material
*****************************************************/
object material(bp::tuple args, bp::dict kw)
{
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : nullptr;
    auto mat = classFactory.createMaterial(aClass.c_str(),number,domain);
    if (!mat) { OOFEM_LOG_ERROR("material: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    mat->initializeFrom(ir);
    return object(mat);
}

object CreateMaterialOfType(const char* type, bp::tuple args, bp::dict kw)
{
    args = len(args)>1? bp::make_tuple(type,args[0],args[1]) : len(args)>0? bp::make_tuple(type,args[0]) : bp::make_tuple(type);
    return material(args, kw);
}
object isoLE(bp::tuple args, bp::dict kw) { return CreateMaterialOfType("isole",args,kw); }


/*****************************************************
* CrossSection
*****************************************************/
object crossSection(bp::tuple args, bp::dict kw)
{
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : nullptr;
    auto cs = classFactory.createCrossSection(aClass.c_str(),number,domain);
    if (!cs) { OOFEM_LOG_ERROR("crossSection: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    cs->initializeFrom(ir);
    return object(cs);
}

object CreateCrossSectionOfType(const char* type, bp::tuple args, bp::dict kw)
{
    args = len(args)>1? bp::make_tuple(type,args[0],args[1]) : len(args)>0? bp::make_tuple(type,args[0]) : bp::make_tuple(type);
    return crossSection(args, kw);
}
object simpleCS(bp::tuple args, bp::dict kw) { return CreateCrossSectionOfType("simplecs",args,kw); }



/*****************************************************
* LoadTimeFunction
*****************************************************/
object loadTimeFunction(bp::tuple args, bp::dict kw)
{
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : nullptr;
    auto ltf = classFactory.createFunction(aClass.c_str(),number,domain);
    if (!ltf) { OOFEM_LOG_ERROR("function: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    ltf->initializeFrom(ir);
    return object(ltf);
}

object CreateLoadTimeFunctionOfType(const char* type, bp::tuple args, bp::dict kw)
{
    args = len(args)>1? bp::make_tuple(type,args[0],args[1]) : len(args)>0? bp::make_tuple(type,args[0]) : bp::make_tuple(type);
    return loadTimeFunction(args, kw);
}
object peakFunction(bp::tuple args, bp::dict kw) { return CreateLoadTimeFunctionOfType("peakfunction",args,kw); }



/*****************************************************
* ExportModule
*****************************************************/
object exportModule(bp::tuple args, bp::dict kw)
{
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    EngngModel *engngm = len(args)>2? extract<EngngModel*>(args[2])() : nullptr;
    auto module = classFactory.createExportModule(aClass.c_str(),number,engngm);
    if (!module) { OOFEM_LOG_ERROR("exportModule: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    module->initializeFrom(ir);
    return object(module);
}

object CreateExportModuleOfType(const char* type, bp::tuple args, bp::dict kw)
{
    args = len(args)>1? bp::make_tuple(type,args[0],args[1]) : len(args)>0? bp::make_tuple(type,args[0]) : bp::make_tuple(type);
    return exportModule(args, kw);
}
object vtkxml(bp::tuple args, bp::dict kw) { return CreateExportModuleOfType("vtkxml",args,kw); }







/*****************************************************
*
* O O F E M L I B   P Y T H O N   M O D U L E
*
*****************************************************/
BOOST_PYTHON_MODULE (liboofem)
{
    pyclass_FloatArray();
    pyclass_FloatMatrix();
    pyclass_IntArray();
    pyclass_SparseMtrx();
    // from-python converters for python tuples/arrays
    custom_AnyArray_from_sequence<FloatArray,double>();
    custom_AnyArray_from_sequence<IntArray,int>();

    converters_cxxVector_pyList_2way<int>();
    converters_cxxVector_pyList_2way<double>();
    // enums are separate types from ints
    converters_cxxVector_pyList_2way<FieldType>();

    // enumerations first, so that they can be used as default values when exposing class methods below
    pyenum_problemMode();
    pyenum_Element_Geometry_Type();
    pyenum_ValueModeType();
    pyenum_DofManTransfType();
    pyenum_DofIDItem();
    pyenum_FieldType();
    pyenum_InternalStateType();
    pyenum_MatResponseMode();
    pyenum_domainType();


    pyclass_SpatialLocalizer();
    pyclass_EngngModelContext();
    pyclass_UnknownNumberingScheme();
    pyclass_EngngModel();
    pyclass_ExportModuleManager();
    pyclass_ExportModule();
    pyclass_DataReader();
    pyclass_OOFEMTXTDataReader();
    pyclass_Domain();
    pyclass_FEMComponent();
    pyclass_DofManager();
    pyclass_Element();
    pyclass_Material();
    pyclass_StructuralMaterial();
    pyclass_StructuralMaterialSettable();
    pyclass_CrossSection();
    pyclass_GeneralBoundaryCondition();
    pyclass_BoundaryCondition();
    pyclass_Load();
    pyclass_Function();
    pyclass_MaterialStatus();
    pyclass_StructuralMaterialStatus();
    pyclass_TimeStep();
    pyclass_MetaStep();
    pyclass_DataStream();
    pyclass_Field();
    pyclass_FieldManager();
    pyclass_DofManValueField();
    pyclass_IntegrationRule();
    pyclass_GaussPoint();
    pyclass_OutputManager();
    pyclass_ClassFactory();

    /*
    Upcast converters which make boost::python accept shared_ptr<Derived> where shared_ptr<Base> is expected.
    Upcasting is automatic with boost::shared_ptr but don't work with std::shared_ptr unfortunately.
    See my post https://mail.python.org/pipermail/cplusplus-sig/2012-April/016521.html
    and the reply https://mail.python.org/pipermail/cplusplus-sig/2012-April/016533.html, and also
    boost issue https://github.com/boostorg/python/issues/29. This is a workaround which could be removed once
    boost is fixed (unlikely?).
    */
    bp::implicitly_convertible<std::shared_ptr<UniformGridField>,std::shared_ptr<Field>>();
    bp::implicitly_convertible<std::shared_ptr<UnstructuredGridField>,std::shared_ptr<Field>>();



    def("InstanciateProblem", InstanciateProblem_1, return_value_policy<manage_new_object>());
    def("createDofManValueFieldPtr", &DofManValueField_create, with_custodian_and_ward_postcall<0,2>());
    // def("FieldManager_registerField", &FieldManager_registerField);

    #ifdef OOFEM_USE_DEPRECATED_AUTOPTR
        implicitly_convertible<std::auto_ptr<DofManValueField>, std::auto_ptr<Field> >();
        register_ptr_to_python< std::auto_ptr<Field> >();
    #endif

    def("material2structuralMaterial", &material2structuralMaterial, return_internal_reference<>());
    def("material2structMatSettable", &material2structMatSettable, return_internal_reference<>());
    // modul variable classFactory
    scope().attr("classFactory") = object(ptr(&classFactory));


    // "constructors" for python interface
    def("engngModel", raw_function(engngModel,1));
    def("linearStatic", raw_function(linearStatic,0));

    def("domain", raw_function(domain,0));

    def("dofManager", raw_function(dofManager,1));
    def("node", raw_function(node,0));

    def("element", raw_function(element,1));
    def("beam2d", raw_function(beam2d,0));

    def("generalBoundaryCondition", raw_function(generalBoundaryCondition,1));
    def("boundaryCondition", raw_function(boundaryCondition,0));
    def("constantEdgeLoad", raw_function(constantEdgeLoad,0));
    def("nodalLoad", raw_function(nodalLoad,0));
    def("structTemperatureLoad", raw_function(structTemperatureLoad,0));

    def("material", raw_function(material,1));
    def("isoLE", raw_function(isoLE,0));

    def("crossSection", raw_function(crossSection,1));
    def("simpleCS", raw_function(simpleCS,0));

    def("loadTimeFunction", raw_function(loadTimeFunction,1));
    def("peakFunction", raw_function(peakFunction,0));

    def("exportModule", raw_function(exportModule,1));
    def("vtkxml", raw_function(vtkxml,0));
}

} // end namespace oofem
