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
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "engngm.h"
#include "domain.h"
#include "sparsemtrx.h"
#include "load.h"
#include "initial.h"
#include "loadtime.h"
#include "material.h"
#include "crosssection.h"
#include "node.h"
#include "conTable.h"
#include "spatiallocalizer.h"
#include "errorestimator.h"
#include "nodalrecoverymodel.h"
#include "oofemtxtdatareader.h"
#include "logger.h"
#include "util.h"
#include "problemmode.h"
#include "dofmanager.h"
#include "elementgeometrytype.h"
#include "field.h"
#include "datastream.h"
#include "dofmanvalfield.h"
#include "dofmantransftype.h"
#include "fieldmanager.h"
#include "generalbc.h"
#include "boundary.h"
#include "integrationrule.h"
#include "gausspnt.h"
#include "internalstatetype.h"
#include "matresponsemode.h"
#include "matresponseform.h"
#include "structuralmaterial.h"
#include "matstatus.h"
#include "structuralms.h"
#include "exportmodulemanager.h"
#include "outputmanager.h"
#include "classfactory.h"


namespace oofem {

/*****************************************************
* OOFEM global objects
*****************************************************/
Logger oofem_logger(Logger :: LOG_LEVEL_INFO, stdout);
Logger oofem_errLogger(Logger :: LOG_LEVEL_WARNING, stderr);
ClassFactory classFactory;





/*****************************************************
*
* O O F E M L I B   C L A S S E S
*
*****************************************************/


/*****************************************************
* FloatArray
*****************************************************/
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(flotarry_overloads_resize, resize, 1, 2)
void (FloatArray::*floatarray_add_1)(const FloatArray&) = &FloatArray::add;
void (FloatArray::*floatarray_add_2)(double, const FloatArray&) = &FloatArray::add;
void (FloatArray::*floatarray_add_3)(double) = &FloatArray::add;
double (FloatArray::*dotProduct_1)(const FloatArray&) const = &FloatArray::dotProduct;
double (FloatArray::*dotProduct_2)(const FloatArray&, int) const = &FloatArray::dotProduct;
void (FloatArray::*beDifferenceOf_1)(const FloatArray&, const FloatArray&) = &FloatArray::beDifferenceOf;
void (FloatArray::*beDifferenceOf_2)(const FloatArray&, const FloatArray&, int) = &FloatArray::beDifferenceOf;

void pyclass_FloatArray()
{
    class_<FloatArray, boost::noncopyable>("FloatArray")
        .def(init< optional<int> >())
        .def(init< FloatArray& >())
        .def("resize", &FloatArray::resize, flotarry_overloads_resize("Checks size of receiver towards requested bounds. If dimension mismatch, size is adjusted accordingly"))
        .def("containsOnlyZeroes", &FloatArray::containsOnlyZeroes, "Returns nonzero if all coefficients of the receiver are 0, else returns zero")
        .def("giveSize", &FloatArray::giveSize, "Returns the size of receiver")
        .def("isNotEmpty", &FloatArray::isNotEmpty, "Returns true if receiver is not empty")
        .def("isEmpty", &FloatArray::isEmpty, "Returns true if receiver is empty")
        .def("negated", &FloatArray::negated, "Switches the sign of every coefficient of receiver")
        .def("printYourself", &FloatArray::printYourself, "Print receiver on stdout")
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

        .def("__len__", &FloatArray::giveSize, "Returns the size of receiver")
        .def("__getitem__", &FloatArray::__getitem__, "Coefficient access function. Provides 0-based indexing access")
        .def("__setitem__", &FloatArray::__setitem__, "Coefficient access function. Provides 0-based indexing access")
        .def("beCopyOf", &FloatArray::beCopyOf, "Modifies receiver to become copy of given parameter")
        .def("copy", &FloatArray::copy, return_value_policy<return_by_value>(), "returns deep copy or receiver")
        ;
}


/*****************************************************
* FloatMatrix
*****************************************************/
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(floatmatrix_overloads_resize, resize, 2, 3)
void (FloatMatrix::*solveForRhs_1)(const FloatArray &b, FloatArray &answer, bool) = &FloatMatrix::solveForRhs;
void (FloatMatrix::*solveForRhs_2)(const FloatMatrix &b, FloatMatrix &answer, bool) = &FloatMatrix::solveForRhs;
void (FloatMatrix::*assemble_1)(const FloatMatrix&, const IntArray&) = &FloatMatrix::assemble;
void (FloatMatrix::*assemble_2)(const FloatMatrix&, const IntArray&, const IntArray&) = &FloatMatrix::assemble;
void (FloatMatrix::*assemble_3)(const FloatMatrix&, const int*, const int*) = &FloatMatrix::assemble;
void (FloatMatrix::*floatmatrix_add_1)(const FloatMatrix&) = &FloatMatrix::add;
void (FloatMatrix::*floatmatrix_add_2)(double, const FloatMatrix&) = &FloatMatrix::add;

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
        .def("beEmptyMtrx", &FloatMatrix::beEmptyMtrx, "Sets size of receiver to be an empty matrix. It will have zero rows and zero columns size")
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
        .def("resize", &FloatMatrix::resize, floatmatrix_overloads_resize("Checks size of receiver towards requested bounds. If dimension mismatch, size is adjusted accordingly"))
        .def("printYourself", &FloatMatrix::printYourself, "Prints matrix to stdout")

        .def("__setitem__", &FloatMatrix::__setitem__, "Coefficient access function. Implements 0-based indexing")
        .def("__getitem__", &FloatMatrix::__getitem__, "Coefficient access function. Implements 0-based indexing")
        .def("beCopyOf", &FloatMatrix::beCopyOf, "Modifies receiver to become copy of given parameter")
        .def("copy", &FloatMatrix::copy, return_value_policy<return_by_value>(), "returns deep copy or receiver")
        ;
}


/*****************************************************
* IntArray
*****************************************************/
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(intarray_overloads_resize, resize, 1, 2)
void pyclass_IntArray()
{
    class_<IntArray, boost::noncopyable>("IntArray")
        .def(init< optional<int> >())
        .def(init< IntArray& >())
        .def("resize", &IntArray::resize, intarray_overloads_resize("Checks size of receiver towards requested bounds. If dimension mismatch, size is adjusted accordingly"))
        .def("giveSize", &IntArray::giveSize, "Returns size of receiver")
        .def("isEmpty", &IntArray::isEmpty, "Checks if receiver is empty (i.e., zero sized)")
        .def("containsOnlyZeroes", &IntArray::containsOnlyZeroes, "Checks if receiver is all zero")
        .def("minimum", &IntArray::minimum, "Finds the minimum component in the array")
        .def("contains", &IntArray::contains, "Retuns true if receiver contains given value")
        .def("add", &IntArray::add, "Adds given scalar to all values of receiver")
        .def("zero", &IntArray::printYourself, "Sets all component to zero")
        .def("printYourself", &IntArray::printYourself, "Prints receiver on stdout")

        .def("__len__", &IntArray::giveSize, "Returns size of receiver")
        .def("__getitem__", &IntArray::__getitem__, "Coefficient access function. Provides 0-based indexing access")
        .def("__setitem__", &IntArray::__setitem__, "Coefficient access function. Provides 0-based indexing access")
        .def("beCopyOf", &IntArray::beCopyOf, "Modifies receiver to become copy of given parameter")
        .def("copy", &IntArray::copy, return_value_policy<return_by_value>(), "returns deep copy or receiver")
        ;
}


/*****************************************************
* SparseMtrx
*****************************************************/
struct PySparseMtrx : SparseMtrx, wrapper<SparseMtrx>
{
    PySparseMtrx() : SparseMtrx() {}
    PySparseMtrx(int n, int m): SparseMtrx(n,m) {}

    SparseMtrx *GiveCopy() const {
        return this->get_override("GiveCopy")();
    }

    void times(const FloatArray &x, FloatArray &answer) const {
        this->get_override("times")();
    }

    void times(double x) {
        this->get_override("times")();
    }


    int buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme& s) {
        return this->get_override("buildInternalStructure")(eModel, di, ut, s);
    }

    int buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme& r, const UnknownNumberingScheme &s) {
        if (override f = this->get_override("buildInternalStructure")) {return f(eModel, di, ut, r, s);}
        return this->get_override("buildInternalStructure")(eModel, di, ut, r, s);
    }

    int assemble(const IntArray &loc, const FloatMatrix &mat) {
        return this->get_override("assemble")();
    }

    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat) {
        return this->get_override("assemble")();
    }

    bool canBeFactorized() const {
        return this->get_override("canBeFactorized")();
    }

    void zero() {
        this->get_override("zero")();
    }

    double &at(int i, int j) {
        return this->get_override("at")();
    }

    double at(int i, int j) const {
        return this->get_override("at")();
    }

    void toFloatMatrix(FloatMatrix &answer) const {
        this->get_override("toFloatMatrix")();
    }

    void printYourself() const {
        this->get_override("printYourself")();
    }

    SparseMtrxType  giveType() const {
        return this->get_override("giveType")();
    }

    bool isAsymmetric() const {
        return this->get_override("isAsymmetric")();
    }

    FloatArray trans_mult(const FloatArray &x) const {
        return this->get_override("trans_mult")();
    }
};

void (PySparseMtrx::*sm_times_1)(const FloatArray &x, FloatArray &ans) const  = &PySparseMtrx::times;
void (PySparseMtrx::*sm_times_2)(double x) = &PySparseMtrx::times;
int (PySparseMtrx::*sm_assemble_1)(const IntArray &loc, const FloatMatrix &mat) = &PySparseMtrx::assemble;
int (PySparseMtrx::*sm_assemble_2)(const IntArray &r, const IntArray& c, const FloatMatrix &mat) = &PySparseMtrx::assemble;

int (PySparseMtrx::*sm_buildInternalStructure_1)(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme& s) = &PySparseMtrx::buildInternalStructure;
int (PySparseMtrx::*sm_buildInternalStructure_2)(EngngModel *eModel, int di, EquationID, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) = &PySparseMtrx::buildInternalStructure;

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
void pyclass_SpatialLocalizer()
{
    class_<SpatialLocalizer, boost::noncopyable>("SpatialLocalizer", no_init)
        .def("giveElementContainingPoint", &SpatialLocalizer::giveElementContainingPoint, return_internal_reference<>())
        .def("giveElementCloseToPoint", &SpatialLocalizer::giveElementCloseToPoint, return_internal_reference<>())
        .def("giveElementClosestToPoint", &SpatialLocalizer::giveElementClosestToPoint, return_internal_reference<>())
        .def("init", &SpatialLocalizer::init)
        .def("giveClassName", &SpatialLocalizer::giveClassName)
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
    //PyEngngModel(int i) : EngngModel(i, NULL) {}
    PyEngngModel(int i, EngngModel *_master = NULL) : EngngModel(i, _master) {}

    void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime) {
        this->get_override("printDofOutputAt")();
    }

    void solveYourself() {
        if (override f = this->get_override("solveYourself")) {f();return;}
        EngngModel::solveYourself();
    }
    void solveYourselfAt(TimeStep *t) {
        if (override f = this->get_override("solveYourselfAt")) {f(t); return;}
        EngngModel::solveYourselfAt(t);
    }
    void terminate(TimeStep *t) {
        if (override f = this->get_override("terminate")) {f(t); return;}
        EngngModel::terminate(t);
    }

    void updateYourself(TimeStep *t) {
        if (override f = this->get_override("updateYourself")) {f(t); return;}
        EngngModel::updateYourself(t);
    }

    TimeStep* giveNextStep() {
        if (override f = this->get_override("giveNextStep")) {return f();}
        return this->get_override("giveNextStep")();
    }

    void default_solveYourself() {return this->EngngModel::solveYourself();}
    void default_solveYourselfAt(TimeStep* t) {return this->EngngModel::solveYourselfAt(t);}
    void default_terminate(TimeStep* t) {return this->EngngModel::terminate(t);}
    void default_updateYourself(TimeStep* t) {return this->EngngModel::updateYourself(t);}
    TimeStep* default_giveNextStep() {return this->EngngModel::giveNextStep();}
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
        .def("checkProblemConsistency", &PyEngngModel::checkProblemConsistency)
        .def("setRenumberFlag", &PyEngngModel::setRenumberFlag)
        .def("giveNumberOfSteps", &EngngModel::giveNumberOfSteps)
        .add_property("numberOfSteps",&PyEngngModel::giveNumberOfSteps)
        .def("giveCurrentStep", &EngngModel::giveCurrentStep, return_internal_reference<>())
        .add_property("currentStep",make_function(&PyEngngModel::giveCurrentStep, return_internal_reference<>()))
        .def("giveNextStep",&EngngModel::giveNextStep, &PyEngngModel::default_giveNextStep, return_internal_reference<>())
        .def("giveExportModuleManager",&EngngModel::giveExportModuleManager, return_internal_reference<>())
        .add_property("exportModuleManager",make_function(&PyEngngModel::giveExportModuleManager, return_internal_reference<>()))
        .def("giveContext", &PyEngngModel::giveContext, return_internal_reference<>())
        .add_property("context", make_function(&PyEngngModel::giveContext, return_internal_reference<>()))
        ;
}

EngngModel *InstanciateProblem_1 (DataReader *dr, problemMode mode, int contextFlag)
{
    return InstanciateProblem (dr, mode, contextFlag, 0);
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
class PyExportModule : public ExportModule, public wrapper<ExportModule>
{
public:
    PyExportModule(int i, EngngModel * e) : ExportModule(i, e) {}

    void doOutput(TimeStep *tStep) {
        this->get_override("doOutput")();
    }

    void default_doOutput(TimeStep *tStep) {this->ExportModule::doOutput(tStep);}
};

void pyclass_ExportModule()
{
    class_<PyExportModule, boost::noncopyable>("ExportModule", no_init)
        .def("doOutput", pure_virtual(&ExportModule::doOutput))
        ;
}


/*****************************************************
* DataReader
*****************************************************/
struct PyDataReader : DataReader, wrapper<DataReader>
{
    InputRecord *giveInputRecord(InputRecordType irType, int recordId)
    {
        return this->get_override("giveInputRecord")();
    }

    void finish() {
        this->get_override("finish")();
    }

    const char *giveDataSourceName() const {
        return this->get_override("giveDataSourceName")();
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
        .def("giveLoadTimeFunction", &Domain::giveLoadTimeFunction, return_internal_reference<>())
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
        .def("giveNumberOfLoadTimeFunctions", &Domain::giveNumberOfLoadTimeFunctions)
        .add_property("numberOfLoadTimeFunctions",&Domain::giveNumberOfLoadTimeFunctions)
        .def("giveNumberOfRegions", &Domain::giveNumberOfRegions)
        .add_property("numberOfRegions",&Domain::giveNumberOfRegions)

        .def("resizeDofManagers", &Domain::resizeDofManagers)
        .def("resizeElements", &Domain::resizeElements)
        .def("resizeCrossSectionModels", &Domain::resizeCrossSectionModels)
        .def("resizeMaterials", &Domain::resizeMaterials)
        .def("resizeBoundaryConditions", &Domain::resizeBoundaryConditions)
        .def("resizeInitialConditions", &Domain::resizeInitialConditions)
        .def("resizeLoadTimeFunctions", &Domain::resizeLoadTimeFunctions)

        .def("setDofManager", &Domain::setDofManager)
        .def("setElement", &Domain::setElement)
        .def("setCrossSection", &Domain::setCrossSection)
        .def("setMaterial", &Domain::setMaterial)
        .def("setBoundaryCondition", &Domain::setBoundaryCondition)
        .def("setInitialCondition", &Domain::setInitialCondition)
        .def("setLoadTimeFunction", &Domain::setLoadTimeFunction)

        .def("checkConsistency", &Domain::checkConsistency)
        .def("giveConnectivityTable", &Domain::giveConnectivityTable, return_internal_reference<>())
        .def("giveSpatialLocalizer", &Domain::giveSpatialLocalizer, return_internal_reference<>())
        .def("giveErrorEstimator", &Domain::giveErrorEstimator, return_internal_reference<>())
        .def("giveSmoother", &Domain::giveSmoother, return_internal_reference<>())

        .def("giveNumberOfSpatialDimensions", &Domain::giveNumberOfSpatialDimensions)
        .add_property("numberOfSpatialDimensions", &Domain::giveNumberOfSpatialDimensions)
        .def("giveDomainType", &Domain::giveDomainType)
        .def("setDomainType", &Domain::setDomainType)
        .add_property("domainType", &Domain::giveDomainType, &Domain::setDomainType)

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

    void giveUnknownVector(FloatArray &answer, const IntArray &dofMask,
                EquationID type, ValueModeType mode, TimeStep *stepN) {
        if ( override f = this->get_override("giveUnknownVector") ) {
        f (answer, dofMask, type, mode, stepN); return;}
        DofManager::giveUnknownVector(answer, dofMask, type, mode, stepN);
    }
    bool computeL2GTransformation(FloatMatrix &answer, const IntArray &dofIDArry) {
        if (override f = this->get_override("computeL2GTransformation")) { return f(answer, dofIDArry); }
        return DofManager::computeL2GTransformation(answer, dofIDArry);
    }

    void default_giveUnknownVector(FloatArray &answer, const IntArray &dofMask,
                    EquationID type, ValueModeType mode, TimeStep *stepN)
    { return this->DofManager::giveUnknownVector(answer, dofMask, type, mode, stepN); }
};

void (DofManager::*giveUnknownVector_1)(FloatArray &answer, const IntArray &dofMask,
                    EquationID type, ValueModeType mode, TimeStep *stepN) = &DofManager::giveUnknownVector;

void pyclass_DofManager()
{
    class_<DofManager, bases<FEMComponent>, boost::noncopyable>("DofManager", no_init)
        .def("hasCoordinates", &DofManager::hasCoordinates)
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

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        this->get_override("giveCharacteristicMatrix")();
    }

    void  giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) {
        this->get_override("giveCharacteristicVector")();
    }

    double giveCharacteristicValue(CharType, TimeStep *) {
        return this->get_override("giveCharacteristicValue")();
    }

    void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const {
        this->get_override("giveDofManDofIDMask")();
    }
};

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
        .def("giveLocationArray", &PyElement::giveLocationArray)
        .def("invalidateLocationArray", &PyElement::invalidateLocationArray)
        .def("giveNumberOfDofManagers", &PyElement::giveNumberOfDofManagers)
        .add_property("numberOfDofManagers", &PyElement::giveNumberOfDofManagers)
        .def("computeNumberOfDofs", &PyElement::computeNumberOfDofs)
        .def("giveNumberOfNodes", &PyElement::giveNumberOfNodes)
        .def("giveGeometryType", &PyElement::giveGeometryType)
        .add_property("geometryType", &PyElement::giveGeometryType)
        .def("giveDofManagerNumber", &Element::giveDofManagerNumber)
        .def("giveDofManager", &Element::giveDofManager, return_internal_reference<>())
        .def("giveNode", &Element::giveNode, return_internal_reference<>())
        // TODO return type (copy rather than pointer?)
        .def("giveDofManArray", &Element::giveDofManager, return_value_policy<reference_existing_object>())
        .def("giveNumberOfIntegrationRules", &Element::giveNumberOfIntegrationRules)
        .add_property("numberOfIntegrationRules", &PyElement::giveNumberOfIntegrationRules)
        .def("giveDefaultIntegrationRulePtr", &Element::giveDefaultIntegrationRulePtr, return_internal_reference<>())
        .add_property("defaultIntegrationRule",make_function(&PyElement::giveDefaultIntegrationRulePtr, return_internal_reference<>()))
        .def("giveMaterial", &Element::giveMaterial, return_internal_reference<>())
        .def("setMaterial", &Element::setMaterial)
        .add_property("material",make_function(&PyElement::giveMaterial, return_internal_reference<>()), &PyElement::setMaterial)
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
        .def("giveCharacteristicMatrix", &Material::giveCharacteristicMatrix)
        .def("giveStatus", &Material::giveStatus, return_internal_reference<>())
        ;
}


/*****************************************************
* StructuralMaterial
*****************************************************/
struct PyStructuralMaterial : StructuralMaterial , wrapper<StructuralMaterial>
{
    PyStructuralMaterial(int i, Domain *d) : StructuralMaterial(i,d) {}
    void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) {
        if (override f = this->get_override("giveRealStressVector")) { f(answer,form,gp,reducedStrain,tStep);}
        this->get_override("giveRealStressVector")(answer,form,gp,reducedStrain,tStep);
    }
};
void pyclass_StructuralMaterial()
{
    class_<PyStructuralMaterial, bases<Material>, boost::noncopyable >("StructuralMaterial", no_init)
        .def("giveRealStressVector", pure_virtual( &StructuralMaterial::giveRealStressVector))
        ;
}

StructuralMaterial* material2structuralMaterial(Material *mat) { return (StructuralMaterial*)mat; }



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
        .def("give", &BoundaryCondition::give)
        .def("isImposed", &BoundaryCondition::isImposed)
        ;
}


/*****************************************************
* Load
*****************************************************/
struct PyLoad : Load , wrapper<Load>
{
    PyLoad(int i, Domain *d) : Load(i,d) {}
    void computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode) {
        if (override f = this->get_override("computeValueAt")) { f(answer,tStep,coords,mode);}
        this->get_override("computeValueAt")(answer,tStep,coords,mode);
    }
};
void pyclass_Load()
{
    class_<PyLoad, bases<GeneralBoundaryCondition>, boost::noncopyable >("Load", init<int, Domain*>())
        .def("setComponentArray", &Load::setComponentArray)
        .def("giveCopyOfComponentArray", &Load::giveCopyOfComponentArray)
        .def("computeValueAt", pure_virtual( &Load::computeValueAt))
        ;
}


/*****************************************************
* LoadTimeFunction
*****************************************************/
void pyclass_LoadTimeFunction()
{
    class_<LoadTimeFunction, bases<FEMComponent>, boost::noncopyable >("LoadTimeFunction", no_init)
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
    class_<TimeStep, boost::noncopyable>("TimeStep", no_init)
        .def("giveTargetTime", &TimeStep::giveTargetTime)
        .def("setTargetTime", &TimeStep::setTargetTime)
        .add_property("targetTime",&TimeStep::giveTargetTime, &TimeStep::setTargetTime)
        .def("giveIntrinsicTime", &TimeStep::giveIntrinsicTime)
        .def("setIntrinsicTime", &TimeStep::setIntrinsicTime)
        .add_property("intrinsicTime",&TimeStep::giveIntrinsicTime, &TimeStep::setIntrinsicTime)
        .def("giveTimeIncrement", &TimeStep::giveTimeIncrement)
        .def("setTimeIncrement", &TimeStep::setTimeIncrement)
        .add_property("timeIncrement",&TimeStep::giveTimeIncrement, &TimeStep::setTimeIncrement)
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
class PyField : public Field, public wrapper<Field>
{
public:
    PyField (FieldType b): Field(b) {}
    int evaluateAt(FloatArray &answer, FloatArray &coords,
            ValueModeType mode, TimeStep *atTime) {
        return this->get_override("evaluateAt")(answer, coords,mode,atTime);
    }
    int evaluateAt(FloatArray &answer, DofManager* dman,
            ValueModeType mode, TimeStep *atTime) {
        return this->get_override("evaluateAt")(answer,dman,mode,atTime);
    }
    contextIOResultType saveContext(DataStream *stream, ContextMode mode) {
        return this->get_override("saveContext")(stream, mode);
    }
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode) {
        return this->get_override("restoreContext")(stream, mode);
    }
};

int (PyField::*evaluateAtPos)(FloatArray &answer, FloatArray &coords, ValueModeType mode, TimeStep *atTime) = &PyField::evaluateAt;
int (PyField::*evaluateAtDman)(FloatArray &answer, DofManager* dman, ValueModeType mode, TimeStep *atTime) = &PyField::evaluateAt;

std::auto_ptr<Field> DofManValueField_create (FieldType t, Domain* d) { return std::auto_ptr<Field> ( new DofManValueField(t, d) ); }

void pyclass_Field()
{
    //this class is held by auto_ptr:  to allow for taking ownership of a raw pointer
    //see http://www.boost.org/doc/libs/1_47_0/libs/python/doc/v2/faq.html#ownership for detail
    // also see http://stackoverflow.com/questions/4112561/boost-python-ownership-of-pointer-variables
    class_<Field, PyField, boost::noncopyable>("Field", no_init)
        .def("evaluateAtPos", evaluateAtPos)
        .def("evaluateAtDman", evaluateAtDman)
        .def("giveType", &PyField::giveType)
        ;
}


/*****************************************************
* FieldManager
*****************************************************/
class PyFieldManager;
//void FieldManager_registerField (FieldManager* fm, std::auto_ptr<Field>& a, FieldType key, bool managedFlag)
/*
void FieldManager_registerField (PyFieldManager* fm, DofManValueField* a, FieldType key, bool managedFlag)
{
    //fm->registerField(a.get(), key, managedFlag);
    //a.release();
    printf("Registering field\n");
};
*/
//void FieldManager_registerField (FieldManager* fm, Field* a, FieldType key, bool managedFlag ) // ok
void FieldManager_registerField (FieldManager* fm, std::auto_ptr<Field>& a, FieldType key, bool managedFlag ) // ok
{
    fm->registerField(a.get(), key, managedFlag);
    a.release();
    printf("Registering field\n");
};

class PyFieldManager : public FieldManager, public wrapper<FieldManager>
{
public:
    PyFieldManager (): FieldManager() {}
    //void myRegisterField(std::auto_ptr<Field>& a, FieldType key, bool managedFlag = false) {
    void myRegisterField(DofManValueField* a, FieldType key, bool managedFlag = false) {
        //FieldManager_registerField(this, a, key, managedFlag);
    }
    //void registerField(std::auto_ptr<PyField> a, FieldType key, bool managedFlag = false) {
    //  FieldManager_registerField(*this, a, key, managedFlag);
    //}

};

void pyclass_FieldManager()
{
    class_<FieldManager, PyFieldManager, boost::noncopyable>("FieldManager", no_init)
        .def("registerField", &PyFieldManager::myRegisterField)
        .def("isFieldRegistered", &FieldManager::isFieldRegistered)
        .def("giveField", &FieldManager::giveField, return_internal_reference<>())
        .def("unregisterField", &FieldManager::unregisterField)
        ;
}


/*****************************************************
* DofManValueField
*****************************************************/
void pyclass_DofManValueField()
{
    class_<DofManValueField, bases<Field>, boost::noncopyable >("DofManValueField", no_init)
        .def("evaluateAt", evaluateAtPos)
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
        .def("getNumberOfIntegrationPoints", &IntegrationRule::getNumberOfIntegrationPoints)
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
* EquationID
*****************************************************/
void pyenum_EquationID()
{
    enum_<EquationID>("EquationID")
        .value("EID_MomentumBalance", EID_MomentumBalance)
        .value("EID_AuxMomentumBalance", EID_AuxMomentumBalance)
        .value("EID_ConservationEquation", EID_ConservationEquation)
        .value("EID_MomentumBalance_ConservationEquation", EID_MomentumBalance_ConservationEquation)
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
* MatResponseFrom
*****************************************************/
void pyenum_MatResponseForm()
{
    enum_<MatResponseForm>("MatResponseForm")
        .value("ReducedForm", ReducedForm)
        .value("FullForm", FullForm)
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

// constuct auxiliary global object (for eval and exec functions) only once
dict temp_global(import("__main__").attr("__dict__"));

/*
make oofem input line from **kw arguments. This function is used by "constructor" methods
e.g. in function:
  engngModel("nonLinearStatic",nSteps=10) # kw = {'nSteps':10}
*/
OOFEMTXTInputRecord makeOOFEMTXTInputRecordFrom(dict &kw) {
    temp_global["kw"] = kw;
    str command =
        "ret = ''\n"
        "for key,val in kw.iteritems():\n" // iterate over key,val pairs
        "  if key.lower()=='number' or key.lower()=='domain': continue\n" // do not include "number" and "domain" kws
        "  if key=='f_t': key='f(t)'\n" // handle f(t) loadTimeFunction field name
        "  ret += ' %s'%key\n" // add key
        "  if isinstance(val,(int,float,str)): ret += ' %s'%val\n" // add val if it is int, float or str
        "  elif isinstance(val,(list,tuple)):\n" // if val is tuple or list
        "    ret += ' %d'%len(val)\n" // add its length first
        "    for v in val:\n" // and then its values
        "      if isinstance(v,(int,float,str)): ret += ' %s'%v\n" // int, float, str
        "      else: ret += ' %d'%(v.giveNumber())\n" // or any other object, that have giveNumber method (e.g. nodes can be passed to elemnt function as they are, without the need of extracting their numbers first
        "  else: ret += ' %d'%(val.giveNumber())\n" // add arbitrary object with giveNumber method
        "ret = ret.lower()\n" // finally make it lower case
        "print ret\n"
        ;
    exec(command,temp_global,temp_global);
    // extract string from globals["ret"], convert it to char* and return OOFEMTXTInputRecord from it
    return OOFEMTXTInputRecord( ( extract<string>(temp_global["ret"])() ).c_str() );
}

OOFEMTXTInputRecord makeOutputManagerOOFEMTXTInputRecordFrom(dict &kw) {
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
void makeDictKeysLowerCase(dict &kw) {
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
object engngModel(tuple args, dict kw) {
    //args
    string aClass = extract<string>(args[0])();
    int number = len(args)>1? extract<int>(args[1])() : 0;
    EngngModel* master = len(args)>2? extract<EngngModel*>(args[2])() : NULL;
    EngngModel *engngm = classFactory.createEngngModel(aClass.c_str(),number,master);
    if (engngm==NULL) { LOG_ERROR(oofem_errLogger,"engngModel: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    engngm->initializeFrom(&ir);
    // instanciateYourself
    if ( ir.hasField(IFT_EngngModel_nmsteps, "nmsteps") ) {
        LOG_ERROR(oofem_errLogger,"engngModel: simulation with metasteps is not (yet) supported in Python");
    } else {
        engngm->instanciateDefaultMetaStep(&ir);
    }
    string outFile;
    if ( ir.hasField(IFT_EngngModel_nmsteps, "outfile") ) {
       ir.giveField(outFile, IFT_EngngModel_outfile, "outfile");
    } else {
       outFile = "oofem.out.XXXXXX";
    }
    engngm->Instanciate_init(outFile.c_str(), engngm->giveNumberOfDomains());
    //
    object ret = object(ptr(engngm));
    /* ????????????????????
    // sets the last created engngModel as default one for further script
    temp_global["defaultEngngModel"] = ret;
    ???????????????????? */
    return ret;
}

object createEngngModelOfType(const char* type, tuple args, dict kw) {
    args = len(args)>1? make_tuple(type,args[0],args[1]) : len(args)>0? make_tuple(type,args[0]) : make_tuple(type);
    return engngModel(args,kw);
}
object linearStatic(tuple args, dict kw) { return createEngngModelOfType("linearstatic", args, kw); }


/*****************************************************
* Domain
*****************************************************/
// domain(number=1,serialNumber=1,engngModel*=None,dType=_unknownMode,**kw)
object domain(tuple args, dict kw) {
    // args
    int number =             len(args)>0? extract<int>(args[0])() : 0;
    int serialNumber =       len(args)>1? extract<int>(args[1])() : 0;
    EngngModel *engngModel = len(args)>2? extract<EngngModel*>(args[2])() : NULL;
    domainType dType =       len(args)>3? extract<domainType>(args[3])() : _unknownMode;
    Domain *d = new Domain(number,serialNumber,engngModel);
    d->setDomainType(dType);
    // output manager record
    OOFEMTXTInputRecord omir = makeOutputManagerOOFEMTXTInputRecordFrom(kw);
    d->giveOutputManager()->initializeFrom(&omir);
    object ret = object(ptr(d));
    /* ????????????????????
    // sets the last created domain as default one for furtherscript
    temp_global["defaultDomain"] = ret;
    ???????????????????? */
    return object(ptr(d));
}


/*****************************************************
* Element
*****************************************************/
// element(aClass,domain=defaultDomain,**kw)
object element(tuple args, dict kw) {
    // extracts first python argument (string element type)
    string aClass = extract<string>(args[0])();
    // extracts values from args if they are specified
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : NULL;
    /* ????????????????
    // if no material is specified, set it to 1
    if (!kw.has_key("mat")) { kw["mat"] = 1; }
    // if no cross section is specified, set it to 1
    if (!kw.has_key("crosssect")) { kw["crosssect"] = 1; }
    // if no domain is specified and one already exists in the script, use that one
    if (domain == NULL && temp_global.has_key("defaultDomain")) { domain = extract<Domain*>(temp_global["defaultDomain"])(); }
    if (domain==NULL) { LOG_ERROR(oofem_errLogger,"wrong Domain"); }
    ???????????????????? */
    // create Element (convert aClass to char* - expected by classFactory.createElement)
    Element *elem = classFactory.createElement(aClass.c_str(),number,domain);
    // if elem==NULL, something was wrong
    if (elem==NULL) { LOG_ERROR(oofem_errLogger,"element: wrong input data"); }
    // sets globalNumber == number befor initializeFrom
    elem->setGlobalNumber(number);
    // construct OOFEMTXTInputRecord from dict **kw
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    // pass input record to elem
    elem->initializeFrom(&ir);
    // convert element to PyObject (expected by raw_function, which enables fun(*args,**kw) syntax in python)
    return object(ptr(elem));
}

// auxiliary constructor for specific element type (name is passed as first argument)
object createElementOfType(const char* type, tuple args, dict kw) {
    args = len(args)>1? make_tuple(type,args[0],args[1]) : len(args)>0? make_tuple(type,args[0]) : make_tuple(type);
    return element(args, kw);
}
// specific elements
object beam2d(tuple args, dict kw) { return createElementOfType("beam2d",args,kw); }



/*****************************************************
* DofManager
*****************************************************/
// dofManager(aClass,domain=defaultDomain,**kw)
object dofManager(tuple args, dict kw) {
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : NULL;
    DofManager *dofMan = classFactory.createDofManager(aClass.c_str(),number,domain);
    if (dofMan==NULL) { LOG_ERROR(oofem_errLogger,"dofManager: wrong input data"); }
    dofMan->setGlobalNumber(number);
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    dofMan->initializeFrom(&ir);
    return object(ptr(dofMan));
}

object createDofManagerOfType(const char* type, tuple args, dict kw) {
    args = len(args)>1? make_tuple(type,args[0],args[1]) : len(args)>0? make_tuple(type,args[0]) : make_tuple(type);
    return dofManager(args, kw);
}
object node(tuple args, dict kw) { return createDofManagerOfType("node",args,kw); }


/*****************************************************
* GeneralBoundaryCondition
*****************************************************/
// generalBoundaryCondition(aClass,domain=defaultDomain,**kw)
object generalBoundaryCondition(tuple args, dict kw) {
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : NULL;
    GeneralBoundaryCondition *bc = classFactory.createBoundaryCondition(aClass.c_str(),number,domain);
    if (bc==NULL) { LOG_ERROR(oofem_errLogger,"generalBoundaryCondition: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    bc->initializeFrom(&ir);
    return object(ptr(bc));
}

object CreateBCOfType(const char* type, tuple args, dict kw) {
    args = len(args)>1? make_tuple(type,args[0],args[1]) : len(args)>0? make_tuple(type,args[0]) : make_tuple(type);
    return generalBoundaryCondition(args, kw);
}
object boundaryCondition(tuple args, dict kw) { return CreateBCOfType("boundarycondition",args,kw); }
object constantEdgeLoad(tuple args, dict kw) { return CreateBCOfType("constantedgeload",args,kw); }
object nodalLoad(tuple args, dict kw) { return CreateBCOfType("nodalload",args,kw); }
object structTemperatureLoad(tuple args, dict kw) { return CreateBCOfType("structtemperatureload",args,kw); }


/*****************************************************
* Material
*****************************************************/
object material(tuple args, dict kw) {
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : NULL;
    Material *mat = classFactory.createMaterial(aClass.c_str(),number,domain);
    if (mat==NULL) { LOG_ERROR(oofem_errLogger,"material: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    mat->initializeFrom(&ir);
    return object(ptr(mat));
}

object CreateMaterialOfType(const char* type, tuple args, dict kw) {
    args = len(args)>1? make_tuple(type,args[0],args[1]) : len(args)>0? make_tuple(type,args[0]) : make_tuple(type);
    return material(args, kw);
}
object isoLE(tuple args, dict kw) { return CreateMaterialOfType("isole",args,kw); }


/*****************************************************
* CrossSection
*****************************************************/
object crossSection(tuple args, dict kw) {
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : NULL;
    CrossSection *cs = classFactory.createCrossSection(aClass.c_str(),number,domain);
    if (cs==NULL) { LOG_ERROR(oofem_errLogger,"crossSection: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    cs->initializeFrom(&ir);
    return object(ptr(cs));
}

object CreateCrossSectionOfType(const char* type, tuple args, dict kw) {
    args = len(args)>1? make_tuple(type,args[0],args[1]) : len(args)>0? make_tuple(type,args[0]) : make_tuple(type);
    return crossSection(args, kw);
}
object simpleCS(tuple args, dict kw) { return CreateCrossSectionOfType("simplecs",args,kw); }



/*****************************************************
* LoadTimeFunction
*****************************************************/
object loadTimeFunction(tuple args, dict kw) {
    string aClass = extract<string>(args[0])();
    int number =     len(args)>1? extract<int>(args[1])() : 0;
    Domain *domain = len(args)>2? extract<Domain*>(args[2])() : NULL;
    LoadTimeFunction *ltf = classFactory.createLoadTimeFunction(aClass.c_str(),number,domain);
    if (ltf==NULL) { LOG_ERROR(oofem_errLogger,"loadTimeFunction: wrong input data"); }
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    ltf->initializeFrom(&ir);
    return object(ptr(ltf));
}

object CreateLoadTimeFunctionOfType(const char* type, tuple args, dict kw) {
    args = len(args)>1? make_tuple(type,args[0],args[1]) : len(args)>0? make_tuple(type,args[0]) : make_tuple(type);
    return loadTimeFunction(args, kw);
}
object peakFunction(tuple args, dict kw) { return CreateLoadTimeFunctionOfType("peakfunction",args,kw); }







/*****************************************************
*
* O O F E M L I B   P Y T H O N   M O D U L E
*
*****************************************************/
BOOST_PYTHON_MODULE (oofemlib)
{
    pyclass_FloatArray();
    pyclass_FloatMatrix();
    pyclass_IntArray();
    pyclass_SparseMtrx();
    pyclass_SpatialLocalizer();
    pyclass_EngngModelContext();
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
    pyclass_CrossSection();
    pyclass_GeneralBoundaryCondition();
    pyclass_BoundaryCondition();
    pyclass_Load();
    pyclass_LoadTimeFunction();
    pyclass_MaterialStatus();
    pyclass_StructuralMaterialStatus();
    pyclass_TimeStep();
    pyclass_DataStream();
    pyclass_Field();
    pyclass_FieldManager();
    pyclass_DofManValueField();
    pyclass_IntegrationRule();
    pyclass_GaussPoint();
    pyclass_OutputManager();
    pyclass_ClassFactory();

    pyenum_problemMode();
    pyenum_Element_Geometry_Type();
    pyenum_EquationID();
    pyenum_ValueModeType();
    pyenum_DofManTransfType();
    pyenum_DofIDItem();
    pyenum_FieldType();
    pyenum_InternalStateType();
    pyenum_MatResponseMode();
    pyenum_MatResponseForm();
    pyenum_domainType();


    //def("make_foo", make_foo, return_value_policy<manage_new_object>())
    def("InstanciateProblem", InstanciateProblem_1, return_value_policy<manage_new_object>());
    def("createDofManValueFieldPtr", &DofManValueField_create, with_custodian_and_ward_postcall<0,2>());
    def("FieldManager_registerField", &FieldManager_registerField);

    implicitly_convertible<std::auto_ptr<DofManValueField>, std::auto_ptr<Field> >();
    register_ptr_to_python< std::auto_ptr<Field> >();

    def("material2structuralMaterial", &material2structuralMaterial, return_internal_reference<>());
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

}
} // end namespace oofem
