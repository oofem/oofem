#include <boost/python.hpp>
//#include <boost/python/return_internal_reference.hpp>
//#include <boost/python/manage_new_object.hpp>
//#include <boost/python/return_value_policy.hpp>
using namespace boost::python;
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

namespace oofem {
 
/* DefaulT oofem loggers */
Logger oofem_logger(Logger :: LOG_LEVEL_INFO, stdout);
Logger oofem_errLogger(Logger :: LOG_LEVEL_WARNING, stderr);

EngngModel *InstanciateProblem_1 (DataReader *dr, problemMode mode, int contextFlag) {
  return oofem::InstanciateProblem (dr, mode, contextFlag, 0);
}
 
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(flotarry_overloads_resize, resize, 1, 2)

void (FloatArray::*add_1)(const FloatArray &src) = &FloatArray::add;
void (FloatArray::*subtract_1)(const FloatArray &src) = &FloatArray::subtract;

class PyFloatMatrix : public FloatMatrix, public wrapper<FloatMatrix>
 {
public:
   PyFloatMatrix () : FloatMatrix () {}
   PyFloatMatrix (int a, int b) : FloatMatrix (a,b) {}
   void __setitem__ (object t, double val) {
    this->at(extract<int>(t[0]), extract<int>(t[1]) ) = val;}
   double __getitem__ (object t, double val) {
     return this->at(extract<int>(t[0]), extract<int>(t[1]) );}
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(floatmatrix_overloads_resize, resize, 2, 3)

void (PyFloatMatrix::*solveForRhs_1)(const FloatArray &b, FloatArray &answer) = &PyFloatMatrix::solveForRhs;
void (PyFloatMatrix::*solveForRhs_2)(const FloatMatrix &b, FloatMatrix &answer) = &PyFloatMatrix::solveForRhs;
void (PyFloatMatrix::*assemble_1)(const FloatMatrix &src, const IntArray &loc) = &PyFloatMatrix::assemble;
void (PyFloatMatrix::*assemble_2)(const FloatMatrix &src, const IntArray &row, const IntArray&col) = &PyFloatMatrix::assemble;


BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(intarray_overloads_resize, resize, 1, 2)


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

  int buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &s) {
    return this->get_override("buildInternalStructure")();
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


class PyEngngModel : public EngngModel, public wrapper<EngngModel>
{
public:
  //PyEngngModel(int i) : EngngModel(i, NULL) {}
  PyEngngModel(int i, EngngModel *_master = NULL) : EngngModel(i, _master) {}
  PyEngngModel(int i, char *s, EngngModel *_master = NULL): EngngModel(i,s,_master) {}

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

    
  void default_solveYourself() {return this->EngngModel::solveYourself();}
  void default_solveYourselfAt(TimeStep*t) {return this->EngngModel::solveYourselfAt(t);}
  void default_terminate(TimeStep* t) {return this->EngngModel::terminate(t);}
  void default_updateYourself(TimeStep* t) {return this->EngngModel::updateYourself(t);}
};



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


BOOST_PYTHON_MODULE (oofemlib)
{

  class_<FloatArray>("FloatArray")
    .def(init< optional<int> >())
    //.def("set", at2, return_internal_reference<>())
    .def("__len__", &FloatArray::giveSize)
    .def("__getitem__", &FloatArray::__getItem__) 
    .def("__setitem__", &FloatArray::__setItem__)
    .def("resize", &FloatArray::resize, flotarry_overloads_resize())
    .def("giveSize", &FloatArray::giveSize)
    .def("isNotEmpty", &FloatArray::isNotEmpty)
    .def("isEmpty", &FloatArray::isEmpty)  
    //.def("negated", &FloatArray::negated)
    .def("printYourself", &FloatArray::printYourself)
    .def("zero", &FloatArray::zero)
    .def("add", add_1)
    .def("subtract", subtract_1)
    //.def("times", &FloatArray::times)
    //.def("normalize", &FloatArray::normalize)
    .def("computeNorm", &FloatArray::computeNorm)
    .def("beProductOf", &FloatArray::beProductOf)
    .def("beTProductOf", &FloatArray::beTProductOf)
    //.def("beCopyOf", &FloatArray::beCopyOf)
    .def("assemble", &FloatArray::assemble)
    ;


  class_<PyFloatMatrix, boost::noncopyable>("FloatMatrix")
    .def(init<int, int>())
    .def("zero", &PyFloatMatrix::zero)
    .def("beUnitMatrix", &PyFloatMatrix::beUnitMatrix)
    .def("printYourself", &PyFloatMatrix::printYourself)
    .def("__setitem__", &PyFloatMatrix::__setitem__)
    .def("__getitem__", &PyFloatMatrix::__getitem__)
    .def("giveDeterminant", &PyFloatMatrix::giveDeterminant)
    .def("beTranspositionOf", &PyFloatMatrix::beTranspositionOf)
    .def("beProductOf", &PyFloatMatrix::beProductOf)
    .def("beTProductOf", &PyFloatMatrix::beTProductOf)
    .def("beProductTOf", &PyFloatMatrix::beProductTOf)
    .def("beInverseOf", &PyFloatMatrix::beInverseOf)
    .def("solveForRhs", solveForRhs_1)
    .def("plusProductSymmUpper", &PyFloatMatrix::plusProductSymmUpper)
    .def("plusProductUnsym", &PyFloatMatrix::plusProductUnsym)
    .def("add", &PyFloatMatrix::add)
    .def("symmetrized", &PyFloatMatrix::symmetrized)
    .def("rotatedWith", &PyFloatMatrix::rotatedWith)
    .def("resize", &PyFloatMatrix::resize, floatmatrix_overloads_resize())
    .def("assemble", assemble_1)
    .def("assembleu", assemble_2)
    ;

  class_<IntArray>("IntArray")
    .def(init< optional<int> >())
    .def("__len__", &FloatArray::giveSize)
    .def("__getitem__", &IntArray::__getItem__) 
    .def("__setitem__", &IntArray::__setItem__)
    .def("resize", &IntArray::resize, intarray_overloads_resize())
    .def("isEmpty", &IntArray::isEmpty)
    .def("containsOnlyZeroes", &IntArray::containsOnlyZeroes)
    .def("contains", &IntArray::contains)
    .def("printYourself", &IntArray::printYourself)
    .def("zero", &IntArray::printYourself)
    ;

  class_<PySparseMtrx, boost::noncopyable>("SparseMtrx")
    //.def(init<>())
    .def("times", pure_virtual(sm_times_1))
    .def("timesd", pure_virtual(sm_times_2))
    .def("buildInternalStructure", &PySparseMtrx::buildInternalStructure)
    .def("assemble", sm_assemble_1)
    .def("assembleu", sm_assemble_2)
    //.def("zero", zero)
    ;
  
  class_<PyEngngModel, boost::noncopyable>("EngngModel", no_init)
    //.def(init<int, optional<EngngModel*> > ())
    //.def(init<int, char*, optional<EngngModel*> >())
    .def("giveDomain", &PyEngngModel::giveDomain, return_internal_reference<>())
    .def("giveNumberOfDomains", &PyEngngModel::giveNumberOfDomains)
    .def("giveDomainErrorEstimator", &PyEngngModel::giveDomainErrorEstimator, return_internal_reference<>())
    .def("terminateAnalysis", &PyEngngModel::terminateAnalysis)
    .def("solveYourself", &EngngModel::solveYourself, &PyEngngModel::default_solveYourself)
    .def("solveYourselfAt", &EngngModel::solveYourselfAt, &PyEngngModel::default_solveYourselfAt)
    .def("terminate", &EngngModel::terminate, &PyEngngModel::default_terminate)
    .def("updateYourself", &EngngModel::updateYourself, &PyEngngModel::default_updateYourself)
    .def("initializeYourself", &PyEngngModel::initializeYourself)
    .def("printDofOutputAt", pure_virtual(&EngngModel::printDofOutputAt))
    .def("checkProblemConsistency", &PyEngngModel::checkProblemConsistency)
    .def("setRenumberFlag", &PyEngngModel::setRenumberFlag)
    ;

  class_<Domain>("Domain", init<int, int, EngngModel* >())
  .def("giveNumber", &Domain::giveNumber)
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
  
  .def("giveNumberOfDofManagers", &Domain::giveNumberOfDofManagers)
  .def("giveNumberOfElements", &Domain::giveNumberOfElements)
  .def("giveNumberOfMaterialModels", &Domain::giveNumberOfMaterialModels)
  .def("giveNumberOfCrossSectionModels", &Domain::giveNumberOfCrossSectionModels)
  .def("giveNumberOfBoundaryConditions", &Domain::giveNumberOfBoundaryConditions)
  .def("giveNumberOfInitialConditions", &Domain::giveNumberOfInitialConditions)
  .def("giveNumberOfLoadTimeFunctions", &Domain::giveNumberOfLoadTimeFunctions)
  
  .def("checkConsistency", &Domain::checkConsistency)
  .def("giveConnectivityTable", &Domain::giveConnectivityTable, return_internal_reference<>())
  .def("giveSpatialLocalizer", &Domain::giveSpatialLocalizer, return_internal_reference<>())
  .def("giveErrorEstimator", &Domain::giveErrorEstimator, return_internal_reference<>())
  .def("giveSmoother", &Domain::giveSmoother, return_internal_reference<>())
  ;
  
  
  class_<PyDataReader, boost::noncopyable>("DataReader")
    ;
  class_<OOFEMTXTDataReader, bases<DataReader> >("OOFEMTXTDataReader", init<char* >())
    ;
  enum_<problemMode>("problemMode")
    .value("_processor", _processor)
    .value("_postProcessor", _postProcessor)
    ;

  def("InstanciateProblem", InstanciateProblem_1, return_value_policy<manage_new_object>());
  //def("make_foo", make_foo, return_value_policy<manage_new_object>())
    
  class_<PyElement, boost::noncopyable>("Element", init<int, Domain* >())
    .def("giveLocationArray", &PyElement::giveLocationArray)
    .def("invalidateLocationArray", &PyElement::invalidateLocationArray)
    .def("giveCharacteristicMatrix", &PyElement::giveCharacteristicMatrix)
    .def("giveCharacteristicVector", &PyElement::giveCharacteristicVector)
    .def("giveCharacteristicValue", &PyElement::giveCharacteristicValue)
    .def("computeNumberOfDofs", &PyElement::computeNumberOfDofs)
    .def("giveGeometryType", &PyElement::giveGeometryType)
    ;


}
} // end namespace oofem
