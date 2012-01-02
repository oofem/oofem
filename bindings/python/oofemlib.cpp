#include <boost/python.hpp>
//#include <boost/python/return_internal_reference.hpp>
//#include <boost/python/manage_new_object.hpp>
//#include <boost/python/return_value_policy.hpp>
using namespace boost::python;
using namespace std;

#include <memory> // for std::auto_ptr<>
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
#include "fieldmanager.h"
#include "generalbc.h"
#include "boundary.h"
#include "integrationrule.h"
#include "gausspnt.h"
#include "internalstatetype.h"
#include "matresponsemode.h"
#include "matresponseform.h"

 
namespace oofem {

/* DefaulT oofem loggers */
Logger oofem_logger(Logger :: LOG_LEVEL_INFO, stdout);
Logger oofem_errLogger(Logger :: LOG_LEVEL_WARNING, stderr);

EngngModel *InstanciateProblem_1 (DataReader *dr, problemMode mode, int contextFlag) {
  return InstanciateProblem (dr, mode, contextFlag, 0);
}

/*****************************************************
* FloatArray
*****************************************************/
class PyFloatArray : public FloatArray , public wrapper<FloatArray>
{
public:
   PyFloatArray () : FloatArray () {}
   PyFloatArray (int i) : FloatArray(i) {}
   void __setitem__(int i, double val) { this->operator()(i) = val; }
   double __getitem__(int i) { return this->operator()(i); }
};
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(flotarry_overloads_resize, resize, 1, 2)
void (FloatArray::*add_1)(const FloatArray &src) = &FloatArray::add;
void (FloatArray::*subtract_1)(const FloatArray &src) = &FloatArray::subtract;

void pyclass_FloatArray()
{
  class_<PyFloatArray, boost::noncopyable>("FloatArray")
    .def(init< optional<int> >())
    //.def("set", at2, return_internal_reference<>())
    .def("__len__", &FloatArray::giveSize)
    .def("__getitem__", &PyFloatArray::__getitem__) 
    .def("__setitem__", &PyFloatArray::__setitem__)
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
}


/*****************************************************
* FloatMatrix
*****************************************************/
class PyFloatMatrix : public FloatMatrix, public wrapper<FloatMatrix>
 {
public:
   PyFloatMatrix () : FloatMatrix () {}
   PyFloatMatrix (int a, int b) : FloatMatrix (a,b) {}
   // 0-index access
   void __setitem__ (object t, double val) {
     this->operator()(extract<int>(t[0]), extract<int>(t[1]) ) = val;}
   double __getitem__ (object t) {
     return this->operator()(extract<int>(t[0]), extract<int>(t[1]) );}
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(floatmatrix_overloads_resize, resize, 2, 3)

void (PyFloatMatrix::*solveForRhs_1)(const FloatArray &b, FloatArray &answer, bool) = &PyFloatMatrix::solveForRhs;
  void (PyFloatMatrix::*solveForRhs_2)(const FloatMatrix &b, FloatMatrix &answer, bool) = &PyFloatMatrix::solveForRhs;
void (PyFloatMatrix::*assemble_1)(const FloatMatrix &src, const IntArray &loc) = &PyFloatMatrix::assemble;
void (PyFloatMatrix::*assemble_2)(const FloatMatrix &src, const IntArray &row, const IntArray&col) = &PyFloatMatrix::assemble;

void pyclass_FloatMatrix()
{
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
}


/*****************************************************
* IntArray
*****************************************************/
class PyIntArray : public IntArray , public wrapper<IntArray>
{
public:
   PyIntArray () : IntArray () {}
   PyIntArray (int i) : IntArray(i) {}
   void __setitem__(int i, double val) { this->operator()(i) = val; }
   int __getitem__(int i) { return this->operator()(i); }
};
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(intarray_overloads_resize, resize, 1, 2)
void pyclass_IntArray()
{
  class_<PyIntArray, boost::noncopyable>("IntArray")
    .def(init< optional<int> >())
    .def("__len__", &IntArray::giveSize)
    .def("__getitem__", &PyIntArray::__getitem__) 
    .def("__setitem__", &PyIntArray::__setitem__)
    .def("resize", &IntArray::resize, intarray_overloads_resize())
    .def("isEmpty", &IntArray::isEmpty)
    .def("containsOnlyZeroes", &IntArray::containsOnlyZeroes)
    .def("contains", &IntArray::contains)
    .def("printYourself", &IntArray::printYourself)
    .def("zero", &IntArray::printYourself)
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

void pyclass_EngngModel()
{
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
    .def("giveCurrentStep", &EngngModel::giveCurrentStep, return_internal_reference<>())
    .def("giveContext", &EngngModel::giveContext, return_internal_reference<>())
    .def("giveNextStep",&EngngModel::giveNextStep, return_internal_reference<>())
    ;
}


/*****************************************************
* Domain
*****************************************************/
void pyclass_Domain()
{
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
  .def("giveMaterial", &Domain::giveMaterial, return_internal_reference<>())
  
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

  .def("giveNumberOfSpatialDimensions", &Domain::giveNumberOfSpatialDimensions)
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
  class_<PyElement, boost::noncopyable>("Element", init<int, Domain* >())
    .def("giveLabel", &Element::giveLabel)
    .def("giveNumber", &Element::giveNumber)
    .def("giveLocationArray", &PyElement::giveLocationArray)
    .def("invalidateLocationArray", &PyElement::invalidateLocationArray)
    .def("giveCharacteristicMatrix", &PyElement::giveCharacteristicMatrix)
    .def("giveCharacteristicVector", &PyElement::giveCharacteristicVector)
    .def("giveCharacteristicValue", &PyElement::giveCharacteristicValue)
    .def("computeNumberOfDofs", &PyElement::computeNumberOfDofs)
    .def("giveGeometryType", &PyElement::giveGeometryType)
    .def("giveDofManagerNumber", &Element::giveDofManagerNumber)
    .def("giveNumberOfIntegrationRules", &Element::giveNumberOfIntegrationRules)
    .def("giveDefaultIntegrationRulePtr", &Element::giveDefaultIntegrationRulePtr, return_internal_reference<>())
    ;
}


/*****************************************************
* DofManager
*****************************************************/
class PyDofManager : public DofManager, public wrapper<DofManager>
{
public:
  PyDofManager (int n, Domain *d) : DofManager (n,d) {}
  
  /*
  void giveUnknownVector(FloatArray &answer, const IntArray &dofMask,
			 EquationID type, ValueModeType mode, TimeStep *stepN) {
    if (override f = this->get_override("giveUnknownVector")) {
      f (answer, dofMask, type, mode, stepN); return;}
    DofManager::giveUnknownVector(answer, dofMask, type, mode, stepN);
  }
  void computeDofTransformation(FloatMatrix &answer, IntArray &dofIDArry, DofManTransfType mode) {
    if (override f = this->get_override("computeDofTransformation")) {f(answer, &dofIDArry, mode); return;}
    DofManager::computeDofTransformation(answer, &dofIDArry, mode);
  }

  void default_giveUnknownVector(FloatArray &answer, const IntArray &dofMask,
				 EquationID type, ValueModeType mode, TimeStep *stepN) 
  {return this->DofManager::giveUnknownVector(answer, dofMask, type, mode, stepN);}
  */
};

void (DofManager::*giveUnknownVector_1)(FloatArray &answer, const IntArray &dofMask,
					EquationID type, ValueModeType mode, TimeStep *stepN) = &DofManager::giveUnknownVector;

void pyclass_DofManager()
{
  class_<DofManager, boost::noncopyable>("DofManager", init<int, Domain* >())
    .def("hasCoordinates", &DofManager::hasCoordinates)
    .def("giveCoordinates", &DofManager::giveCoordinates, return_internal_reference<>())
    .def("giveCoordinate", &DofManager::giveCoordinate)
    .def("giveLabel", &DofManager::giveLabel)
    .def("giveUnknownVector", giveUnknownVector_1)
    .def("computeDofTransformation", &DofManager::computeDofTransformation)
    ;
}


/*****************************************************
* TimeStep
*****************************************************/
void pyclass_TimeStep()
{
  class_<TimeStep, boost::noncopyable>("TimeStep", no_init)
    .def("giveTargetTime", &TimeStep::giveTargetTime)
    .def("giveIntrinsicTime", &TimeStep::giveIntrinsicTime)
    .def("setTargetTime", &TimeStep::setTargetTime)
    .def("setIntrinsicTime", &TimeStep::setIntrinsicTime)
    .def("setTimeIncrement", &TimeStep::setTimeIncrement)
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
* GeneralBoundaryCondition
*****************************************************/
void pyclass_GeneralBoundaryCondition()
{
  class_<GeneralBoundaryCondition, boost::noncopyable >("GeneralBoundaryCondition", init<int, Domain*>())
    ;
}


/*****************************************************
* BoundaryCondition
*****************************************************/
void pyclass_BoundaryCondition()
{
  class_<BoundaryCondition, bases<GeneralBoundaryCondition>, boost::noncopyable >("BoundaryCondition", init<int, Domain*>())
    .def("setPrescribedValue", &BoundaryCondition::setPrescribedValue)
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
    ;
}


/*****************************************************
* Material
*****************************************************/
void pyclass_Material()
{
  class_<Material, boost::noncopyable >("Material", init<int, Domain*>())
    .def("giveIPValue", &Material::giveIPValue)
    .def("giveCharacteristicMatrix", &Material::giveCharacteristicMatrix)
    ;
}






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
    .value("EGT_line_2", EGT_line_2) /* line element with three nodes 1---2---3 */	\
    .value("EGT_triangle_1",EGT_triangle_1) /* triangle element with three nodes */	\
    .value("EGT_triangle_2", EGT_triangle_2) /* triangle element with 6 nodes */ \
    .value("EGT_quad_1",EGT_quad_1)   /* quadrialateral with 4 nodes */		\
    .value("EGT_quad_2", EGT_quad_2)   /* quadratic quadrialateral with 8 nodes */	\
    .value("EGT_tetra_1", EGT_tetra_1)  /* tetrahedron with 4 nodes */			\
    .value("EGT_hexa_1", EGT_hexa_1)   /* hexahedron with 8 nodes */			\
    .value("EGT_hexa_2", EGT_hexa_2)   /* hexahedron with 20 nodes */			\
    .value("EGT_Composite", EGT_Composite)/* Composite" geometry, vtk export supported by individual elements */ \
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




BOOST_PYTHON_MODULE (oofemlib)
{
  pyclass_FloatArray();
  pyclass_FloatMatrix();
  pyclass_IntArray();
  pyclass_SparseMtrx();
  pyclass_SpatialLocalizer();
  pyclass_EngngModelContext();
  pyclass_EngngModel();
  pyclass_Domain();
  pyclass_DataReader();
  pyclass_OOFEMTXTDataReader();
  pyclass_Element();
  pyclass_DofManager();
  pyclass_TimeStep();
  pyclass_DataStream();
  pyclass_Field();
  pyclass_FieldManager();
  pyclass_DofManValueField();
  pyclass_GeneralBoundaryCondition();
  pyclass_BoundaryCondition();
  pyclass_IntegrationRule();
  pyclass_GaussPoint();
  pyclass_Material();

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

  def("InstanciateProblem", InstanciateProblem_1, return_value_policy<manage_new_object>());
  //def("make_foo", make_foo, return_value_policy<manage_new_object>())
  def("createDofManValueFieldPtr", &DofManValueField_create, with_custodian_and_ward_postcall<0,2>());
  def("FieldManager_registerField", &FieldManager_registerField);

  implicitly_convertible<std::auto_ptr<DofManValueField>, std::auto_ptr<Field> >(); 
  register_ptr_to_python< std::auto_ptr<Field> >();
}
} // end namespace oofem
