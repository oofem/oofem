#ifndef pythonmaterial_h
#define pythonmaterial_h

#include "material.h"
#include "matstatus.h"
#include <string>
#include <memory>

#ifdef _USE_NANOBIND
    #include <nanobind/nanobind.h>
    namespace nb = nanobind;
#elif defined(_PYBIND_BINDINGS)
    #include <pybind11/pybind11.h>
    namespace py = pybind11;
#endif

#define _IFT_PythonMaterial_Name "pythonmaterial"
#define _IFT_PythonMaterial_module "module"
#define _IFT_PythonMaterial_object "object"

namespace oofem {

class OOFEM_EXPORT PythonMaterialStatus : public MaterialStatus
{
protected:
#ifdef _USE_NANOBIND
    nb::dict stateDict, tempStateDict;
#elif defined(_PYBIND_BINDINGS)
    py::dict stateDict, tempStateDict;
#endif

public:
    PythonMaterialStatus(GaussPoint * gp);
    ~PythonMaterialStatus() override;

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

#ifdef _USE_NANOBIND
    nb::dict giveStateDictionary() const { return stateDict; }
    nb::dict giveTempStateDictionary() const { return tempStateDict; }
#elif defined(_PYBIND_BINDINGS)
    py::dict giveStateDictionary() const { return stateDict; }
    py::dict giveTempStateDictionary() const { return tempStateDict; }
#endif

    const char *giveClassName() const override { return "PythonMaterialStatus"; }
};

class OOFEM_EXPORT PythonMaterial : public Material
{
protected:
    std::string moduleName;
    std::string objectName;

#ifdef _USE_NANOBIND
    nb::object pyObject;
    nb::object pyHasMaterialModeCapability;
    nb::object pyGiveCharacteristicMatrix;
    nb::object pyGiveCharacteristicVector;
    nb::object pyGiveCharacteristicValue;
    nb::object pyPrintOutputAt;
    nb::object pyUpdateTempState;
    nb::object pyGiveStateVariableIDs;
#elif defined(_PYBIND_BINDINGS)
    py::object pyObject;
    py::object pyHasMaterialModeCapability;
    py::object pyGiveCharacteristicMatrix;
    py::object pyGiveCharacteristicVector;
    py::object pyGiveCharacteristicValue;
    py::object pyPrintOutputAt;
    py::object pyUpdateTempState;
    py::object pyGiveStateVariableIDs;
#endif

public:
    PythonMaterial(int n, Domain *d);
    ~PythonMaterial() override;

    const char *giveClassName() const override { return "PythonMaterial"; }
    const char *giveInputRecordName() const override { return "PythonMaterial"; }

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;
    void postInitialize() override;
    
    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    
    void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override;
    void giveCharacteristicVector(FloatArray &answer, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override;
    double giveCharacteristicValue(MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override;

    /**
     * Forwards the state push to the python object's updateTempState.
     *
     * The material decomposes the state itself into the entries of the temporary state dictionary
     * it uses -- that dictionary is a python material's state store. postInitialize requires the
     * method to exist, since the state is no longer passed to giveCharacteristicVector.
     */
    void updateTempState(const FloatArray &stateVector, GaussPoint *gp, TimeStep *tStep) override;

    /**
     * Forwards to the python object's giveStateVariableIDs, if it defines one; returns an empty
     * layout otherwise, meaning the material does not advertise one. The python side returns a
     * sequence of (FieldType, StateOperator) pairs.
     */
    StateVariableLayout giveStateVariableIDs(MaterialMode mmode) const override;

    void printOutputAt(FILE *file, TimeStep *tStep, const PythonMaterialStatus *status) const;
};

} // end namespace oofem
#endif // pythonmaterial_h
