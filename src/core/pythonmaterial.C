#include "pythonmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "timestep.h"
#include "error.h"
#include "classfactory.h"

#ifdef _USE_NANOBIND
    #include <nanobind/nanobind.h>
    #include "../../bindings/python/oofemarray-nanobind.h"
    namespace nb = nanobind;
#elif defined(_PYBIND_BINDINGS)
    #include <pybind11/embed.h>
    #include <pybind11/numpy.h>
    namespace py = pybind11;
#endif

namespace oofem {

REGISTER_Material(PythonMaterial);

#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
namespace {
class PythonInterpreterManager {
public:
    static void initialize() {
        static PythonInterpreterManager instance;
    }
private:
    PythonInterpreterManager() {
        if (!Py_IsInitialized()) {
            Py_Initialize();
            PyEval_SaveThread(); // Release GIL
        }
    }
    ~PythonInterpreterManager() = default;
};
} // end anonymous namespace
#endif

PythonMaterialStatus::PythonMaterialStatus(GaussPoint * gp) : MaterialStatus(gp)
{
}

PythonMaterialStatus::~PythonMaterialStatus()
{
#if defined(_USE_NANOBIND)
    nb::gil_scoped_acquire gil;
    stateDict.reset();
    tempStateDict.reset();
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    stateDict.release().dec_ref();
    tempStateDict.release().dec_ref();
#endif
}

void PythonMaterialStatus::initTempStatus()
{
    MaterialStatus::initTempStatus();
#ifdef _USE_NANOBIND
    nb::gil_scoped_acquire gil;
    if (stateDict.ptr() != nullptr) {
        tempStateDict = nb::cast<nb::dict>(stateDict.attr("copy")());
    } else {
        tempStateDict = nb::dict();
    }
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    if (stateDict.ptr() != nullptr) {
        tempStateDict = py::cast<py::dict>(stateDict.attr("copy")());
    } else {
        tempStateDict = py::dict();
    }
#endif
}

void PythonMaterialStatus::updateYourself(TimeStep *tStep)
{
    MaterialStatus::updateYourself(tStep);
#ifdef _USE_NANOBIND
    nb::gil_scoped_acquire gil;
    if (tempStateDict.ptr() != nullptr) {
        stateDict = nb::cast<nb::dict>(tempStateDict.attr("copy")());
    } else {
        stateDict = nb::dict();
    }
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    if (tempStateDict.ptr() != nullptr) {
        stateDict = py::cast<py::dict>(tempStateDict.attr("copy")());
    } else {
        stateDict = py::dict();
    }
#endif
}

PythonMaterial::PythonMaterial(int n, Domain *d) : Material(n, d)
{
}

PythonMaterial::~PythonMaterial()
{
#if defined(_USE_NANOBIND)
    nb::gil_scoped_acquire gil;
    pyObject.reset();
    pyHasMaterialModeCapability.reset();
    pyGiveCharacteristicMatrix.reset();
    pyGiveCharacteristicVector.reset();
    pyGiveCharacteristicValue.reset();
    pyPrintOutputAt.reset();
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    pyObject.release().dec_ref();
    pyHasMaterialModeCapability.release().dec_ref();
    pyGiveCharacteristicMatrix.release().dec_ref();
    pyGiveCharacteristicVector.release().dec_ref();
    pyGiveCharacteristicValue.release().dec_ref();
    pyPrintOutputAt.release().dec_ref();
#endif
}

void PythonMaterial::initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    Material::initializeFrom(ir);
    IR_GIVE_FIELD(ir, moduleName, _IFT_PythonMaterial_module);
    IR_GIVE_FIELD(ir, objectName, _IFT_PythonMaterial_object);
}

void PythonMaterial::postInitialize()
{
    Material::postInitialize();

#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    PythonInterpreterManager::initialize();
#ifdef _USE_NANOBIND
    nb::gil_scoped_acquire gil;
    try {
        nb::module_ calc = nb::module_::import_(moduleName.c_str());
        if (!nb::hasattr(calc, objectName.c_str())) {
            OOFEM_ERROR("PythonMaterial: module '%s' does not have object '%s'.", moduleName.c_str(), objectName.c_str());
        }
        pyObject = calc.attr(objectName.c_str());
        
        if (!nb::hasattr(pyObject, "hasMaterialModeCapability") ||
            !nb::hasattr(pyObject, "giveCharacteristicMatrix") ||
            !nb::hasattr(pyObject, "giveCharacteristicVector") ||
            !nb::hasattr(pyObject, "giveCharacteristicValue")) {
            OOFEM_ERROR("PythonMaterial: object '%s' is missing required methods.", objectName.c_str());
        }

        pyHasMaterialModeCapability = pyObject.attr("hasMaterialModeCapability");
        pyGiveCharacteristicMatrix = pyObject.attr("giveCharacteristicMatrix");
        pyGiveCharacteristicVector = pyObject.attr("giveCharacteristicVector");
        pyGiveCharacteristicValue = pyObject.attr("giveCharacteristicValue");
        if (nb::hasattr(pyObject, "printOutputAt")) {
            pyPrintOutputAt = pyObject.attr("printOutputAt");
        }
        // Push half of the interface. Required: the state is no longer handed to
        // giveCharacteristicVector, so a material that still expects it there would be called
        // with the wrong arguments. Report it here, with what to do about it, rather than let it
        // surface as a python TypeError deep in the assembly loop.
        if (!nb::hasattr(pyObject, "updateTempState")) {
            OOFEM_ERROR("PythonMaterial: object '%s' has no updateTempState. The material state is "
                        "now pushed once per iteration instead of being deposited inside "
                        "giveCharacteristicVector: move the state deposit into "
                        "updateTempState(self, stateVector, gp, tStep, stateDict, tempStateDict) "
                        "and drop the leading flux argument of "
                        "giveCharacteristicVector(self, mode, gp, tStep, stateDict, tempStateDict).",
                        objectName.c_str());
        }
        pyUpdateTempState = pyObject.attr("updateTempState");
        if (nb::hasattr(pyObject, "giveStateVariableIDs")) {
            pyGiveStateVariableIDs = pyObject.attr("giveStateVariableIDs");
        }
    } catch (const std::exception &e) {
        OOFEM_ERROR("PythonMaterial: initialization failed: %s", e.what());
    }
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    try {
        py::module calc = py::module::import(moduleName.c_str());
        if (!py::hasattr(calc, objectName.c_str())) {
            OOFEM_ERROR("PythonMaterial: module '%s' does not have object '%s'.", moduleName.c_str(), objectName.c_str());
        }
        pyObject = calc.attr(objectName.c_str());
        
        if (!py::hasattr(pyObject, "hasMaterialModeCapability") ||
            !py::hasattr(pyObject, "giveCharacteristicMatrix") ||
            !py::hasattr(pyObject, "giveCharacteristicVector") ||
            !py::hasattr(pyObject, "giveCharacteristicValue")) {
            OOFEM_ERROR("PythonMaterial: object '%s' is missing required methods.", objectName.c_str());
        }

        pyHasMaterialModeCapability = pyObject.attr("hasMaterialModeCapability");
        pyGiveCharacteristicMatrix = pyObject.attr("giveCharacteristicMatrix");
        pyGiveCharacteristicVector = pyObject.attr("giveCharacteristicVector");
        pyGiveCharacteristicValue = pyObject.attr("giveCharacteristicValue");
        if (py::hasattr(pyObject, "printOutputAt")) {
            pyPrintOutputAt = pyObject.attr("printOutputAt");
        }
        // Push half of the interface. Required: the state is no longer handed to
        // giveCharacteristicVector, so a material that still expects it there would be called
        // with the wrong arguments. Report it here, with what to do about it, rather than let it
        // surface as a python TypeError deep in the assembly loop.
        if (!py::hasattr(pyObject, "updateTempState")) {
            OOFEM_ERROR("PythonMaterial: object '%s' has no updateTempState. The material state is "
                        "now pushed once per iteration instead of being deposited inside "
                        "giveCharacteristicVector: move the state deposit into "
                        "updateTempState(self, stateVector, gp, tStep, stateDict, tempStateDict) "
                        "and drop the leading flux argument of "
                        "giveCharacteristicVector(self, mode, gp, tStep, stateDict, tempStateDict).",
                        objectName.c_str());
        }
        pyUpdateTempState = pyObject.attr("updateTempState");
        if (py::hasattr(pyObject, "giveStateVariableIDs")) {
            pyGiveStateVariableIDs = pyObject.attr("giveStateVariableIDs");
        }
    } catch (const std::exception &e) {
        OOFEM_ERROR("PythonMaterial: initialization failed: %s", e.what());
    }
#endif

#endif
}

std::unique_ptr<MaterialStatus> PythonMaterial::CreateStatus(GaussPoint *gp) const
{
#if defined(_USE_NANOBIND)
    nb::gil_scoped_acquire gil;
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
#endif
    return std::make_unique<PythonMaterialStatus>(gp);
}

int PythonMaterial::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    auto ms = static_cast<PythonMaterialStatus *>(this->giveStatus(gp));
    std::string key = std::to_string(type);

#ifdef _USE_NANOBIND
    nb::gil_scoped_acquire gil;
    if (ms->giveStateDictionary().contains(nb::str(key.c_str()))) {
        nb::object val = ms->giveStateDictionary()[nb::str(key.c_str())];
        try {
            answer.resize(1);
            answer.at(1) = nb::cast<double>(val);
            return 1;
        } catch (...) {}
        try {
            answer = nb::cast<FloatArray>(val);
            return 1;
        } catch (...) {}
        OOFEM_WARNING("Dictionary entry of material state not double or FloatArray");
    }
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    if (ms->giveStateDictionary().contains(key)) {
        py::object val = ms->giveStateDictionary()[key.c_str()];
        try {
            answer.resize(1);
            answer.at(1) = val.cast<double>();
            return 1;
        } catch (...) {}
        try {
            answer = val.cast<FloatArray>();
            return 1;
        } catch (...) {}
        OOFEM_WARNING("Dictionary entry of material state not double or FloatArray");
    }
#endif

#endif

    return Material::giveIPValue(answer, gp, type, tStep);
}

bool PythonMaterial::hasMaterialModeCapability(MaterialMode mode) const
{
#ifdef _USE_NANOBIND
    nb::gil_scoped_acquire gil;
    nb::object result = pyHasMaterialModeCapability(nb::cast(mode));
    return nb::cast<bool>(result);
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    py::object result = pyHasMaterialModeCapability(mode);
    return result.cast<bool>();
#else
    OOFEM_ERROR("Not compiled with python support.");
    return false;
#endif
}

void PythonMaterial::giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    auto ms = static_cast<PythonMaterialStatus *>(this->giveStatus(gp));
#endif

#ifdef _USE_NANOBIND
    nb::gil_scoped_acquire gil;
    nb::object result = pyGiveCharacteristicMatrix(nb::cast(type), nb::cast(gp), nb::cast(tStep), ms->giveStateDictionary(), ms->giveTempStateDictionary());
    answer = nb::cast<FloatMatrix>(result);
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    py::object result = pyGiveCharacteristicMatrix(type, gp, tStep, ms->giveStateDictionary(), ms->giveTempStateDictionary());
    answer = result.cast<FloatMatrix>();
#else
    OOFEM_ERROR("Not compiled with python support.");
#endif
}

void PythonMaterial::giveCharacteristicVector(FloatArray &answer, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    auto ms = static_cast<PythonMaterialStatus *>(this->giveStatus(gp));
#endif

#ifdef _USE_NANOBIND
    nb::gil_scoped_acquire gil;
    nb::object result = pyGiveCharacteristicVector(nb::cast(type), nb::cast(gp), nb::cast(tStep), ms->giveStateDictionary(), ms->giveTempStateDictionary());
    answer = nb::cast<FloatArray>(result);
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    py::object result = pyGiveCharacteristicVector(type, gp, tStep, ms->giveStateDictionary(), ms->giveTempStateDictionary());
    answer = result.cast<FloatArray>();
#else
    OOFEM_ERROR("Not compiled with python support.");
#endif
}

void PythonMaterial::updateTempState(const FloatArray &stateVector, GaussPoint *gp, TimeStep *tStep)
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    auto ms = static_cast<PythonMaterialStatus *>(this->giveStatus(gp));
#endif

    // The material decomposes the state itself, into whatever entries of the temporary state
    // dictionary it uses -- that dictionary is a python material's state store, so there is
    // nothing for this side to unpack or cache.
#ifdef _USE_NANOBIND
    nb::gil_scoped_acquire gil;
    pyUpdateTempState(nb::cast(stateVector), nb::cast(gp), nb::cast(tStep), ms->giveStateDictionary(), ms->giveTempStateDictionary());
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    pyUpdateTempState(stateVector, gp, tStep, ms->giveStateDictionary(), ms->giveTempStateDictionary());
#else
    OOFEM_ERROR("Not compiled with python support.");
#endif
}

StateVariableLayout PythonMaterial::giveStateVariableIDs(MaterialMode mmode) const
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    if ( pyGiveStateVariableIDs ) {
#ifdef _USE_NANOBIND
        nb::gil_scoped_acquire gil;
        nb::object result = pyGiveStateVariableIDs(nb::cast(mmode));
#else
        py::gil_scoped_acquire gil;
        py::object result = pyGiveStateVariableIDs(mmode);
#endif
        // The python side returns a sequence of (field, operator) pairs. Converted element by
        // element: handing the sequence straight to a caster would ask pybind to materialize a
        // temporary, which it refuses outside a bound function.
        StateVariableLayout layout;
        for ( auto item : result ) {
            auto pair = item;
            layout.push_back( { (FieldType) pair[ py::int_(0) ].template cast<int>(),
                                (StateOperator) pair[ py::int_(1) ].template cast<int>() } );
        }
        return layout;
    }
#endif
    return StateVariableLayout();
}

double PythonMaterial::giveCharacteristicValue(MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    auto ms = static_cast<PythonMaterialStatus *>(this->giveStatus(gp));
#endif

#ifdef _USE_NANOBIND
    nb::gil_scoped_acquire gil;
    nb::object result = pyGiveCharacteristicValue(nb::cast(type), nb::cast(gp), nb::cast(tStep), ms->giveStateDictionary(), ms->giveTempStateDictionary());
    return nb::cast<double>(result);
#elif defined(_PYBIND_BINDINGS)
    py::gil_scoped_acquire gil;
    py::object result = pyGiveCharacteristicValue(type, gp, tStep, ms->giveStateDictionary(), ms->giveTempStateDictionary());
    return result.cast<double>();
#else
    OOFEM_ERROR("Not compiled with python support.");
    return 0.0;
#endif
}

void PythonMaterial::printOutputAt(FILE *file, TimeStep *tStep, const PythonMaterialStatus *status) const
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    if (pyPrintOutputAt) {
        // Flush the C++ stream buffer to ensure output order is maintained
        // before the Python part writes to the same file descriptor.
        fflush(file);

#ifdef _USE_NANOBIND
        nb::gil_scoped_acquire gil;
        int fd = fileno(file);
        try {
            pyPrintOutputAt(fd, nb::cast(tStep), status->giveStateDictionary());
        } catch (const std::exception &e) {
            OOFEM_ERROR("PythonMaterial::printOutputAt failed: %s", e.what());
        }
#elif defined(_PYBIND_BINDINGS)
        py::gil_scoped_acquire gil;
        int fd = fileno(file);
        try {
            pyPrintOutputAt(fd, tStep, status->giveStateDictionary());
        } catch (const std::exception &e) {
            OOFEM_ERROR("PythonMaterial::printOutputAt failed: %s", e.what());
        }
#endif
    }
#endif
}

void PythonMaterialStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    const Material *mat_base = gp->giveCrossSection()->giveMaterial(gp);
    const PythonMaterial *mat = static_cast<const PythonMaterial *>(mat_base);
    mat->printOutputAt(file, tStep, this);
}

} // end namespace oofem
