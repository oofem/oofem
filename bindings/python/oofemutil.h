#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
namespace py = pybind11;
#include <string>
#include <iostream>
// for alternative tokens like 'or' 
#include <ciso646>

#include "oofemtxtdatareader.h"
#include "classfactory.h"
#include "logger.h"
#include "domain.h"
#include "outputmanager.h"
#include "modulemanager.h"
#include "set.h"


#ifdef _MSC_VER 
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

using namespace oofem;

oofem::OOFEMTXTInputRecord makeOOFEMTXTInputRecordFrom(py::kwargs &kw)
{
    py::dict tmp;
    std::string rec;
    std::string space= " ";
    
    for (auto item: kw) {
        std::string key = std::string(py::str(item.first));
        py::handle value = item.second;
        if ((strcasecmp(std::string(key).c_str(), "number") == 0) or (strcasecmp(std::string(key).c_str(), "domain") == 0))
            continue;
        if (key=="f_t") key="f(t)"; // handle f(t) loadTimeFunction field name
        rec.append(space);
        rec.append(key);
        if (PyBool_Check(value.ptr())) {
            // printf ("huhu\n");
            ;
        } else if (PyNumber_Check(value.ptr()) or PyUnicode_Check(value.ptr())) {
            PyObject* repr = PyObject_Str(value.ptr()); //PyObject_Repr puts single quotes around
            PyObject* str = PyUnicode_AsEncodedString(repr, "utf-8", "~E~");
            const char*bytes = PyBytes_AS_STRING(str) ;
            rec.append(space);
            rec.append(bytes);
            Py_XDECREF(repr);
            Py_XDECREF(str);
        } else if (PySequence_Check(value.ptr())) {
            PyObject *iterator = PyObject_GetIter(value.ptr());
            PyObject *item;
            if (iterator == NULL) {
                    ; /* propagate error */
            }
            unsigned int size = PySequence_Length(value.ptr());
            rec.append(space);
            rec.append(std::to_string(size));

            while ((item = PyIter_Next(iterator))) {
                /* do something with item */
                
                if (PyNumber_Check(item) || PyUnicode_Check(item)) {
                    PyObject* repr = PyObject_Repr(item);
                    PyObject* str = PyUnicode_AsEncodedString(repr, "utf-8", "~E~");
                    const char*bytes = PyBytes_AS_STRING(str) ;
                    rec.append(space);
                    rec.append(bytes);
                    Py_XDECREF(repr);
                    Py_XDECREF(str);
 
                } else {
                    // todo: error handling when conversion fails
                    py::handle h = item; // handle case when reference to object given; check reference counting
                    auto o = h.cast<py::object>();
                    oofem::FEMComponent* f = o.cast<oofem::FEMComponent*> ();
                    rec.append(space);
                    rec.append(std::to_string(f->giveNumber()));


                    //makeOOFEMTXTValueentryFrom (item);
                    
                }
               Py_DECREF(item);
            }
      
            Py_DECREF(iterator);
        } else { // handle case when reference to object given; check reference counting
            // todo: error handling when conversion fails
            py::handle h = value;
            auto o = h.cast<py::object>();
            oofem::FEMComponent* f = o.cast<oofem::FEMComponent*> ();
            rec.append(space);
            rec.append(std::to_string(f->giveNumber()));
        }
        /*
        if (PyErr_Occurred()) {
            // propagate error 
        } else {
            // continue doing useful work 
        }
        */  
   } 
    transform(rec.begin(), rec.end(), rec.begin(), ::tolower); // convert to lowercase, text probably in "" should not be converted (like filenames)
//     std::cout << rec << std::endl;
    oofem::OOFEMTXTInputRecord answer(0, rec) ;
    return answer;
}

oofem::OOFEMTXTInputRecord makeOutputManagerOOFEMTXTInputRecordFrom(py::kwargs kw)
{
    py::kwargs kw2;

    for (auto item: kw) {
        std::string key = std::string(py::str(item.first));
        py::handle value = item.second;

        transform(key.begin(), key.end(), key.begin(), ::tolower);
        
        if (key=="tstep_all" || key=="tstep_step" || key=="tsteps_out" || key=="dofman_all" || key=="dofman_output" || key=="dofman_except" || key=="element_all" || key=="element_output" || key=="element_except" || key=="pythonexport") {
            kw2[key.c_str()] = value;
        }
        kw2[""] = "outputmanager";
    }
    return makeOOFEMTXTInputRecordFrom(kw);

}


/*****************************************************
* EngngModel
*****************************************************/

py::object createEngngModelOfType(const char* type, py::args args, py::kwargs kw)
{
    //args
    int number = len(args)>0? PyLong_AsUnsignedLong(args[0].ptr()) : 0;
    oofem::EngngModel* master = len(args)>1? args[1].cast<oofem::EngngModel *>() : nullptr;
    std::unique_ptr<EngngModel> engngm = classFactory.createEngngModel(type,number,master);
    if ( !engngm ) { oofem::OOFEM_LOG_ERROR("engngModel: wrong input data"); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    // instanciateYourself
    ///@todo Output filename isn't stored like this (and has never been!)!?
    std::string outFile;
    if ( ir.hasField("outfile") ) {
        ir.giveField(outFile, "outfile");
    } else {
        outFile = "oofem.out.XXXXXX";
    }
    
    //engngm->Instanciate_init(outFile.c_str(), engngm->giveNumberOfDomains());
    engngm->letOutputBaseFileNameBe(outFile);
    engngm->initializeFrom(ir);

    if ( ir.hasField(_IFT_EngngModel_nmsteps) ) {
      oofem::OOFEM_LOG_ERROR("engngModel: simulation with metasteps is not (yet) supported in Python");
    } else {
      engngm->instanciateDefaultMetaStep(ir);
    }

    engngm->Instanciate_init();
    //
    py::object ret = py::cast(engngm.release());
    /* ????????????????????
    // sets the last created engngModel as default one for further script
    temp_global["defaultEngngModel"] = ret;
    ???????????????????? */
    return ret;
}

py::object linearStatic(py::args args, py::kwargs kw) { return createEngngModelOfType("linearstatic", args, kw); }

py::object staticStructural(py::args args, py::kwargs kw) { return createEngngModelOfType("staticstructural", args, kw); }

py::object transientTransport(py::args args, py::kwargs kw) { return createEngngModelOfType("transienttransport", args, kw); }



/*****************************************************
* Domain
*****************************************************/
py::object domain(py::args args, py::kwargs kw)
{
    // args
    int number =             len(args)>0? PyLong_AsUnsignedLong(args[0].ptr()) : 0;
    int serialNumber =       len(args)>1? PyLong_AsUnsignedLong(args[1].ptr()) : 0;
    EngngModel *engngModel = len(args)>2? args[2].cast<oofem::EngngModel *>() : nullptr;
    domainType dType =       len(args)>3? args[3].cast<domainType>(): oofem::domainType::_unknownMode;
    auto d = std::make_unique<Domain>(number,serialNumber,engngModel);
    d->setDomainType(dType);
    // output manager record
    oofem::OOFEMTXTInputRecord omir = makeOutputManagerOOFEMTXTInputRecordFrom(kw);
    d->giveOutputManager()->initializeFrom(omir);
    py::object ret = py::cast(d.release());
    /* ????????????????????
    // sets the last created domain as default one for furtherscript
    temp_global["defaultDomain"] = ret;
    ???????????????????? */
    return ret;
}

/*****************************************************
* Element
*****************************************************/
py::object createElementOfType(const char* type, py::args args, py::kwargs kw)
{
    
    // extracts values from args if they are specified
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? args[1].cast<oofem::Domain *>() : nullptr;
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
    auto elem = oofem::classFactory.createElement(type,number,domain);
    // if elem==nullptr, something was wrong
    if (!elem) { oofem::OOFEM_LOG_ERROR("element: wrong input data"); }
    // sets globalNumber == number befor initializeFrom
    elem->setGlobalNumber(number);
    // construct OOFEMTXTInputRecord from bp::dict **kw
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    // pass input record to elem
    elem->initializeFrom(ir);
    // convert element to PyObject (expected by raw_function, which enables fun(*args,**kw) syntax in python)
    return py::cast(elem.release());
}

// specific elements
py::object beam2d(py::args args, py::kwargs &kw) { return createElementOfType("beam2d",args,kw); }
py::object truss1d(py::args args, py::kwargs &kw) { return createElementOfType("truss1d",args,kw); }
py::object trPlaneStress2d(py::args args, py::kwargs &kw) { return createElementOfType("trplanestress2d",args,kw); }
py::object planeStress2d(py::args args, py::kwargs &kw) { return createElementOfType("planestress2d",args,kw); }
py::object qBrick1ht(py::args args, py::kwargs &kw) { return createElementOfType("qbrick1ht",args,kw); }
py::object quadAxiSym1ht(py::args args, py::kwargs &kw) { return createElementOfType("quadaxisym1ht",args,kw); }
py::object lspace(py::args args, py::kwargs &kw) { return createElementOfType("lspace",args,kw); }
py::object tr1ht(py::args args, py::kwargs &kw) { return createElementOfType("tr1ht",args,kw); }
py::object quad1ht(py::args args, py::kwargs &kw) { return createElementOfType("quad1ht",args,kw); }
py::object qquad1ht(py::args args, py::kwargs &kw) { return createElementOfType("qquad1ht",args,kw); }




/*****************************************************
* DofManager
*****************************************************/

py::object createDofManagerOfType(const char*type, py::args args, py::kwargs &kw)
{
    unsigned long number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? args[1].cast<oofem::Domain *>() : nullptr;
    auto dofMan = oofem::classFactory.createDofManager(type,number,domain);
    if (!dofMan) { oofem::OOFEM_LOG_ERROR("dofManager: wrong input data"); }
    dofMan->setGlobalNumber(number);
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    dofMan->initializeFrom(ir);
    return py::cast(dofMan.release());
}


py::object node(py::args args, py::kwargs kw) { return createDofManagerOfType("node",args,kw); }


/*****************************************************
* GeneralBoundaryCondition
*****************************************************/
// generalBoundaryCondition(aClass,domain=defaultDomain,**kw)
py::object createGeneralBoundaryConditionOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? args[1].cast<oofem::Domain *>() : nullptr;
    auto bc = oofem::classFactory.createBoundaryCondition(type,number,domain);
    if (!bc) { oofem::OOFEM_LOG_ERROR("generalBoundaryCondition: wrong input data"); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    bc->initializeFrom(ir);
    return py::cast(bc.release());
}

py::object boundaryCondition(py::args args, py::kwargs kw) { return createGeneralBoundaryConditionOfType("boundarycondition",args,kw); }
py::object constantEdgeLoad(py::args args, py::kwargs kw) { return createGeneralBoundaryConditionOfType("constantedgeload",args,kw); }
py::object nodalLoad(py::args args, py::kwargs kw) { return createGeneralBoundaryConditionOfType("nodalload",args,kw); }
py::object structTemperatureLoad(py::args args, py::kwargs kw) { return createGeneralBoundaryConditionOfType("structtemperatureload",args,kw); }
py::object structEigenstrainLoad(py::args args, py::kwargs kw) { return createGeneralBoundaryConditionOfType("structEigenstrainLoad",args,kw); }
py::object constantSurfaceLoad(py::args args, py::kwargs kw) { return createGeneralBoundaryConditionOfType("constantsurfaceload",args,kw); }
py::object deadWeight(py::args args, py::kwargs kw) { return createGeneralBoundaryConditionOfType("DeadWeight",args,kw); }


/*****************************************************
* InitialCondition
*****************************************************/
py::object createInitialConditionOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? args[1].cast<oofem::Domain *>() : nullptr;
    auto ic = oofem::classFactory.createInitialCondition(type,number,domain);
    if (!ic) { oofem::OOFEM_LOG_ERROR("initialCondition: wrong input data"); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    ic->initializeFrom(ir);
    return py::cast(ic.release());
}

py::object initialCondition(py::args args, py::kwargs kw) { return createInitialConditionOfType("initialcondition",args,kw); }


/*****************************************************
* Material
*****************************************************/
py::object createMaterialOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? args[1].cast<oofem::Domain *>() : nullptr;
    auto mat = oofem::classFactory.createMaterial(type,number,domain);
    if (!mat) { oofem::OOFEM_LOG_ERROR("material: wrong input data"); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    mat->initializeFrom(ir);
    return py::cast(mat.release());
}

py::object isoLE(py::args args, py::kwargs kw) { return createMaterialOfType("isole",args,kw); }
py::object idm1(py::args args, py::kwargs kw) { return createMaterialOfType("idm1",args,kw); }
py::object isoHeat(py::args args, py::kwargs kw) { return createMaterialOfType("isoheat",args,kw); }
py::object j2mat(py::args args, py::kwargs kw) { return createMaterialOfType("j2mat",args,kw); }
py::object steel1(py::args args, py::kwargs kw) { return createMaterialOfType("steel1",args,kw); }


/*****************************************************
* CrossSection
*****************************************************/
py::object createCrossSectionOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? args[1].cast<oofem::Domain *>() : nullptr;
    auto cs = oofem::classFactory.createCrossSection(type,number,domain);
    if (!cs) { oofem::OOFEM_LOG_ERROR("crossSection: wrong input data"); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    cs->initializeFrom(ir);
    return py::cast(cs.release());
}

py::object simpleCS(py::args args, py::kwargs kw) { return createCrossSectionOfType("simplecs",args,kw); }
py::object simpleTransportCS(py::args args, py::kwargs kw) { return createCrossSectionOfType("simpletransportcs",args,kw); }



/*****************************************************
* LoadTimeFunction
*****************************************************/
py::object createLoadTimeFunctionOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? args[1].cast<oofem::Domain *>() : nullptr;
    auto ltf = oofem::classFactory.createFunction(type,number,domain);
    if (!ltf) { oofem::OOFEM_LOG_ERROR("function: wrong input data"); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    ltf->initializeFrom(ir);
    return py::cast(ltf.release());
}

py::object peakFunction(py::args args, py::kwargs kw) { return createLoadTimeFunctionOfType("peakfunction",args,kw); }

py::object constantFunction(py::args args, py::kwargs kw) { return createLoadTimeFunctionOfType("constantfunction",args,kw); }

py::object piecewiseLinFunction(py::args args, py::kwargs kw) { return createLoadTimeFunctionOfType("piecewiselinfunction",args,kw); }


/*****************************************************
* ExportModule
*****************************************************/
py::object createExportModuleOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::EngngModel *engngModel = len(args)>1? args[1].cast<oofem::EngngModel *>() : nullptr;
    auto module = oofem::classFactory.createExportModule(type,number,engngModel);
    if (!module) { oofem::OOFEM_LOG_ERROR("exportModule: wrong input data"); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    module->initializeFrom(ir);
    module->initialize();
    engngModel->giveExportModuleManager()->registerModule(module);
    return py::cast(engngModel->giveExportModuleManager()->giveModule(engngModel->giveExportModuleManager()->giveNumberOfModules()));
}

py::object vtkxml(py::args args, py::kwargs kw) { return createExportModuleOfType("vtkxml",args,kw); }
py::object homExport(py::args args, py::kwargs kw) { return createExportModuleOfType("hom",args,kw); }
py::object vtkmemory(py::args args, py::kwargs kw) { return createExportModuleOfType("vtkmemory",args,kw); }

/*****************************************************
* Sets
*****************************************************/
py::object createSetOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? args[1].cast<oofem::Domain *>() : nullptr;
    std::unique_ptr<Set> setP = std::make_unique<Set>(number, domain);
    if ( !setP ) { oofem::OOFEM_LOG_ERROR("Couldn't create set: %d", number); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    setP->initializeFrom(ir);
    return py::cast(setP.release());
}

py::object createSet(py::args args, py::kwargs kw) { return createSetOfType("set",args,kw); }
