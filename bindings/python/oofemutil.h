#ifdef _USE_NANOBIND
    #include <nanobind/nanobind.h>
    #include <nanobind/operators.h>
    #include <nanobind/stl/string.h>
    namespace py = nanobind;
    #define PY_STD_STRING(s) std::string(s.c_str())
    #define PY_CAST(a,b) py::cast<a>(b)
#else
    #include <pybind11/pybind11.h>
    #include <pybind11/operators.h>
    namespace py = pybind11;
    #define PY_STD_STRING(s) std::string(s)
    #define PY_CAST(a,b) b.cast<a>()
#endif
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


#include "feinterpol.h"
#include "fei2dquadlin.h"
#include "skyline.h"
#include "ldltfact.h"

#ifdef __MPM_MODULE
#include "prototype2.h"
#endif

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
        std::string key = PY_STD_STRING(py::str(item.first));
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
                    #if _USE_NANOBIND
                        oofem::FEMComponent* f = py::cast<oofem::FEMComponent*>(h);
                    #else
                        auto o = h.cast<py::object>();
                        oofem::FEMComponent* f = o.cast<oofem::FEMComponent*> ();
                    #endif
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
            #if _USE_NANOBIND
                oofem::FEMComponent* f = py::cast<oofem::FEMComponent*>(h);
            #else
                auto o = h.cast<py::object>();
                oofem::FEMComponent* f = o.cast<oofem::FEMComponent*> ();
            #endif
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
        std::string key = PY_STD_STRING(py::str(item.first));
        py::handle value = item.second;

        transform(key.begin(), key.end(), key.begin(), ::tolower);
        
        if (key=="tstep_all" || key=="tstep_step" || key=="tsteps_out" || key=="dofman_all" || key=="dofman_output" || key=="dofman_except" || key=="element_all" || key=="element_output" || key=="element_except" || key=="pythonexport") {
            kw2[key.c_str()] = value;
        }
        kw2[""] = "outputmanager";
    }
    return makeOOFEMTXTInputRecordFrom(kw);

}


#define OOFEM_RAISE(msg) { oofem::OOFEM_LOG_ERROR("%s",msg); throw std::runtime_error(msg); }


/*****************************************************
* EngngModel
*****************************************************/

py::object createEngngModelOfType(const char* type, py::args args, py::kwargs kw)
{
    //args
    int number = len(args)>0? PyLong_AsUnsignedLong(args[0].ptr()) : 0;
    oofem::EngngModel* master = len(args)>1? PY_CAST(oofem::EngngModel *,args[1]) : nullptr;
    std::unique_ptr<EngngModel> engngm = classFactory.createEngngModel(type,number,master);
    if ( !engngm ) { OOFEM_RAISE("engngModel: wrong input data"); }
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
        OOFEM_RAISE("engngModel: simulation with metasteps is not (yet) supported in Python");
    } else {
      //engngm->instanciateDefaultMetaStep(ir);
        engngm->giveTimeStepController()->instanciateDefaultMetaStep(ir);
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

py::object dummyProblem(py::args args, py::kwargs kw) { return createEngngModelOfType("dummy", args, kw); }
py::object mpmProblem(py::args args, py::kwargs kw) { return createEngngModelOfType("mpmproblem", args, kw); }


/*****************************************************
* Domain
*****************************************************/
py::object domain(py::args args, py::kwargs kw)
{
    // args
    int number =             len(args)>0? PyLong_AsUnsignedLong(args[0].ptr()) : 0;
    int serialNumber =       len(args)>1? PyLong_AsUnsignedLong(args[1].ptr()) : 0;
    EngngModel *engngModel = len(args)>2? PY_CAST(oofem::EngngModel *,args[2]) : nullptr;
    domainType dType =       len(args)>3? PY_CAST(domainType,args[3]): oofem::domainType::_unknownMode;
    auto d = std::make_unique<Domain>(number,serialNumber,engngModel);
    d->setDomainType(dType);
    // output manager record
    oofem::OOFEMTXTInputRecord omir = makeOutputManagerOOFEMTXTInputRecordFrom(kw);
    if ( !engngModel->giveSuppressOutput() ) {
        d->giveOutputManager()->initializeFrom(omir);
    }
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
    oofem::Domain* domain = len(args)>1? PY_CAST(oofem::Domain*,args[1]) : nullptr;
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
    if (!elem) { OOFEM_RAISE("element: wrong input data"); }
    // sets globalNumber == number befor initializeFrom
    elem->setGlobalNumber(number);
    // construct OOFEMTXTInputRecord from bp::dict **kw
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    // pass input record to elem
    elem->initializeFrom(ir, 1);
    // convert element to PyObject (expected by raw_function, which enables fun(*args,**kw) syntax in python)
    return py::cast(elem.release());
}

// specific elements
py::object truss1d(py::args args, py::kwargs &kw) { return createElementOfType("truss1d",args,kw); }
py::object truss2d(py::args args, py::kwargs &kw) { return createElementOfType("truss2d",args,kw); }
py::object truss3d(py::args args, py::kwargs &kw) { return createElementOfType("truss3d",args,kw); }
py::object beam2d(py::args args, py::kwargs &kw) { return createElementOfType("beam2d",args,kw); }
py::object beam3d(py::args args, py::kwargs &kw) { return createElementOfType("beam3d",args,kw); }
py::object lattice2d(py::args args, py::kwargs &kw) { return createElementOfType("lattice2d",args,kw); }
py::object lattice2dboundary(py::args args, py::kwargs &kw) { return createElementOfType("lattice2dboundary",args,kw); }
py::object lattice3d(py::args args, py::kwargs &kw) { return createElementOfType("lattice3d",args,kw); }
py::object lattice3dboundary(py::args args, py::kwargs &kw) { return createElementOfType("lattice3dboundary",args,kw); }
py::object latticelink3d(py::args args, py::kwargs &kw) { return createElementOfType("latticelink3d",args,kw); }
py::object latticelink3dboundary(py::args args, py::kwargs &kw) { return createElementOfType("latticelink3dboundary",args,kw); }
py::object planeStress2d(py::args args, py::kwargs &kw) { return createElementOfType("planestress2d",args,kw); }
py::object linquad3dplanestress(py::args args, py::kwargs &kw) { return createElementOfType("linquad3dplanestress",args,kw); }
py::object qPlaneStress2d(py::args args, py::kwargs &kw) { return createElementOfType("qplanestress2d",args,kw); }
py::object trPlaneStress2d(py::args args, py::kwargs &kw) { return createElementOfType("trplanestress2d",args,kw); }
py::object qTrPlStr(py::args args, py::kwargs &kw) { return createElementOfType("qtrplstr",args,kw); }
py::object trPlaneStrRot(py::args args, py::kwargs &kw) { return createElementOfType("trplanestrrot",args,kw); }
py::object trPlaneStressRotAllman(py::args args, py::kwargs &kw) { return createElementOfType("trplanestressrotallman",args,kw); }
py::object trPlanestressrotallman3d(py::args args, py::kwargs &kw) { return createElementOfType("trplanestressrotallman3d",args,kw); }
py::object quad1PlaneStrain(py::args args, py::kwargs &kw) { return createElementOfType("quad1planestrain",args,kw); }
py::object trPlaneStrain(py::args args, py::kwargs &kw) { return createElementOfType("trplanestrain",args,kw); }
py::object axisymm3d(py::args args, py::kwargs &kw) { return createElementOfType("axisymm3d",args,kw); }
py::object q4axisymm(py::args args, py::kwargs &kw) { return createElementOfType("q4axisymm",args,kw); } 
py::object l4axisymm(py::args args, py::kwargs &kw) { return createElementOfType("l4axisymm",args,kw); }
py::object lspace(py::args args, py::kwargs &kw) { return createElementOfType("lspace",args,kw); }
py::object lspaceBB(py::args args, py::kwargs &kw) { return createElementOfType("lspacebb",args,kw); }
py::object qspace(py::args args, py::kwargs &kw) { return createElementOfType("qspace",args,kw); }
py::object ltrspace(py::args args, py::kwargs &kw) { return createElementOfType("ltrspace",args,kw); }
py::object qTRSpace(py::args args, py::kwargs &kw) { return createElementOfType("qtrspace",args,kw); } 
py::object lWedge(py::args args, py::kwargs &kw) { return createElementOfType("lwedge",args,kw); }
py::object qWedge(py::args args, py::kwargs &kw) { return createElementOfType("qwedge",args,kw); } 
//transport elements
py::object line1ht(py::args args, py::kwargs &kw) { return createElementOfType("line1ht",args,kw); }
py::object line1mt(py::args args, py::kwargs &kw) { return createElementOfType("line1mt",args,kw); }
py::object line1hmt(py::args args, py::kwargs &kw) { return createElementOfType("line1hmt",args,kw); }
py::object tr1ht(py::args args, py::kwargs &kw) { return createElementOfType("tr1ht",args,kw); }
py::object tr1mt(py::args args, py::kwargs &kw) { return createElementOfType("tr1mt",args,kw); }
py::object tr1hmt(py::args args, py::kwargs &kw) { return createElementOfType("tr1hmt",args,kw); }
py::object quad1ht(py::args args, py::kwargs &kw) { return createElementOfType("quad1ht",args,kw); }
py::object quad1mt(py::args args, py::kwargs &kw) { return createElementOfType("quad1mt",args,kw); }
py::object quad1hmt(py::args args, py::kwargs &kw) { return createElementOfType("quad1hmt",args,kw); }
py::object qquad1ht(py::args args, py::kwargs &kw) { return createElementOfType("qquad1ht",args,kw); }
py::object qquad1mt(py::args args, py::kwargs &kw) { return createElementOfType("qquad1mt",args,kw); }
py::object qquad1hmt(py::args args, py::kwargs &kw) { return createElementOfType("qquad1hmt",args,kw); }
py::object quadAxiSym1ht(py::args args, py::kwargs &kw) { return createElementOfType("quadaxisym1ht",args,kw); }
py::object traxisym1ht(py::args args, py::kwargs &kw) { return createElementOfType("traxisym1ht",args,kw); }
py::object tetrah1ht(py::args args, py::kwargs &kw) { return createElementOfType("tetrah1ht",args,kw); }
py::object brick1ht(py::args args, py::kwargs &kw) { return createElementOfType("brick1ht",args,kw); }
py::object brick1hmt(py::args args, py::kwargs &kw) { return createElementOfType("brick1hmt",args,kw); }
py::object qBrick1ht(py::args args, py::kwargs &kw) { return createElementOfType("qbrick1ht",args,kw); }
py::object qBrick1mt(py::args args, py::kwargs &kw) { return createElementOfType("qbrick1mt",args,kw); }
py::object qBrick1hmt(py::args args, py::kwargs &kw) { return createElementOfType("qbrick1mt",args,kw); }
// mpm up
py::object upQuad11(py::args args, py::kwargs &kw) { return createElementOfType("upquad11",args,kw); }
// mpm experimental
py::object q1(py::args args, py::kwargs &kw) { return createElementOfType("q1",args,kw); }
py::object l1(py::args args, py::kwargs &kw) { return createElementOfType("l1",args,kw); }





/*****************************************************
* DofManager
*****************************************************/

py::object createDofManagerOfType(const char*type, py::args args, py::kwargs &kw)
{
    unsigned long number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? PY_CAST(oofem::Domain*,args[1]) : nullptr;
    auto dofMan = oofem::classFactory.createDofManager(type,number,domain);
    if (!dofMan) { OOFEM_RAISE("dofManager: wrong input data"); }
    dofMan->setGlobalNumber(number);
    OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    dofMan->initializeFrom(ir);
    return py::cast(dofMan.release());
}


py::object node(py::args args, py::kwargs kw) { return createDofManagerOfType("node",args,kw); }
py::object hangingnode(py::args args, py::kwargs kw) { return createDofManagerOfType("hangingnode",args,kw); }


/*****************************************************
* GeneralBoundaryCondition
*****************************************************/
// generalBoundaryCondition(aClass,domain=defaultDomain,**kw)
py::object createGeneralBoundaryConditionOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? PY_CAST(oofem::Domain*,args[1]) : nullptr;
    auto bc = oofem::classFactory.createBoundaryCondition(type,number,domain);
    if (!bc) { OOFEM_RAISE("generalBoundaryCondition: wrong input data"); }
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
    oofem::Domain* domain = len(args)>1? PY_CAST(oofem::Domain*,args[1]) : nullptr;
    auto ic = oofem::classFactory.createInitialCondition(type,number,domain);
    if (!ic) { OOFEM_RAISE("initialCondition: wrong input data"); }
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
    oofem::Domain* domain = len(args)>1? PY_CAST(oofem::Domain*,args[1]) : nullptr;
    auto mat = oofem::classFactory.createMaterial(type,number,domain);
    if (!mat) { OOFEM_RAISE("material: wrong input data"); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    mat->initializeFrom(ir);
    return py::cast(mat.release());
}

py::object isoLE(py::args args, py::kwargs kw) { return createMaterialOfType("isole",args,kw); }
py::object idm1(py::args args, py::kwargs kw) { return createMaterialOfType("idm1",args,kw); }
py::object isoHeat(py::args args, py::kwargs kw) { return createMaterialOfType("isoheat",args,kw); }
py::object hydratingConcreteMat(py::args args, py::kwargs kw) { return createMaterialOfType("hydratingconcretemat",args,kw); }
py::object j2mat(py::args args, py::kwargs kw) { return createMaterialOfType("j2mat",args,kw); }
py::object steel1(py::args args, py::kwargs kw) { return createMaterialOfType("steel1",args,kw); }
py::object concreteFcmViscoelastic(py::args args, py::kwargs kw) { return createMaterialOfType("concretefcmviscoelastic",args,kw); }
py::object mps(py::args args, py::kwargs kw) { return createMaterialOfType("mps",args,kw); }

py::object upm(py::args args, py::kwargs kw) { return createMaterialOfType("upm",args,kw); }



/*****************************************************
* CrossSection
*****************************************************/
py::object createCrossSectionOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? PY_CAST(oofem::Domain*,args[1]) : nullptr;
    auto cs = oofem::classFactory.createCrossSection(type,number,domain);
    if (!cs) { OOFEM_RAISE("crossSection: wrong input data"); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    cs->initializeFrom(ir);
    return py::cast(cs.release());
}

py::object simpleCS(py::args args, py::kwargs kw) { return createCrossSectionOfType("simplecs",args,kw); }
py::object simpleTransportCS(py::args args, py::kwargs kw) { return createCrossSectionOfType("simpletransportcs",args,kw); }
py::object dummyCS(py::args args, py::kwargs kw) { return createCrossSectionOfType("dummycs",args,kw); }



/*****************************************************
* LoadTimeFunction
*****************************************************/
py::object createLoadTimeFunctionOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::Domain* domain = len(args)>1? PY_CAST(oofem::Domain*,args[1]) : nullptr;
    auto ltf = oofem::classFactory.createFunction(type,number,domain);
    if (!ltf) { OOFEM_RAISE("function: wrong input data"); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    ltf->initializeFrom(ir);
    return py::cast(ltf.release());
}

py::object peakFunction(py::args args, py::kwargs kw) { return createLoadTimeFunctionOfType("peakfunction",args,kw); }

py::object constantFunction(py::args args, py::kwargs kw) { return createLoadTimeFunctionOfType("constantfunction",args,kw); }

py::object piecewiseLinFunction(py::args args, py::kwargs kw) { return createLoadTimeFunctionOfType("piecewiselinfunction",args,kw); }

py::object usrDefFunction(py::args args, py::kwargs kw) { return createLoadTimeFunctionOfType("usrdefltf",args,kw); }




/*****************************************************
* ExportModule
*****************************************************/
py::object createExportModuleOfType(const char* type, py::args args, py::kwargs kw)
{
    int number = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):0;
    oofem::EngngModel *engngModel = len(args)>1? PY_CAST(oofem::EngngModel*,args[1]) : nullptr;
    auto module = oofem::classFactory.createExportModule(type,number,engngModel);
    if (!module) { OOFEM_RAISE("exportModule: wrong input data"); }
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
    oofem::Domain* domain = len(args)>1? PY_CAST(oofem::Domain*,args[1]) : nullptr;
    std::unique_ptr<Set> setP = std::make_unique<Set>(number, domain);
    if ( !setP ) { OOFEM_RAISE(("Couldn't create set: "+std::to_string(number)).c_str()); }
    oofem::OOFEMTXTInputRecord ir = makeOOFEMTXTInputRecordFrom(kw);
    setP->initializeFrom(ir);
    return py::cast(setP.release());
}

py::object createSet(py::args args, py::kwargs kw) { return createSetOfType("set",args,kw); }


/******************************************************
 * Interpolations
*******************************************************/
py::object createInterpolationOfType(std::string type, py::args args, py::kwargs kw)
{
    if (type == "fei2dquadlin") {
        int i1 = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):1;
        int i2 = len(args)>0?PyLong_AsUnsignedLong(args[1].ptr()):2;
        std::unique_ptr<FEInterpolation> interpol = std::make_unique<FEI2dQuadLin>(i1,i2); 
        return py::cast(interpol.release());
    } else if (type == "fei2dlinelin") {
        int i1 = len(args)>0?PyLong_AsUnsignedLong(args[0].ptr()):1;
        int i2 = len(args)>0?PyLong_AsUnsignedLong(args[1].ptr()):2;
        std::unique_ptr<FEInterpolation> interpol = std::make_unique<FEI2dLineLin>(i1,i2); 
        return py::cast(interpol.release());
    } 
#ifdef __MPM_MODULE    
    else if (type == "linearinterpolation") {
        std::unique_ptr<FEInterpolation> interpol = std::make_unique<LinearInterpolation>(); 
        return py::cast(interpol.release());
    } 
#endif
   return py::none();
}

py::object fei2dquadlin(py::args args, py::kwargs kw) { return createInterpolationOfType("fei2dquadlin",args,kw); }
py::object fei2dlinelin(py::args args, py::kwargs kw) { return createInterpolationOfType("fei2dlinelin",args,kw); }
py::object linearinterpolation(py::args args, py::kwargs kw) { return createInterpolationOfType("linearinterpolation",args,kw); }


/******************************************************
 * Terms
*******************************************************/
#ifdef __MPM_MODULE
#include "prototype2.h"
py::object createTermOfType(std::string type, py::args args, py::kwargs kw)
{
    if (type == "BTSigmaTerm") {
        if (len(args)>2) {
            oofem::Variable * f = PY_CAST(oofem::Variable*,args[0]);
            oofem::Variable * tf = PY_CAST(oofem::Variable*,args[1]);
            oofem::MaterialMode m = PY_CAST(oofem::MaterialMode,args[2]);
            std::unique_ptr<Term> t = std::make_unique<BTSigmaTerm2>(f, tf, m);
            return py::cast(t.release());
        }
    } else if (type == "NTfTerm") {
        if (len(args)>=3) {
            oofem::Variable * f = PY_CAST(oofem::Variable*,args[0]);
            oofem::Variable * tf = PY_CAST(oofem::Variable*,args[1]);
            oofem::MaterialMode m = PY_CAST(oofem::MaterialMode,args[2]);
            std::unique_ptr<Term> t;
            if (len(args)==3) {
                std::unique_ptr<Term> t = std::make_unique<NTfTerm>(f, tf, m);
                return py::cast(t.release());

            } else {
                oofem::FloatArray* flux = PY_CAST(oofem::FloatArray*,args[3]);
                std::unique_ptr<Term> t = std::make_unique<NTfTerm>(f, tf, m, flux);
                return py::cast(t.release());
            }
        }
    }
    return py::none();
}

py::object BTSigma_Term(py::args args, py::kwargs kw) { return createTermOfType("BTSigmaTerm",args,kw); }
py::object NTf_Term(py::args args, py::kwargs kw) { return createTermOfType("NTfTerm",args,kw); }

#endif

/************************************************************
 * Sparse matrices
*************************************************************/
py::object skyline() {
    std::unique_ptr<SparseMtrx> t = std::make_unique<Skyline>();
    return py::cast(t.release()); 
}

/************************************************************
 * Linear sparse solvers
*************************************************************/
py::object ldltFactorization(py::args args, py::kwargs kw) {
    oofem::Domain* domain = len(args)>0? PY_CAST(oofem::Domain*,args[0]) : nullptr;
    oofem::EngngModel *emodel = len(args)>1? PY_CAST(oofem::EngngModel*,args[1]) : nullptr;
    std::unique_ptr<SparseLinearSystemNM> t = std::make_unique<LDLTFactorization>(domain, emodel);
    return py::cast(t.release()); 
}
