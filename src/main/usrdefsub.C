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
 *               Copyright (C) 1993 - 2010   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "usrdefsub.h"
#include "compiler.h" // to supply missing strncasecmp on some platforms
#ifndef __MAKEDEPEND
 #include <string.h>
 #ifdef HAVE_STRINGS_H
  #include <strings.h>
 #endif
#endif


// __OOFEMLIB_MODULE
#include "node.h"
#include "element.h"
#include "engngm.h"
#include "xfemmanager.h"
#include "load.h"
#include "loadtime.h"
#include "material.h"
#include "gaussintegrationrule.h"

#include "sparsemtrx.h"
#include "skyline.h"
#include "skylineu.h"
#include "compcol.h"
#include "dyncompcol.h"
#include "symcompcol.h"
#include "dyncomprow.h"
#include "spoolessparsemtrx.h"
#include "petscsparsemtrx.h"

#include "ldltfact.h"
#include "imlsolver.h"
#include "spoolessolver.h"
#include "petscsolver.h"
#include "dsssolver.h"
#ifdef __PARALLEL_MODE
 #include "fetisolver.h"
#endif

#include "subspaceit.h"
#include "inverseit.h"
#include "slepcsolver.h"

// Nonlinear solvers
#include "nrsolver.h"
#include "nrsolver2.h"
#include "calmls.h"

// materials in oofemlib
#include "dummymaterial.h"

// general loads in OOFEMLIB
#include "linearedgeload.h"
#include "constantedgeload.h"
#include "constantsurfaceload.h"
#include "pointload.h"
#include "prescribedgradient.h"
#include "mixedgradientpressurebc.h"
#include "surfacetensionbc.h"

// ltf in OOFEMLIB
#include "peak.h"
#include "piecewis.h"
#include "piecewisper.h"
#include "heavisideltf.h"
#include "usrdeftimefunct.h"

// export modules
#include "vtkexportmodule.h"
#include "vtkxmlexportmodule.h"

// nodal recovery models
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

//#include "particletopologydescription.h" // Soon

// end __OOFEMLIB_MODULE


#ifdef __SM_MODULE

// Elements of SM module
 #include "truss2d.h"
 #include "trplanstrss.h"
 #include "trplanrot.h"
 #include "trplanrot3d.h"
 #include "libeam2d.h"
 #include "planstrss.h"
 #include "quad1planestrain.h"
 #include "qplanestrain.h"
 #include "qplanestraingrad.h"
 #include "qtrplanestraingrad.h"
 #include "qtrplanestrain.h"
 #include "qtruss1d.h"
 #include "misesmatgrad.h"
 #include "misesmatnl.h"
 #include "qtruss1dgrad.h"
 #include "qspacegrad.h"
 #include "qplanstrss.h"
 #include "qplanestressgrad.h"
 #include "qtrplstr.h"
 #include "qtrplstrgrad.h"
 #include "lspace.h"
 #include "lspacebb.h"
 #include "qspace.h"
 #include "axisymm3d.h"
 #include "q4axisymm.h"
 #include "l4axisymm.h"
 #include "ltrspace.h"
 #include "beam2d.h"
 #include "beam3d.h"
 #include "libeam2dnl.h"
 #include "libeam3dnl.h"
 #include "truss3d.h"
 #include "trplanestrain.h"
 #include "libeam3dnl2.h"
 #include "libeam3d.h"
 #include "libeam3d2.h"
 #include "truss1d.h"
 #include "cct.h"
 #include "cct3d.h"
 #include "tr_shell01.h"
 #include "rershell.h"
 #include "interfaceelem2dquad.h"
 #include "interfaceelement1d.h"
 #include "interfaceelem3dtrlin.h"
 #include "macrolspace.h"
 #include "planstrssxfem.h"
 #include "cohsur3d.h"
 #include "lumpedmasselement.h"
 #include "springelement.h"
 #include "particle.h"
 #include "lattice2d.h"

// iga elements
 #include "igaelements.h"

// Emodels of SM module
 #include "nlinearstatic.h"
 #include "nlineardynamic.h"
 #include "eigenvaluedynamic.h"
 #include "deidynamic.h"
 #include "nldeidynamic.h"
 #include "pnldeidynamic.h"
 #include "diidynamic.h"
 #include "incrementallinearstatic.h"
 #include "linearstability.h"
 #include "plinearstatic.h"
 #include "linearstatic.h"
 #include "stationaryflow.h"
 #include "adaptnlinearstatic.h"
 #include "adaptlinearstatic.h"


// loads of SM module
 #include "structtemperatureload.h"
 #include "structeigenstrainload.h"
 #include "usrdeftempfield.h"
 #include "tf1.h"


// crosssections of SM module
 #include "layeredcrosssection.h"
 #include "fiberedcs.h"


// materials of SM module
 #include "ortholinearelasticmaterial.h"
 #include "perfectlyplasticmaterial.h"
 #include "steel1.h"
 #include "concrete2.h"
 #include "concrete3.h"
 #include "cebfip78.h"
 #include "doublepowerlaw.h"
 #include "b3mat.h"
 #include "b3solidmat.h"
 #include "mps.h"
 #include "j2plasticmaterial.h"
 #include "rcsd.h"
 #include "rcsde.h"
 #include "rcsdnl.h"
 #include "m4.h"
 #include "idm1.h"
 #include "idmnl1.h"
 #include "mazarsmodel.h"
 #include "mazarsmodelnl.h"
 #include "druckerPragerPlasticitySM.h"
 #include "j2mplasticmaterial.h"
 #include "rankinepm.h"
 #include "masonry02.h"
 #include "isointerfacedamage01.h"
 #include "j2mat.h"
 #include "mat_cebfip90.h"
 #include "hellmat.h"
 #include "mdm.h"
 #include "compodamagemat.h"
 #include "micromaterial.h"
 #include "hyperelasticmaterial.h"
 #include "misesmat.h"
 #include "rankinematnl.h"
 #include "rankinematgrad.h"
 #include "trabbonematerial.h"
 #include "trabbonenl.h"
 #include "trabbone3d.h"
 #include "trabboneembed.h"
 #include "trabbonenlembed.h"
 #include "trabbonenl3d.h"
 #include "concretedpm.h"
 #include "concretedpm2.h"
 #include "cohint.h"
 #include "latticedamage2d.h"

 #include "scalarerrorindicator.h"
 #include "zzerrorestimator.h"
 #include "combinedzzsiee.h"
 #include "huertaerrorestimator.h"
 #include "simpleinterfacemat.h"

// export modules
 #include "poiexportmodule.h"
 #include "homexportmodule.h"
 #include "dmexportmodule.h"
 #include "gpexportmodule.h"

// init modules
 #include "gpinitmodule.h"

// nonlocal barriers
 #include "polylinenonlocalbarrier.h"
 #include "symmetrybarrier.h"

// random generators
 #include "localgaussianrandomgenerator.h"
 #include "externalfieldgenerator.h"

// mesher interfaces
 #include "t3dinterface.h"
 #include "targe2interface.h"
 #include "freeminterface.h"
 #include "subdivision.h"

 #include "dss.h"
#endif //__SM_MODULE



#ifdef __TM_MODULE
// Emodels of SM module
 #include "stationarytransportproblem.h"
 #include "nonstationarytransportproblem.h"
 #include "nltransienttransportproblem.h"
 #include "staggeredproblem.h"
// Elements of TM module
 #include "quad1_ht.h"
 #include "tr1_ht.h"
 #include "quadaxisym1_ht.h"
 #include "traxisym1_ht.h"
 #include "brick1_ht.h"
 #include "tetrah1_ht.h"
// materials of TM module
 #include "isoheatmat.h"
 #include "hemotkmat.h"
 #include "hydratingisoheatmat.h"
 #include "hydratinghemomat.h"
 #include "hydratingconcretemat.h"
 #ifdef __CEMHYD_MODULE
  #include "cemhydmat.h"
 #endif
#endif //__TM_MODULE


#ifdef __FM_MODULE
// Emodels
 #include "cbs.h"
 #include "stokesflow.h"
 #include "stokesflowstresshomogenization.h"
// Elements
 #include "tr1_2d_cbs.h"
 #include "tr21stokes.h"
 #include "linesurfacetension.h"
 #include "line2surfacetension.h"
 #include "line2boundaryelement.h"
// materials
 #include "newtonianfluid.h"
 #include "fe2sinteringmaterial.h"
// boundary conditions
 #include "tractionpressurebc.h"

 #include "supg.h"
 #include "tr1_2d_supg.h"
 #include "tr1_2d_supg2.h"
 #include "tr1_2d_supg_axi.h"
 #include "tr1_2d_supg2_axi.h"
 #include "tet1_3d_supg.h"
 #include "tr21_2d_supg.h"
 #include "quad10_2d_supg.h"
 #include "twofluidmaterial.h"
 #include "binghamfluid2.h"

 #include "surfacetensionmaterial.h"

#endif // __FM_Module

// GENERAL
#include "masterdof.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "activedof.h"

// XFEM
#include "enrichmentfunction.h"
#include "enrichmentitem.h"

#ifdef __PARALLEL_MODE
 #include "loadbalancer.h"
 #include "parmetisloadbalancer.h"
#endif

#include <string>

// Comparison operator for strings. Just strcasecmp here?
struct CaseComp
{
    int operator() (std::string a, std::string b) const { return strncasecmp(a.c_str(), b.c_str(), b.length()) < 0; }
};

namespace oofem {

// Template to wrap constructors into functions
template < typename T > Element* elemCreator( int n, Domain *d ) { return new T(n, d); }
std::map< std::string, Element*(*)(int,Domain*), CaseComp > elemList;

Element *CreateUsrDefElementOfType(const char *aClass, int number, Domain *domain)
{
    if (elemList.size() == 0) {
#ifdef __SM_MODULE
        elemList["planestress2dxfem"]  = elemCreator< PlaneStress2dXfem >;
        elemList["planestress2d"]      = elemCreator< PlaneStress2d >;
        elemList["quad1planestrain"]   = elemCreator< Quad1PlaneStrain >;
        elemList["qtrplanestraingrad"] = elemCreator< QTrPlaneStrainGrad >;
        elemList["qtrplanestrain"]     = elemCreator< QTrPlaneStrain >;
        elemList["qplanestraingrad"]   = elemCreator< QPlaneStrainGrad >;
        elemList["qplanestrain"]       = elemCreator< QPlaneStrain >;
        elemList["trplanestress2d"]    = elemCreator< TrPlaneStress2d >;
        elemList["trplanestrrot3d"]    = elemCreator< TrPlaneStrRot3d >;
        elemList["trplanestrrot"]      = elemCreator< TrPlaneStrRot >;
        elemList["qplanestressgrad"]   = elemCreator< QPlaneStressGrad >;
        elemList["qplanestress2d"]     = elemCreator< QPlaneStress2d >;
        elemList["qtrplstrgrad"]       = elemCreator< QTrPlaneStressGrad >;
        elemList["qtrplstr"]           = elemCreator< QTrPlaneStress2d >;
        elemList["axisymm3d"]          = elemCreator< Axisymm3d >;
        elemList["q4axisymm"]          = elemCreator< Q4Axisymm >;
        elemList["l4axisymm"]          = elemCreator< L4Axisymm >;
        elemList["lspacebb"]           = elemCreator< LSpaceBB >;
        elemList["lspace"]             = elemCreator< LSpace >;
        elemList["qspacegrad"]         = elemCreator< QSpaceGrad >;
        elemList["qspace"]             = elemCreator< QSpace >;
        elemList["cctplate3d"]         = elemCreator< CCTPlate3d >;
        elemList["cctplate"]           = elemCreator< CCTPlate >;
        //elemList["ltrspaceec"]       = elemCreator< LTRSpaceWithEmbeddedCrack >;
        elemList["ltrspace"]           = elemCreator< LTRSpace >;
        elemList["truss2d"]            = elemCreator< Truss2d >;
        elemList["rershell"]           = elemCreator< RerShell >;
        elemList["tr_shell01"]         = elemCreator< TR_SHELL01 >;
        elemList["beam2d"]             = elemCreator< Beam2d >;
        elemList["beam3d"]             = elemCreator< Beam3d >;
        elemList["libeam2dNL"]         = elemCreator< LIBeam2dNL >;
        elemList["libeam2d"]           = elemCreator< LIBeam2d >;
        elemList["libeam3dnl2"]        = elemCreator< LIBeam3dNL2 >;
        elemList["libeam3dNL"]         = elemCreator< LIBeam3dNL >;
        elemList["truss3d"]            = elemCreator< Truss3d >;
        elemList["trplanestrain"]      = elemCreator< TrPlaneStrain >;
        elemList["libeam3d2"]          = elemCreator< LIBeam3d2 >;
        elemList["libeam3d"]           = elemCreator< LIBeam3d >;
        elemList["qtruss1dgrad"]       = elemCreator< QTruss1dGrad >;
        elemList["qtruss1d"]           = elemCreator< QTruss1d >;
        elemList["truss1d"]            = elemCreator< Truss1d >;
        elemList["interface2dquad"]    = elemCreator< InterfaceElem2dQuad >;
        elemList["interface3dtrlin"]   = elemCreator< InterfaceElement3dTrLin >;
        elemList["interface1d"]        = elemCreator< InterfaceElem1d >;
        elemList["macrolspace"]        = elemCreator< MacroLSpace >;
        elemList["lumpedmass"]         = elemCreator< LumpedMassElement >;
        elemList["spring"]             = elemCreator< SpringElement >;
        elemList["cohsur3d"]           = elemCreator< CohesiveSurface3d >;
        elemList["bsplineplanestresselement"] = elemCreator< BsplinePlaneStressElement >;
        elemList["nurbsplanestresselement"]   = elemCreator< NURBSPlaneStressElement >;
        elemList["tsplineplanestresselement"] = elemCreator< TSplinePlaneStressElement >;
        elemList["nurbs3delement"]            = elemCreator< NURBSSpace3dElement >;
        elemList["lattice2d"]            = elemCreator< Lattice2d >;

#endif
#ifdef __FM_MODULE
        elemList["tr1cbs"]         = elemCreator< TR1_2D_CBS >;
        elemList["tr1supgaxi"]     = elemCreator< TR1_2D_SUPG_AXI >;
        elemList["tr1supg2axi"]    = elemCreator< TR1_2D_SUPG2_AXI >;
        elemList["tr1supg2"]       = elemCreator< TR1_2D_SUPG2 >;
        elemList["tr1supg"]        = elemCreator< TR1_2D_SUPG >;
        elemList["tet1supg"]       = elemCreator< Tet1_3D_SUPG >;
        elemList["tr21supg"]       = elemCreator< TR21_2D_SUPG >;
        elemList["quad1supg"]      = elemCreator< Quad10_2D_SUPG >;
        elemList["tr21stokes"]     = elemCreator< Tr21Stokes >;
        elemList["line2boundary"]  = elemCreator< Line2BoundaryElement >;
        elemList["linesurfacetension"]  = elemCreator< LineSurfaceTension >;
        elemList["line2surfacetension"] = elemCreator< Line2SurfaceTension >;
#endif
#ifdef __TM_MODULE
        elemList["quad1ht"]        = elemCreator< Quad1_ht >;
        elemList["quad1hmt"]       = elemCreator< Quad1_hmt >;
        elemList["quadaxisym1ht"]  = elemCreator< QuadAxisym1_ht >;
        elemList["quadaxisym1hmt"] = elemCreator< QuadAxisym1_hmt >;
        elemList["brick1ht"]       = elemCreator< Brick1_ht >;
        elemList["brick1hmt"]      = elemCreator< Brick1_hmt >;
        elemList["tetrah1ht"]      = elemCreator< Tetrah1_ht >;
        elemList["tetrah1hmt"]     = elemCreator< Tetrah1_hmt >;
        elemList["tr1ht"]          = elemCreator< Tr1_ht >;
        elemList["traxisym1ht"]    = elemCreator< TrAxisym1_ht >;
#endif
    }
    return (elemList.count(aClass) == 1) ? elemList[aClass](number, domain) : NULL;
}


template < typename T > DofManager* dofmanCreator( int n, Domain *d ) { return new T(n, d); }
std::map< std::string, DofManager*(*)(int,Domain*), CaseComp > dofmanList;

DofManager *CreateUsrDefDofManagerOfType(const char *aClass, int number, Domain *domain)
{
    if ( dofmanList.size() == 0 ) {
#ifdef __SM_MODULE
        dofmanList["particle"] = dofmanCreator< Particle >;
#endif //__SM_MODULE
    }
    return (dofmanList.count(aClass) == 1) ? dofmanList[aClass](number, domain) : NULL;
}


template < typename T > CrossSection* csCreator( int n, Domain *d ) { return new T(n, d); }
std::map< std::string, CrossSection*(*)(int,Domain*), CaseComp > csList;

CrossSection *CreateUsrDefCrossSectionOfType(const char *aClass, int number, Domain *domain)
{
    if ( dofmanList.size() == 0 ) {
#ifdef __SM_MODULE
        csList["layeredcs"] = csCreator< LayeredCrossSection >;
        csList["fiberedcs"] = csCreator< FiberedCrossSection >;
#endif //__SM_MODULE
    }
    return (csList.count(aClass) == 1) ? csList[aClass](number, domain) : NULL;
}


// Template to wrap constructors into functions
template < typename T > EngngModel* engngCreator( int n, EngngModel *m ) { return ( new T(n, m) ); }
std::map< std::string, EngngModel*(*)(int, EngngModel*), CaseComp > engngList;

EngngModel *CreateUsrDefEngngModelOfType(const char *aClass, int number, EngngModel *master)
{
    if (engngList.size() == 0) {
#ifdef __SM_MODULE
        engngList["linearstatic"]       = engngCreator< LinearStatic >;
        engngList["stationaryflow"]     = engngCreator< StationaryFlow >;
        engngList["eigenvaluedynamic"]  = engngCreator< EigenValueDynamic >;
        engngList["nonlinearstatic"]    = engngCreator< NonLinearStatic >;
        engngList["nonlineardynamic"]   = engngCreator< NonLinearDynamic >;
        engngList["nldeidynamic"]       = engngCreator< NlDEIDynamic >;
        engngList["pnldeidynamic"]      = engngCreator< PNlDEIDynamic >;
        engngList["deidynamic"]         = engngCreator< DEIDynamic >;
        engngList["diidynamic"]         = engngCreator< DIIDynamic >;
        engngList["incrlinearstatic"]   = engngCreator< IncrementalLinearStatic >;
        engngList["linearstability"]    = engngCreator< LinearStability >;
        engngList["adaptnlinearstatic"] = engngCreator< AdaptiveNonLinearStatic >;
        engngList["adaptlinearstatic"]  = engngCreator< AdaptiveLinearStatic >;
 #ifdef __PARALLEL_MODE
        engngList["plinearstatic"]      = engngCreator< PLinearStatic >;
 #endif
#endif //__SM_MODULE
#ifdef __TM_MODULE
        engngList["stationaryproblem"]              = engngCreator< StationaryTransportProblem >;
        engngList["nonstationaryproblem"]           = engngCreator< NonStationaryTransportProblem >;
        engngList["nltransienttransportproblem"]    = engngCreator< NLTransientTransportProblem >;
        engngList["staggeredproblem"]               = engngCreator< StaggeredProblem >;
#endif //__TM_MODULE
#ifdef __FM_MODULE
        engngList["cbs"]        = engngCreator< CBS >;
        engngList["supg"]       = engngCreator< SUPG >;
        engngList["stokesflow"] = engngCreator< StokesFlow >;
        engngList["stokesflowstresshomogenization"] = engngCreator< StokesFlowStressHomogenization >;
#endif //__FM_MODULE
    }
    return (engngList.count(aClass) == 1) ? engngList[aClass](number, master) : NULL;
}


template < typename T > GeneralBoundaryCondition* bcCreator( int n, Domain *d ) { return new T(n, d); }
std::map< std::string, GeneralBoundaryCondition*(*)(int, Domain*), CaseComp > bcList;

GeneralBoundaryCondition *CreateUsrDefBoundaryConditionOfType(const char *aClass, int number, Domain *domain)
{
    if (bcList.size() == 0) {
        bcList["prescribedgradient"]    = bcCreator< PrescribedGradient >;
        bcList["mixedgradientpressure"] = bcCreator< MixedGradientPressureBC >;
        bcList["linearedgeload"]        = bcCreator< LinearEdgeLoad >;
        bcList["constantedgeload"]      = bcCreator< ConstantEdgeLoad >;
        bcList["constantsurfaceload"]   = bcCreator< ConstantSurfaceLoad >;
        bcList["pointload"]             = bcCreator< PointLoad >;
        bcList["surfacetension"]        = bcCreator< SurfaceTensionBoundaryCondition >;
#ifdef __SM_MODULE
        bcList["structtemperatureload"] = bcCreator< StructuralTemperatureLoad >;
        bcList["structeigenstrainload"] = bcCreator< StructuralEigenstrainLoad >;
        bcList["usrdeftempfield"]       = bcCreator< UserDefinedTemperatureField >;
        bcList["tf1"]                   = bcCreator< TF1 >;
#endif //__SM_MODULE
#ifdef __FM_MODULE
        bcList["prescribedtractionpressurebc"] = bcCreator< TractionPressureBC >;
#endif
    }
    return (bcList.count(aClass) == 1) ? bcList[aClass](number, domain) : NULL;
}


template < typename T > LoadTimeFunction* ltfCreator( int n, Domain *d ) { return new T(n, d); }
std::map< std::string, LoadTimeFunction*(*)(int, Domain*), CaseComp > ltfList;

LoadTimeFunction *CreateUsrDefLoadTimeFunctionOfType(const char *aClass, int number, Domain *domain)
{
    if (ltfList.size() == 0) {
        ltfList["peakfunction"]                 = ltfCreator< PeakFunction >;
        ltfList["piecewiselinfunction"]         = ltfCreator< PiecewiseLinFunction >;
        ltfList["periodicpiecewiselinfunction"] = ltfCreator< PeriodicPiecewiseLinFunction >;
        ltfList["heavisideltf"]                 = ltfCreator< HeavisideLTF >;
        ltfList["usrdefltf"]                    = ltfCreator< UserDefinedLoadTimeFunction >;
    }
    return (ltfList.count(aClass) == 1) ? ltfList[aClass](number, domain) : NULL;
}


template < typename T > Material* matCreator( int n, Domain *d ) { return new T(n, d); }
std::map< std::string, Material*(*)(int,Domain*), CaseComp > matList;

Material *CreateUsrDefMaterialOfType(const char *aClass, int number, Domain *domain)
{
    if (matList.size() == 0) {
#ifdef __SM_MODULE
        matList["dummymat"]         = matCreator< DummyMaterial >;
        matList["orthole"]          = matCreator< OrthotropicLinearElasticMaterial >;
        matList["steel1"]           = matCreator< Steel1 >;
        matList["concrete2"]        = matCreator< Concrete2 >;
        matList["concrete3"]        = matCreator< Concrete3 >;
        matList["cebfip78"]         = matCreator< CebFip78Material >;
        matList["doublepowerlaw"]   = matCreator< DoublePowerLawMaterial >;
        matList["b3mat"]            = matCreator< B3Material >;
        matList["b3solidmat"]       = matCreator< B3SolidMaterial >;
        matList["mps"]              = matCreator< MPSMaterial >;
        matList["j2mat"]            = matCreator< J2plasticMaterial >;
        matList["rcsdnl"]           = matCreator< RCSDNLMaterial >;
        matList["rcsde"]            = matCreator< RCSDEMaterial >;
        matList["rcsd"]             = matCreator< RCSDMaterial >;
        matList["microplane_m4"]    = matCreator< M4Material >;
        matList["idm1"]             = matCreator< IsotropicDamageMaterial1 >;
        matList["idmnl1"]           = matCreator< IDNLMaterial >;
        matList["mazarsmodelnl"]    = matCreator< MazarsNLMaterial >;
        matList["mazarsmodel"]      = matCreator< MazarsMaterial >;
        matList["druckerprager"]    = matCreator< DruckerPragerPlasticitySM >;
        matList["j2mmat"]           = matCreator< J2MPlasticMaterial >;
        matList["rankine"]          = matCreator< RankinePlasticMaterial >;
        matList["masonry02"]        = matCreator< Masonry02 >;
        matList["isointrfdm01"]     = matCreator< IsoInterfaceDamageMaterial >;
        matList["j22mat"]           = matCreator< J2Mat >;
        matList["cebfipslip90"]     = matCreator< CebFipSlip90Material >;
        matList["hellmat"]          = matCreator< HellmichMaterial >;
        matList["mdm"]              = matCreator< MDM >;
        matList["compdammat"]       = matCreator< CompoDamageMat >;
        matList["micromat"]         = matCreator< MicroMaterial >;
        matList["hyperelmat"]       = matCreator< HyperElasticMaterial >;
        matList["misesmatgrad"]     = matCreator< MisesMatGrad >;
        matList["misesmatnl"]       = matCreator< MisesMatNl >;
        matList["misesmat"]         = matCreator< MisesMat >;
        matList["rankmatgrad"]      = matCreator< RankineMatGrad >;
        matList["rankmatnl"]        = matCreator< RankineMatNl >;
        matList["rankmat"]          = matCreator< RankineMat >;
        matList["trabbonenl3d"]     = matCreator< TrabBoneNL3D >;
        matList["trabboneembed"]    = matCreator< TrabBoneEmbed >;
        matList["trabbonenlembed"]  = matCreator< TrabBoneNLEmbed >;
        matList["trabbonenl"]       = matCreator< TrabBoneNL >;
        matList["trabbone3d"]       = matCreator< TrabBone3D >;
        matList["trabbone"]         = matCreator< TrabBoneMaterial >;
        matList["concretedpm"]      = matCreator< ConcreteDPM >;
        matList["concreteidm"]      = matCreator< ConcreteDPM >; // for compatibility with old inputfiles
        matList["cohint"]           = matCreator< CohesiveInterfaceMaterial >;
        matList["simpleintermat"]   = matCreator< SimpleInterfaceMaterial >;
        matList["con2dpm"]          = matCreator< ConcreteDPM2 >;
        matList["latticedamage2d"]  = matCreator< LatticeDamage2d >;
#endif //__SM_MODULE
#ifdef __TM_MODULE
        matList["isoheat"]              = matCreator< IsotropicHeatTransferMaterial >;
        matList["hemotk"]               = matCreator< HeMoTKMaterial >;
        matList["hydratingconcretemat"] = matCreator< HydratingConcreteMat >;
#ifdef __CEMHYD_MODULE
        matList["cemhydmat"]            = matCreator< CemhydMat >;
 #endif //__CEMHYD_MODULE
#endif //__TM_MODULE
#if defined( __SM_MODULE ) && defined ( __TM_MODULE )
        matList["hisoheat"] = matCreator< HydratingIsoHeatMaterial >;
        matList["hhemotk"] = matCreator< HydratingHeMoMaterial >;
#endif //defined(__SM_MODULE) && defined (__TM_MODULE)
#ifdef __FM_MODULE
        matList["newtonianfluid"]       = matCreator< NewtonianFluidMaterial >;
        matList["twofluidmat"]          = matCreator< TwoFluidMaterial >;
        matList["binghamfluid2"]        = matCreator< BinghamFluidMaterial2 >;
        matList["binghamfluid"]         = matCreator< BinghamFluidMaterial2 >;
        matList["fe2sinteringmaterial"] = matCreator< FE2SinteringMaterial >;
        matList["surfacetension"]       = matCreator< SurfaceTensionMaterial >;
#endif // __FM_MODULE
    }
    return (matList.count(aClass) == 1) ? matList[aClass](number, domain) : NULL;
}


template < typename T > TopologyDescription* topologyCreator( Domain *d ) { return new T(d); }
std::map< std::string, TopologyDescription*(*)(Domain*), CaseComp > topologyList;

TopologyDescription *CreateUsrDefTopologyOfType(const char *aClass, Domain *domain)
{
    if (topologyList.size() == 0) {
        //topologyList["particletopology"] = topologyCreator< ParticleTopologyDescription >; // Soon..
    }
    return (topologyList.count(aClass) == 1) ? topologyList[aClass](domain) : NULL;
}

SparseMtrx *CreateUsrDefSparseMtrx(SparseMtrxType type)
{
    SparseMtrx *answer = NULL;

    if ( type == SMT_Skyline ) {
        answer = new Skyline();
    } else if ( type == SMT_SkylineU ) {
        answer = new SkylineUnsym();
    }

#ifdef __IML_MODULE
    else if ( type == SMT_CompCol ) {
        answer = new CompCol();
    } else if ( type == SMT_DynCompCol ) {
        answer = new DynCompCol();
    } else if ( type == SMT_SymCompCol ) {
        answer = new SymCompCol();
    } else if ( type == SMT_DynCompRow ) {
        answer = new DynCompRow();
    }
#endif
#ifdef __SPOOLES_MODULE
    else if ( type == SMT_SpoolesMtrx ) {
        answer = new SpoolesSparseMtrx();
    }
#endif
#ifdef __PETSC_MODULE
    else if ( type == SMT_PetscMtrx ) {
        answer = new PetscSparseMtrx();
    }
#endif
#ifdef __DSS_MODULE
    else if ( type == SMT_DSS_sym_LDL ) {
        answer = new DSSMatrix(DSSMatrix :: sym_LDL);
    } else if ( type == SMT_DSS_sym_LL ) {
        answer = new DSSMatrix(DSSMatrix :: sym_LL);
    } else if ( type == SMT_DSS_unsym_LU ) {
        answer = new DSSMatrix(DSSMatrix :: unsym_LU);
    }
#endif
    else {
        OOFEM_ERROR("CreateUsrDefSparseMtrx: Unknown mtrx type\n");
    }

    return answer;
}

SparseLinearSystemNM *CreateUsrDefSparseLinSolver(LinSystSolverType st, int i, Domain *d, EngngModel *m)
{
    SparseLinearSystemNM *nm = NULL;
    if ( st == ST_Direct ) {
        nm = ( SparseLinearSystemNM * ) new LDLTFactorization(i, d, m);
        return nm;
    } else if ( st == ST_IML ) {
        nm = ( SparseLinearSystemNM * ) new IMLSolver(i, d, m);
        return nm;
    } else if ( st == ST_Spooles ) {
        nm = ( SparseLinearSystemNM * ) new SpoolesSolver(i, d, m);
        return nm;
    } else if ( st == ST_Petsc ) {
        nm = ( SparseLinearSystemNM * ) new PetscSolver(i, d, m);
        return nm;

#ifdef __PARALLEL_MODE
    } else if ( st == ST_Feti ) {
        nm = ( SparseLinearSystemNM * ) new FETISolver(i, d, m);
        return nm;

#endif
#ifdef __DSS_MODULE
    } else if ( st == ST_DSS ) {
        nm = ( SparseLinearSystemNM * ) new DSSSolver(i, d, m);
        return nm;

#endif
    } else {
        OOFEM_ERROR("CreateUsrDefSparseLinSolver: Unknown solver type\n");
    }

    return nm;
}

SparseGeneralEigenValueSystemNM *CreateUsrDefGeneralizedEigenValueSolver(GenEigvalSolverType st, int i, Domain *d, EngngModel *m)
{
    SparseGeneralEigenValueSystemNM *nm = NULL;
    if ( st == GES_SubspaceIt ) {
        nm = ( SparseGeneralEigenValueSystemNM * ) new SubspaceIteration(i, d, m);
        return nm;
    } else if ( st == GES_InverseIt ) {
        nm = ( SparseGeneralEigenValueSystemNM * ) new InverseIteration(i, d, m);
        return nm;
    } else if ( st == GES_SLEPc ) {
        nm = ( SparseGeneralEigenValueSystemNM * ) new SLEPcSolver(i, d, m);
        return nm;
    } else {
        OOFEM_ERROR("CreateUsrDefGeneralizedEigenValueSolver: Unknown solver type\n");
    }

    return nm;
}


template < typename T > SparseNonLinearSystemNM* nonlinCreator( int n, Domain *d, EngngModel *m, EquationID eid ) { return ( new T(n, d, m, eid) ); }
std::map< std::string, SparseNonLinearSystemNM*(*)(int,Domain*,EngngModel*,EquationID), CaseComp > nonlinList;

SparseNonLinearSystemNM *CreateUsrDefNonLinearSolver(const char *aClass, int number, Domain *d, EngngModel *emodel, EquationID eid)
{
    if ( nonlinList.size() == 0 ) {
        //nonlinList["snes"]       = nonlinCreator< PETScSNES >;
        nonlinList["nrsolver"]   = nonlinCreator< NRSolver >;
        nonlinList["nrsolver2"]  = nonlinCreator< NRSolver2 >;
        nonlinList["calm"]       = nonlinCreator< CylindricalALM >;
    }
    return (nonlinList.count(aClass) == 1) ? nonlinList[aClass](number, d, emodel, eid) : NULL;
}

ErrorEstimator *CreateUsrDefErrorEstimator(ErrorEstimatorType type, int number, Domain *d)
{
    ErrorEstimator *answer = NULL;

#ifdef __SM_MODULE
    if ( type == EET_SEI ) {
        answer = new ScalarErrorIndicator(number, d);
    } else if ( type == EET_ZZEE ) {
        answer = new ZZErrorEstimator(number, d);
    } else if ( type == EET_CZZSI ) {
        answer = new CombinedZZSIErrorEstimator(number, d);
    } else if ( type == EET_HEE ) {
        answer = new HuertaErrorEstimator(number, d);
    }

#endif //__SM_MODULE

    if ( answer == NULL ) {
        OOFEM_ERROR("CreateUsrDefErrorEstimator: Unknown error estimator type\n");
    }

    return answer;
}


template < typename T > ExportModule* exportCreator( int n, EngngModel *e ) { return ( new T(n, e) ); }
std::map< std::string, ExportModule*(*)(int,EngngModel*), CaseComp > exportList;

ExportModule *CreateUsrDefExportModuleOfType(const char *aClass, int number, EngngModel *emodel)
{
    if (exportList.size() == 0) {
        exportList["vtkxml"]    = exportCreator< VTKXMLExportModule >;
        exportList["vtk"]       = exportCreator< VTKExportModule >;
#ifdef __SM_MODULE
        exportList["poi"]       = exportCreator< POIExportModule >;
        exportList["hom"]       = exportCreator< HOMExportModule >;
        exportList["dm"]        = exportCreator< DofManExportModule >;
        exportList["gp"]        = exportCreator< GPExportModule >;
#endif //__SM_MODULE
    }
    return (exportList.count(aClass) == 1) ? exportList[aClass](number, emodel) : NULL;
}


template < typename T > InitModule* initCreator( int n, EngngModel *e ) { return ( new T(n, e) ); }
std::map< std::string, InitModule*(*)(int,EngngModel*), CaseComp > initList;

InitModule *CreateUsrDefInitModuleOfType(const char *aClass, int number, EngngModel *emodel)
{
    if ( initList.size() == 0 ) {
#ifdef __SM_MODULE
        initList["gpinitmodule"] = initCreator< GPInitModule >;
#endif //__SM_MODULE
    }
    return (initList.count(aClass) == 1) ? initList[aClass](number, emodel) : NULL;
}


template < typename T > NonlocalBarrier* barrierCreator( int n, Domain *d ) { return ( new T(n, d) ); }
std::map< std::string, NonlocalBarrier*(*)(int,Domain*), CaseComp > barrierList;

NonlocalBarrier *CreateUsrDefNonlocalBarrierOfType(const char *aClass, int number, Domain *domain)
{
    if ( barrierList.size() == 0 ) {
#ifdef __SM_MODULE
        barrierList["polylinebarrier"] = barrierCreator< PolylineNonlocalBarrier >;
        barrierList["symmetrybarrier"] = barrierCreator< SymmetryBarrier >;
#endif //__SM_MODULE
    }
    return (barrierList.count(aClass) == 1) ? barrierList[aClass](number, domain) : NULL;
}


template < typename T > RandomFieldGenerator* randomFieldCreator( int n, Domain *d ) { return ( new T(n, d) ); }
std::map< std::string, RandomFieldGenerator*(*)(int,Domain*), CaseComp > randomFieldList;

RandomFieldGenerator *CreateUsrDefRandomFieldGenerator(const char *aClass, int number, Domain *domain)
{
    if ( randomFieldList.size() == 0 ) {
#ifdef __SM_MODULE
        randomFieldList["localgaussrandomgenerator"] = randomFieldCreator< LocalGaussianRandomGenerator >;
	randomFieldList["externalfieldgenerator"] = randomFieldCreator< ExternalFieldGenerator >;
#endif
    }
    return (randomFieldList.count(aClass) == 1) ? randomFieldList[aClass](number, domain) : NULL;
}


IntegrationRule *CreateUsrDefIRuleOfType(classType type, int number, Element *e)
{
    IntegrationRule *answer = NULL;
    if ( type == GaussIntegrationRuleClass ) {
        answer = new GaussIntegrationRule(number, e);
    }

    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefIRuleOfType: Unknown integration rule type [%d]", type);
    }

    return answer;
}


Element *CreateUsrDefElementOfType(classType type, int number, Domain *domain)
{
    Element *answer = NULL;

#ifdef __SM_MODULE
    if ( type == PlaneStress2dClass ) {
        answer = new PlaneStress2d(number, domain);
    } else if ( type == TrPlaneStress2dClass ) {
        answer = new TrPlaneStress2d(number, domain);
    } else if ( type == LTRSpaceClass ) {
        answer = new LTRSpace(number, domain);
    } else if ( type == TrPlaneStrainClass ) {
        answer = new TrPlaneStrain(number, domain);
    }

#endif
#ifdef __FM_MODULE
    if ( type == Tr21StokesElementClass ) {
        answer = new Tr21Stokes(number, domain);
    } else if ( type == LineSurfaceTensionElementClass ) {
        answer = new LineSurfaceTension(number, domain);
    } else if ( type == Line2SurfaceTensionElementClass ) {
        answer = new Line2SurfaceTension(number, domain);
    } else if ( type == Line2BoundaryElementClass ) {
        answer = new Line2BoundaryElement(number, domain);
    }

#endif

    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefElementOfType: Unknown element type [%d]", type);
    }

    return answer;
}

DofManager *CreateUsrDefDofManagerOfType(classType type, int number, Domain *domain)
{
    DofManager *answer = NULL;
    if ( type == NodeClass ) {
        answer = new Node(number, domain);
    }

    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefDofManagerOfType: Unknown dofman type [%d]", type);
    }

    return answer;
}

Dof *CreateUsrDefDofOfType(classType type, int number, DofManager *dman)
{
    Dof *answer = NULL;
    if  ( type == MasterDofClass ) {
        answer = new MasterDof(number, dman);
    } else if ( type == SimpleSlaveDofClass ) {
        answer = new SimpleSlaveDof(number, dman);
    } else if ( type == SlaveDofClass ) {
        answer = new SlaveDof(number, dman);
    } else if ( type == ActiveDofClass ) {
        answer = new ActiveDof(number, dman);
    }

    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefDofOfType: Unknown dof type [%d]", type);
    }

    return answer;
}

MaterialMappingAlgorithm *CreateUsrDefMaterialMappingAlgorithm(MaterialMappingAlgorithmType type)
{
    MaterialMappingAlgorithm *answer = NULL;

#ifdef __SM_MODULE
    if ( type == MMA_ClosestPoint ) {
        answer = new MMAClosestIPTransfer();
    } else if ( type == MMA_LeastSquareProjection ) {
        answer = new MMALeastSquareProjection();
    } else if ( type == MMA_ShapeFunctionProjection ) {
        answer = new MMAShapeFunctProjection();
    }

#endif
    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefMaterialMappingAlgorithm: Unknown mma type [%d]", type);
    }

    return answer;
}

MesherInterface *CreateUsrDefMesherInterface(MeshPackageType type, Domain *d)
{
    MesherInterface *answer = NULL;
#ifdef __SM_MODULE
    if ( type == MPT_T3D ) {
        answer = new T3DInterface(d);
    } else if ( type == MPT_TARGE2 ) {
        answer = new Targe2Interface(d);
    } else if ( type == MPT_FREEM ) {
        answer = new FreemInterface(d);
    } else if ( type == MPT_SUBDIVISION ) {
        answer = new Subdivision(d);
    }

#endif
    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefMesherInterface: Unknown MI type [%d]", type);
    }

    return answer;
}


template < typename T > EnrichmentItem* enrichItemCreator( int n, XfemManager *x, Domain *d ) { return new T(n, x, d); }
std::map< std::string, EnrichmentItem*(*)(int,XfemManager*,Domain*), CaseComp > enrichItemList;

EnrichmentItem *CreateUsrDefEnrichmentItem(const char *aClass, int number, XfemManager *xm, Domain *domain)
{
    if (enrichItemList.size() == 0) {
        enrichItemList["cracktip"]      = enrichItemCreator< CrackTip >;
        enrichItemList["crackinterior"] = enrichItemCreator< CrackInterior >;
        enrichItemList["inclusion"]     = enrichItemCreator< Inclusion >;
    }
    return (enrichItemList.count(aClass) == 1) ? enrichItemList[aClass](number, xm, domain) : NULL;
}


template < typename T > EnrichmentFunction* enrichFuncCreator( int n, Domain *d ) { return new T(n, d); }
std::map< std::string, EnrichmentFunction*(*)(int,Domain*), CaseComp > enrichFuncList;

EnrichmentFunction *CreateUsrDefEnrichmentFunction(const char *aClass, int number, Domain *domain)
{
    if (enrichFuncList.size() == 0) {
        enrichFuncList["discontinuousfunction"] = enrichFuncCreator< DiscontinuousFunction >;
        enrichFuncList["branchfunction"] = enrichFuncCreator< BranchFunction >;
        enrichFuncList["rampfunction"] = enrichFuncCreator< RampFunction >;
    }
    return (enrichFuncList.count(aClass) == 1) ? enrichFuncList[aClass](number, domain) : NULL;
}


template < typename T > BasicGeometry* geometryCreator() { return new T(); }
std::map< std::string, BasicGeometry*(*)(), CaseComp > geometryList;

BasicGeometry *CreateUsrDefGeometry(const char *aClass)
{
    if (geometryList.size() == 0) {
        geometryList["line"] = geometryCreator< Line >;
        geometryList["circle"] = geometryCreator< Circle >;
    }
    return (geometryList.count(aClass) == 1) ? geometryList[aClass]() : NULL;
}

Patch *CreateUsrDefPatch(Patch :: PatchType ptype, Element *e)
{
    Patch *answer = NULL;
    if ( ptype == Patch :: PT_TrianglePatch ) {
        answer = new TrianglePatch(e);
    } else {
        OOFEM_ERROR2("CreateUsrDefPatch: Unknown PatchType [%d]", ptype);
    }

    return answer;
}

NodalRecoveryModel *
CreateUsrDefNodalRecoveryModel(NodalRecoveryModel :: NodalRecoveryModelType type, Domain *d)
{
    NodalRecoveryModel *answer = NULL;
    if ( type == NodalRecoveryModel :: NRM_NodalAveraging ) {
        answer = new NodalAveragingRecoveryModel(d);
    } else if ( type == NodalRecoveryModel :: NRM_ZienkiewiczZhu ) {
        answer = new ZZNodalRecoveryModel(d);
    } else if ( type == NodalRecoveryModel :: NRM_SPR ) {
        answer = new SPRNodalRecoveryModel(d);
    } else {
        OOFEM_ERROR2("CreateUsrDefNodalRecoveryModel: unsupported NodalRecoveryModelType [%d]", type);
    }

    return answer;
}


#ifdef __PARALLEL_MODE
LoadBalancerMonitor *CreateUsrDefLoadBalancerMonitorOfType(classType type, EngngModel *e)
{
    LoadBalancerMonitor *answer = NULL;
    if ( type == WallClockLoadBalancerMonitorClass ) {
        answer = new WallClockLoadBalancerMonitor(e);
    }

    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefLoadBalancerMonitorOfType: Unknown type [%d]", type);
    }

    return answer;
}

LoadBalancer *CreateUsrDefLoadBalancerOfType(classType type, Domain *d)
{
    LoadBalancer *answer = NULL;
    if ( type == ParmetisLoadBalancerClass ) {
        answer = new ParmetisLoadBalancer(d);
    }

    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefLoadBalancerOfType: Unknown type [%d]", type);
    }

    return answer;
}
#endif
} // end namespace oofem
