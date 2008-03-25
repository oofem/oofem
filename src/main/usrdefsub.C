/* $Header: /home/cvs/bp/oofem/main/src/usrdefsub.C,v 1.8.4.1 2004/04/05 15:19:41 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

// file: usrdefsub.C

#include "usrdefsub.h"
#ifndef __MAKEDEPEND
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#endif

#include "node.h"
#include "element.h"
#include "engngm.h"
#include "load.h"
#include "loadtime.h"
#include "material.h"

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



#ifdef __SM_MODULE

// Elements of SM module
#include "truss2d.h"
#include "trplanstrss.h"
#include "trplanrot.h"
#include "libeam2d.h"
#include "planstrss.h"
#include "quad1planestrain.h"
#include "qplanstrss.h"
#include "qtrplstr.h"
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
#include "rershell.h"
#include "interfaceelem2dquad.h"
#include "interfaceelement1d.h"

// Emodels of SM module
#include "nlinearstatic.h"
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
#include "linearedgeload.h"
#include "constantedgeload.h"
#include "constantsurfaceload.h"
#include "usrdeftempfield.h"
#include "tf1.h"
#include "pointload.h"

// ltf of SM module
#include "peak.h"
#include "piecewis.h"
#include "piecewisper.h"
#include "heavisideltf.h"
#include "usrdeftimefunct.h" 

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

#include "scalarerrorindicator.h"
#include "zzerrorestimator.h"
#include "combinedzzsiee.h"
#include "huertaerrorestimator.h"

// export modules
#include "vtkexportmodule.h"
#include "poiexportmodule.h"


// nonlocal barriers
#include "polylinenonlocalbarrier.h"
#include "symmetrybarrier.h"

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
#endif //__TM_MODULE


#ifdef __FM_MODULE
// Emodels 
#include "cbs.h"
// Elements
#include "tr1_2d_cbs.h"
// materials
#include "newtonianfluid.h"
// boundary conditions
#include "tractionpressurebc.h"

#include "supg.h"
#include "tr1_2d_supg.h"
#include "tr1_2d_supg2.h"
#include "tr1_2d_supg_axi.h"
#include "tr1_2d_supg2_axi.h"
#include "py1_3d_supg.h"
//#include "tr1_2d_supg99.h"
#include "twofluidmaterial.h"
#include "binghamfluid2.h"

#endif // __FM_Module

// GENERAL 
#include "masterdof.h"
#ifdef __PARALLEL_MODE
#include "loadbalancer.h"
#include "parmetisloadbalancer.h"
#endif

Element* CreateUsrDefElementOfType (char* aClass, int number, Domain* domain) 
{
   Element* newElement = NULL;
#ifdef __SM_MODULE
   if (! strncasecmp(aClass,"planestress2d",12))
     newElement = new PlaneStress2d(number,domain) ;
   if (! strncasecmp(aClass,"quad1planestrain",16))
     newElement = new Quad1PlaneStrain(number,domain) ;
   else if (! strncasecmp(aClass,"trplanestress2d",12))
     newElement = new TrPlaneStress2d(number,domain) ;
   else if (! strncasecmp(aClass,"trplanestrrot",12))
     newElement = new TrPlaneStrRot(number,domain) ;
   else if (! strncasecmp(aClass,"qplanestress2d",12))
     newElement = new QPlaneStress2d(number,domain) ; 
   else if (! strncasecmp(aClass,"qtrplstr",8))
     newElement = new QTrPlaneStress2d(number,domain) ;
   else if (! strncasecmp(aClass,"axisymm3d",9))
     newElement = new Axisymm3d (number,domain) ;
   else if (! strncasecmp(aClass,"q4axisymm",9))
     newElement = new Q4Axisymm (number,domain) ;
   else if (! strncasecmp(aClass,"l4axisymm",9))
     newElement = new L4Axisymm (number,domain) ;
   else if (! strncasecmp(aClass,"lspacebb",8))
     newElement = new LSpaceBB(number,domain) ; 
   else if (! strncasecmp(aClass,"lspace",6))
     newElement = new LSpace(number,domain) ; 
   else if (! strncasecmp(aClass,"qspace",6))
     newElement = new QSpace(number,domain) ; 
   else if (! strncasecmp(aClass,"cctplate",8))
     newElement = new CCTPlate(number,domain) ; 
//   else if (! strncasecmp(aClass,"ltrspaceec",10))
//     newElement = new LTRSpaceWithEmbeddedCrack (number,domain) ;
   else if (! strncasecmp(aClass,"ltrspace",8))
     newElement = new LTRSpace(number,domain) ; 
   else if (! strncasecmp(aClass,"truss2d",7))
     newElement = new Truss2d(number,domain) ; 
   else if (! strncasecmp(aClass,"rershell",8))
     newElement = new RerShell(number,domain) ;
   else if (! strncasecmp(aClass,"beam2d",12))
     newElement = new Beam2d (number,domain) ;
   else if (! strncasecmp(aClass,"beam3d",12))
     newElement = new Beam3d (number,domain) ;
  else if (! strncasecmp(aClass,"libeam2dNL",10))
     newElement = new LIBeam2dNL (number,domain) ;
   else if (! strncasecmp(aClass,"libeam2d",8))
     newElement = new LIBeam2d(number,domain) ;
  else if (! strncasecmp(aClass,"libeam3dnl2",11))
     newElement = new LIBeam3dNL2 (number,domain) ;
  else if (! strncasecmp(aClass,"libeam3dNL",10))
     newElement = new LIBeam3dNL (number,domain) ;
  else if (! strncasecmp(aClass,"truss3d",7))
     newElement = new Truss3d (number,domain) ;
  else if (! strncasecmp(aClass,"trplanestrain",13))
     newElement = new TrPlaneStrain (number,domain) ;
   else if (! strncasecmp(aClass,"libeam3d2",9))
     newElement = new LIBeam3d2(number,domain) ;
   else if (! strncasecmp(aClass,"libeam3d",8))
     newElement = new LIBeam3d (number,domain) ;
   else if (! strncasecmp(aClass,"truss1d",7))
     newElement = new Truss1d (number,domain) ;
   else if (! strncasecmp(aClass,"interface2dquad",15))
     newElement = new InterfaceElem2dQuad (number,domain) ;
   else if (! strncasecmp(aClass,"interface1d",11))
     newElement = new InterfaceElem1d (number,domain) ;
#endif //__SM_MODULE
#ifdef __TM_MODULE
  if (! strncasecmp(aClass,"quad1ht",7))
   newElement = new Quad1_ht(number,domain) ;
  else if (! strncasecmp(aClass,"tr1ht",5))
   newElement = new Tr1_ht(number,domain) ;
  else if (! strncasecmp(aClass,"quadaxisym1ht",13))
   newElement = new QuadAxisym1_ht(number,domain) ;
  else if (! strncasecmp(aClass,"traxisym1ht",11))
   newElement = new TrAxisym1_ht(number,domain) ;
  else if (! strncasecmp(aClass,"quad1hmt",8))
   newElement = new Quad1_ht(number,domain,Quad1_ht::HeatMass1TransferEM) ;
  else if (! strncasecmp(aClass,"quadaxisym1hmt",14))
   newElement = new QuadAxisym1_ht(number,domain,Quad1_ht::HeatMass1TransferEM) ;
  else if (! strncasecmp(aClass,"brick1ht",8))
   newElement = new Brick1_ht(number,domain) ;
  else if (! strncasecmp(aClass,"brick1hmt",9))
   newElement = new Brick1_ht(number,domain,Brick1_ht::HeatMass1TransferEM) ;
  else if (! strncasecmp(aClass,"tetrah1ht",9))
   newElement = new Tetrah1_ht(number,domain) ;
  else if (! strncasecmp(aClass,"tetrah1hmt",10))
   newElement = new Tetrah1_ht(number,domain,Tetrah1_ht::HeatMass1TransferEM) ;
#endif //__TM_MODULE
#ifdef __FM_MODULE
  else if (! strncasecmp(aClass,"tr1cbs",6))
    newElement = new TR1_2D_CBS (number,domain) ;
  else if (! strncasecmp(aClass,"tr1supgaxi",10))
    newElement = new TR1_2D_SUPG_AXI (number,domain) ;
  else if (! strncasecmp(aClass,"tr1supg2axi",11))
    newElement = new TR1_2D_SUPG2_AXI (number,domain) ;
  else if (! strncasecmp(aClass,"tr1supg2",8))
    newElement = new TR1_2D_SUPG2 (number,domain) ;
  /*
  else if (! strncasecmp(aClass,"tr1supg99",9))
    newElement = new TR1_2D_SUPG99 (number,domain) ;
  */
  else if (! strncasecmp(aClass,"tr1supg",7))
    newElement = new TR1_2D_SUPG (number,domain) ;
  else if (! strncasecmp(aClass,"py1supg",7))
    newElement = new PY1_3D_SUPG (number,domain) ;

#endif //__FM_MODULE
   return newElement ;
}

DofManager* CreateUsrDefDofManagerOfType (char* aClass, int number, Domain* domain) 
{
   DofManager* newDofManager=NULL ;
/*
   if (! strncasecmp(aClass,"planestress2d",12))
  newDofManager = new PlaneStress2d(number,domain) ;
*/
  return newDofManager;
 }

CrossSection* CreateUsrDefCrossSectionOfType(char* aClass, int number, Domain* domain)
{
 CrossSection* newCS = NULL;
#ifdef __SM_MODULE
 if (! strncasecmp(aClass,"layeredcs",14))
  newCS =   new LayeredCrossSection (number,domain) ;
 else if (! strncasecmp(aClass,"fiberedcs",9))
  newCS =   new FiberedCrossSection (number,domain) ;
#endif //__SM_MODULE
 return newCS;
}

EngngModel* CreateUsrDefEngngModelOfType(char* aClass, int number, EngngModel* master)
{

 EngngModel* newEModel = NULL;

#ifdef __SM_MODULE
 if (! strncasecmp(aClass,"linearstatic",12))
  newEModel = new LinearStatic(number,master) ;
 else if (! strncasecmp(aClass,"stationaryflow",14))
  newEModel = new StationaryFlow(number,master) ; 
 else if (! strncasecmp(aClass,"eigenvaluedynamic",14))
  newEModel = new EigenValueDynamic(number,master) ; 
 else if (! strncasecmp(aClass,"nonlinearstatic",13))
  newEModel = new NonLinearStatic(number,master) ; 
 else if (! strncasecmp(aClass,"nldeidynamic",8))
  newEModel = new NlDEIDynamic(number,master) ; 
//#ifdef __PARALLEL_MODE
 else if (! strncasecmp(aClass,"pnldeidynamic",9))
  newEModel = new PNlDEIDynamic(number,master) ; 
//#endif
 else if (! strncasecmp(aClass,"deidynamic",10))
  newEModel = new DEIDynamic(number,master) ; 
 else if (! strncasecmp(aClass,"diidynamic",10))
  newEModel = new DIIDynamic(number,master) ; 
 else if (! strncasecmp(aClass,"incrlinearstatic",16))
  newEModel = new IncrementalLinearStatic(number,master) ; 
 else if (! strncasecmp(aClass,"linearstability",15))
  newEModel = new LinearStability(number,master) ; 
 else if (! strncasecmp(aClass,"adaptnlinearstatic",18))
  newEModel = new AdaptiveNonLinearStatic(number,master) ; 
 else if (! strncasecmp(aClass,"adaptlinearstatic",17))
  newEModel = new AdaptiveLinearStatic(number,master) ; 
#ifdef __PARALLEL_MODE
 else if (! strncasecmp(aClass,"plinearstatic",13))
  newEModel = new PLinearStatic(number,master) ; 
#endif
#endif //__SM_MODULE

#ifdef __TM_MODULE
 if (! strncasecmp(aClass,"stationaryproblem",13))
  newEModel = new StationaryTransportProblem(number,master) ; 
 else if (! strncasecmp(aClass,"nonstationaryproblem",16))
  newEModel = new NonStationaryTransportProblem(number,master) ; 
 else if (! strncasecmp(aClass,"nltransienttransportproblem",23))
  newEModel = new NLTransientTransportProblem(number,master) ; 
 else if (! strncasecmp(aClass,"staggeredproblem",16))
  newEModel = new StaggeredProblem(number,master) ; 
#endif //__TM_MODULE

#ifdef __FM_MODULE
 if (! strncasecmp(aClass,"cbs",3))
   newEModel = new CBS (number,master) ; 
 else if (! strncasecmp(aClass,"supg",4))
   newEModel = new SUPG(number,master) ; 
#endif //__FM_MODULE

 if (newEModel == NULL) {
  printf ("%s : unknown EngngModel type \n",aClass) ;
  exit(0) ;
 }
 return newEModel;
}

GeneralBoundaryCondition* CreateUsrDefBoundaryConditionOfType(char* aClass, int number, Domain* domain)
{
   GeneralBoundaryCondition* newBc = NULL;
  
#ifdef __SM_MODULE
  if (!strncasecmp(aClass,"structtemperatureload",18))
  newBc = new StructuralTemperatureLoad(number,domain) ;
  else if (!strncasecmp(aClass,"linearedgeload",14))
  newBc = new LinearEdgeLoad(number,domain) ;
  else if (!strncasecmp(aClass,"constantedgeload",16))
  newBc = new ConstantEdgeLoad(number,domain) ;
  else if (!strncasecmp(aClass,"constantsurfaceload",19))
  newBc = new ConstantSurfaceLoad(number,domain) ;
  else if (!strncasecmp(aClass,"usrdeftempfield",15))
  newBc = new UserDefinedTemperatureField (number,domain) ;
  else if (!strncasecmp(aClass,"tf1",3))
  newBc = new TF1 (number,domain) ;
  else if (!strncasecmp(aClass,"pointload",9))
  newBc = new PointLoad (number,domain) ;
#endif //__SM_MODULE
#ifdef __FM_MODULE
  if (!strncasecmp(aClass,"prescribedtractionpressurebc",28))
  newBc = new TractionPressureBC(number,domain) ;
#endif  
  return newBc;
 }

LoadTimeFunction* CreateUsrDefLoadTimeFunctionOfType(char* aClass, int number, Domain* domain)
{
   LoadTimeFunction* newLTF = NULL;

#ifdef __SM_MODULE
  if (! strncasecmp(aClass,"peakfunction",5))
      newLTF = new PeakFunction(number,domain) ;
   else if (! strncasecmp(aClass,"piecewiselinfunction",5))
      newLTF = new PiecewiseLinFunction(number,domain) ;
   else if (! strncasecmp(aClass,"periodicpiecewiselinfunction",5))
      newLTF = new PeriodicPiecewiseLinFunction(number,domain) ;
   else if (! strncasecmp(aClass,"heavisideltf",12))
      newLTF = new HeavisideLTF(number,domain) ;
   else if (! strncasecmp(aClass,"usrdefltf",9))
      newLTF = new UserDefinedLoadTimeFunction(number,domain) ;
#endif //__SM_MODULE

  return newLTF;
}

Material* CreateUsrDefMaterialOfType(char* aClass, int number, Domain* domain)
{

  Material* newMaterial = NULL;
  
#ifdef __SM_MODULE
 if (! strncasecmp(aClass,"orthole",7))
  newMaterial = new OrthotropicLinearElasticMaterial(number,domain);
  else if (! strncasecmp(aClass,"steel1",6))
  newMaterial = new Steel1 (number,domain) ;
  else if (! strncasecmp(aClass,"concrete2",9))
  newMaterial = new Concrete2 (number,domain) ;
  else if (! strncasecmp(aClass,"concrete3",9))
  newMaterial = new Concrete3 (number,domain) ;
  else if (! strncasecmp(aClass,"cebfip78",8))
  newMaterial = new CebFip78Material (number,domain) ;
  else if (! strncasecmp(aClass,"doublepowerlaw",14))
  newMaterial = new DoublePowerLawMaterial (number,domain) ;
  else if (! strncasecmp(aClass,"b3mat",5))
  newMaterial = new B3Material (number,domain) ;
  else if (! strncasecmp(aClass,"j2mat",5))
  newMaterial = new J2plasticMaterial (number, domain);
  else if (! strncasecmp(aClass,"rcsdnl",6))
  newMaterial = new RCSDNLMaterial (number, domain);
  else if (! strncasecmp(aClass,"rcsde",5))
  newMaterial = new RCSDEMaterial (number, domain);
  else if (! strncasecmp(aClass,"rcsd",4))
  newMaterial = new RCSDMaterial (number, domain);
 else if (! strncasecmp(aClass,"microplane_m4",13))
    newMaterial = new M4Material (number,domain);
 else if (! strncasecmp(aClass,"idm1",4))
    newMaterial = new IsotropicDamageMaterial1 (number,domain);
 else if (! strncasecmp(aClass,"idmnl1",6))
    newMaterial = new IDNLMaterial (number,domain);
 else if (! strncasecmp(aClass,"mazarsmodelnl",13))
    newMaterial = new MazarsNLMaterial (number,domain);
 else if (! strncasecmp(aClass,"mazarsmodel",11))
    newMaterial = new MazarsMaterial (number,domain);
 else if ( ! strncmp(aClass, "druckerprager", 13 ) )
   newMaterial = new DruckerPragerPlasticitySM(number, domain) ;
 else if (! strncasecmp(aClass,"j2mmat",6))
   newMaterial = new J2MPlasticMaterial (number, domain);
 else if (! strncasecmp(aClass,"rankine",7))
   newMaterial = new RankinePlasticMaterial (number, domain);
 else if (! strncasecmp(aClass,"masonry02",9))
   newMaterial = new Masonry02 (number, domain);
 else if (! strncasecmp(aClass,"isointrfdm01",12))
   newMaterial = new IsoInterfaceDamageMaterial (number, domain);
 else if (! strncasecmp(aClass,"j22mat",6))
   newMaterial = new J2Mat (number, domain);
 else if (! strncasecmp(aClass,"cebfipslip90",12))
   newMaterial = new CebFipSlip90Material (number, domain);
//modified_zh 7.6.2004
 else if ( ! strncmp(aClass, "hellmat", 7 ) )
   newMaterial = new HellmichMaterial(number, domain) ;

 else if (! strncasecmp(aClass,"mdm",3))
    newMaterial = new MDM (number,domain);
#endif //__SM_MODULE

#ifdef __TM_MODULE
  if (! strncasecmp(aClass,"isoheat",7))
   newMaterial = new IsotropicHeatTransferMaterial(number,domain) ;
  else if (! strncasecmp(aClass,"hemotk",6))
   newMaterial = new HeMoTKMaterial (number,domain);
  else if ( ! strncmp(aClass, "hisoheat", 8 ) )
    newMaterial = new HydratingIsoHeatMaterial(number, domain) ;
  else if ( ! strncmp(aClass, "hhemotk", 7 ) )
    newMaterial = new HydratingHeMoMaterial(number, domain) ;
 
#endif //__TM_MODULE

#ifdef __FM_MODULE
  if (! strncasecmp(aClass,"newtonianfluid",14))
    newMaterial = new NewtonianFluidMaterial(number,domain) ;
  else if (! strncasecmp(aClass,"twofluidmat",11))
    newMaterial = new TwoFluidMaterial(number,domain) ;
  else if ( ! strncasecmp(aClass,"binghamfluid2",13))
    newMaterial = new BinghamFluidMaterial2(number,domain) ;
  else if ( ! strncasecmp(aClass,"binghamfluid",12))
    newMaterial = new BinghamFluidMaterial2(number,domain) ;
#endif // __FM_MODULE

 return newMaterial;
}
 
SparseMtrx* CreateUsrDefSparseMtrx(SparseMtrxType type)
{
 SparseMtrx* answer = NULL;

 if (type == SMT_Skyline) answer = new Skyline ();
 else if (type == SMT_SkylineU) answer = new SkylineUnsym();
#ifdef __IML_MODULE
 else if (type == SMT_CompCol)  answer = new CompCol ();
 else if (type == SMT_DynCompCol) answer = new DynCompCol ();
 else if (type == SMT_SymCompCol) answer = new SymCompCol ();
 else if (type == SMT_DynCompRow) answer = new DynCompRow ();
#endif
#ifdef __SPOOLES_MODULE
 else if (type == SMT_SpoolesMtrx) answer = new SpoolesSparseMtrx ();
#endif
#ifdef __PETSC_MODULE
 else if (type == SMT_PetscMtrx) answer = new PetscSparseMtrx ();
#endif
#ifdef __DSS_MODULE
 else if (type == SMT_DSS_sym_LDL)  answer = new DSSMatrix(DSSMatrix::sym_LDL);
 else if (type == SMT_DSS_sym_LL)   answer = new DSSMatrix(DSSMatrix::sym_LL);
 else if (type == SMT_DSS_unsym_LU) answer = new DSSMatrix(DSSMatrix::unsym_LU);
#endif
 else {
   fprintf(stderr, "CreateUsrDefSparseMtrx: Unknown mtrx type\n");
   exit(1);
 }

 return answer;
}

SparseLinearSystemNM* CreateUsrDefSparseLinSolver(LinSystSolverType st, int i, Domain* d, EngngModel* m)
{
  
  SparseLinearSystemNM* nm = NULL;
  if (st == ST_Direct) {
    nm = (SparseLinearSystemNM*) new LDLTFactorization (i,d,m);
    return nm;
  } else if (st == ST_IML) {
    nm = (SparseLinearSystemNM*) new IMLSolver (i,d,m);
    return nm;
  } else if (st == ST_Spooles) {
    nm = (SparseLinearSystemNM*) new SpoolesSolver (i,d,m);
    return nm;
  } else if (st == ST_Petsc) {
    nm = (SparseLinearSystemNM*) new PetscSolver (i,d,m);
    return nm;
#ifdef __PARALLEL_MODE
  } else if (st == ST_Feti) {
    nm = (SparseLinearSystemNM*) new FETISolver (i,d,m);
    return nm;
#endif
#ifdef __DSS_MODULE
  } else if (st == ST_DSS) {
    nm = (SparseLinearSystemNM*) new DSSSolver (i,d,m);
    return nm;
#endif
  } else {
    fprintf(stderr, "CreateUsrDefSparseLinSolver: Unknown solver type\n");
    exit(1);
  }
  return nm;
}

SparseGeneralEigenValueSystemNM* CreateUsrDefGeneralizedEigenValueSolver(GenEigvalSolverType st, int i, Domain* d, EngngModel* m)
{
  
  SparseGeneralEigenValueSystemNM* nm = NULL;
  if (st == GES_SubspaceIt) {
    nm = (SparseGeneralEigenValueSystemNM*) new SubspaceIteration (i,d,m);
    return nm;
  } else if (st == GES_InverseIt) {
    nm = (SparseGeneralEigenValueSystemNM*) new InverseIteration (i,d,m);
    return nm;
  } else {
    fprintf(stderr, "CreateUsrDefGeneralizedEigenValueSolver: Unknown solver type\n");
    exit(1);
  }
  return nm;
}



ErrorEstimator*
CreateUsrDefErrorEstimator(ErrorEstimatorType type, int number, Domain* d)
{
 ErrorEstimator* answer = NULL;

#ifdef __SM_MODULE
 if (type == EET_SEI) answer = new ScalarErrorIndicator (number, d);
 else if (type == EET_ZZEE) answer = new ZZErrorEstimator(number, d);
 else if (type == EET_CZZSI) answer = new CombinedZZSIErrorEstimator (number, d);
 else if (type == EET_HEE) answer = new HuertaErrorEstimator (number, d);
#endif //__SM_MODULE

 if (answer == NULL){
  fprintf(stderr, "CreateUsrDefErrorEstimator: Unknown error estimator type\n");
  exit(1);
 }

 return answer;
}

ExportModule* CreateUsrDefExportModuleOfType (char* aClass, EngngModel* emodel)
{
 ExportModule *answer = NULL;

#ifdef __SM_MODULE
 if (! strncasecmp(aClass,"vtk",3))
   answer = new VTKExportModule (emodel);
 else if (! strncasecmp(aClass,"poi",3))
   answer = new POIExportModule (emodel);
#endif //__SM_MODULE

 return answer;
}


NonlocalBarrier* CreateUsrDefNonlocalBarrierOfType (char* aClass, int num, Domain*d)
{
 NonlocalBarrier *answer = NULL;

#ifdef __SM_MODULE
 if (! strncasecmp(aClass,"polylinebarrier",15))
    answer = new PolylineNonlocalBarrier (num, d);
 if (! strncasecmp(aClass,"symmetrybarrier",15))
    answer = new SymmetryBarrier (num, d);
#endif //__SM_MODULE

 return answer;
}

IntegrationRule* CreateUsrDefIRuleOfType (classType type, int number, Element* e)
{
  IntegrationRule* answer = NULL;
  if (type == GaussIntegrationRuleClass) answer = new GaussIntegrationRule (number, e);

  if (answer == NULL)  OOFEM_ERROR2 ("CreateUsrDefIRuleOfType: Unknown integration rule type [%d]", type);
  return answer;
}


Element* CreateUsrDefElementOfType (classType type, int number, Domain* domain)
{
  Element* answer = NULL;

  if (type == PlaneStress2dClass) answer = new PlaneStress2d(number,domain) ;
  else if (type == TrPlaneStress2dClass) answer = new TrPlaneStress2d (number, domain);
  else if (type == LTRSpaceClass) answer = new LTRSpace (number, domain);

  if (answer == NULL)  OOFEM_ERROR2 ("CreateUsrDefElementOfType: Unknown element type [%d]", type);
  return answer;
}

DofManager* CreateUsrDefDofManagerOfType (classType type, int number, Domain* domain)
{
  DofManager* answer = NULL;
  if (type == NodeClass) answer = new Node (number,domain) ;

  if (answer == NULL)  OOFEM_ERROR2 ("CreateUsrDefDofManagerOfType: Unknown dofman type [%d]", type);
  return answer;
}

Dof* CreateUsrDefDofOfType (classType type, int number, DofManager* dman)
{
  Dof* answer = NULL;
  if  (type == MasterDofClass) answer = new MasterDof (number, dman);

  if (answer == NULL)  OOFEM_ERROR2 ("CreateUsrDefDofOfType: Unknown dof type [%d]", type);
  return answer;
}

MaterialMappingAlgorithm* CreateUsrDefMaterialMappingAlgorithm (MaterialMappingAlgorithmType type)
{
  MaterialMappingAlgorithm* answer=NULL;
  if (type == MMA_ClosestPoint) answer = new MMAClosestIPTransfer ();
  else if (type == MMA_LeastSquareProjection) answer = new MMALeastSquareProjection();
  else if (type == MMA_ShapeFunctionProjection) answer = new MMAShapeFunctProjection();

  if (answer == NULL) OOFEM_ERROR2 ("CreateUsrDefMaterialMappingAlgorithm: Unknown mma type [%d]", type);
  return answer;

}


#ifdef __PARALLEL_MODE
LoadBalancerMonitor* CreateUsrDefLoadBalancerMonitorOfType (classType type, EngngModel* e)
{
  LoadBalancerMonitor *answer = NULL;
  if (type == WallClockLoadBalancerMonitorClass) answer = new WallClockLoadBalancerMonitor (e);
  
  if (answer == NULL) OOFEM_ERROR2 ("CreateUsrDefLoadBalancerMonitorOfType: Unknown type [%d]", type);
  return answer;
}

LoadBalancer* CreateUsrDefLoadBalancerOfType (classType type, Domain* d)
{
  LoadBalancer *answer = NULL;
  if (type == ParmetisLoadBalancerClass) answer = new ParmetisLoadBalancer (d);
  
  if (answer == NULL) OOFEM_ERROR2 ("CreateUsrDefLoadBalancerOfType: Unknown type [%d]", type);
  return answer;
}
#endif
