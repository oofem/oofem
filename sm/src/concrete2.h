/* $Header: /home/cvs/bp/oofem/sm/src/concrete2.h,v 1.5 2003/04/06 14:08:30 bp Exp $ */
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

//   ****************************************************************************
//       CLASS CONCRETE2 - NonLinear elasto-plastic material model with hardening 
//                         Plane stress or uniaxial stress + transverse shear in
//                         concrete layers with transverse stirrups.
//   ****************************************************************************

#ifndef concrete2_h

#include "femcmpnn.h"
#include "dictionr.h"
#include "material.h"
#include "deformationtheorymaterial.h"
#include "isolinearelasticmaterial.h"
#include "flotarry.h"
#include "cltypes.h"
#include "structuralms.h"

#define c2_SCCC  300
#define c2_SCCT  301
#define c2_EPP   302
#define c2_EPU   303
#define c2_EOPP  304
#define c2_EOPU  305
#define c2_SHEARTOL 306
#define c2_E     307
#define c2_n     308
#define stirr_E  309
#define stirr_Ft 310
#define stirr_A  311
#define stirr_TOL 312
#define stirr_EREF 313
#define stirr_LAMBDA 314
#define c2_IS_PLASTIC_FLOW 315
#define c2_IFAD  316




class Concrete2MaterialStatus : public StructuralMaterialStatus
{
/*
 This class implements associated Material Status to Concrete2Material.
 It is atribute of matStatusDictionary at every GaussPoint, for which this material 
 is active.
 DESCRIPTION:
    This class contains state variables of the materil (layer)
  
 VARIABLE DESCRIPTION:
    SCCM      current pressure strenght
  EPM       max. eff. plastic strain
  SCTM      current tension strenght
  E0PM      max. vol. plastic strain
  SRF       current stress in stirrups
  SEZ       current strain in transverse (z) direction.
  
   TASK:
    
    returning and setting variables
  printing
  saving & restoring context
*/
protected:

  FloatArray  plasticStrainVector; // full form
  FloatArray  plasticStrainIncrementVector;

  double tempSCCM, tempEPM, tempSCTM, tempE0PM, tempSRF, tempSEZ;
  double SCCM, EPM, SCTM, E0PM, SRF, SEZ;
 


public:
 Concrete2MaterialStatus (int n, Domain*d, GaussPoint* g);
 ~Concrete2MaterialStatus () ;
 void   printOutputAt (FILE *file, TimeStep* tStep) 
        {StructuralMaterialStatus :: printOutputAt(file, tStep);}

 void givePlasticStrainVector(FloatArray& answer) const {answer = plasticStrainVector;}
 void givePlasticStrainIncrementVector(FloatArray& answer) const 
  {answer = plasticStrainIncrementVector;}
  void         letPlasticStrainVectorBe (FloatArray& v)
  {plasticStrainVector = v;}
  void         letPlasticStrainIncrementVectorBe (FloatArray& v) 
  {plasticStrainIncrementVector = v;}

 double& giveTempCurrentPressureStrength () {return tempSCCM;}
 double& giveTempMaxEffPlasticStrain ()     {return tempEPM ;}
 double& giveTempCurrentTensionStrength ()  {return tempSCTM;}
 double& giveTempCurrentStressInStirrups () {return tempSRF ;}
 double& giveTempCurrentStrainInZDir ()     {return tempSEZ ;}
 double& giveTempMaxVolPlasticStrain ()     {return tempE0PM;}

 // query for non-tem variables (usefull for postprocessing)
 double& giveCurrentPressureStrength () {return SCCM;}
 double& giveMaxEffPlasticStrain ()     {return EPM ;}
 double& giveCurrentTensionStrength ()  {return SCTM;}
 double& giveCurrentStressInStirrups () {return SRF ;}
 double& giveCurrentStrainInZDir ()     {return SEZ ;}
 double& giveMaxVolPlasticStrain ()     {return E0PM;}


 virtual void initTempStatus () ;
 virtual void updateYourself(TimeStep*) ;  // update after new equilibrium state reached

 // saves current context(state) into stream
 contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);
 contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);
  
 // definition
 const char* giveClassName () const { return "Concrete2MaterialStatus" ;}
   classType             giveClassID () const { return Concrete2MaterialStatusClass; }
 
}; 



/* Material constant description:
 SCCC  - pressure strength
 SCCT  - tension strength
 EPP   - treshold eff. plastic strain for softening in compress.
 EPU   - ultimate eff. pl. strain
 EOPP  - threshold vlumetric plastic strain for soft. in tension
 EOPU  - ultimate vol. pl. strain
 SHEARTOL -  threshold value of the relative shear deformation
             (psi**2/eef) at which shear is considered in layers. for
             lower r.s.d. the transverse shear remains elastic decoupled
             from bending. default value SHEARTOL = 0.01 
IS_PLASTIC_FLOW 
   indicates that plastic flow (not deformation theory) 
   is used in pressure.
IFAD 
  IFAD<=0 STATE VARIABLES WILL NOT BE UPDATED
       >0 UPDATE S.V.

STIRR_MAT - material id num for stirrups
*/


class Concrete2 : public DeformationTheoryMaterial
{
/*
   This class implements a Concrete2 material in a finite element problem. A material
   is an attribute of a domain. It is usually also attribute of many elements.
 DESCRIPTION
     Plane stress or uniaxial stress + transverse shear in
     concrete layers with transverse stirrups.

   The attribute 'propertyDictionary' contains all the properties of a mate-
   rial, like its Young modulus, its mass density or poisson ratio.
 TASK
   - Returning standard material marices (like 3dstress-strain, 2d plane , plate,
     3dbeam, 2d beam ..) according to current state determined by using data stored 
     in Gausspoint.
   - Returning a material property (method 'give'). Only for non-standard elements.
     In this case such element is fully responsible to compute material matrix, handle 
     material nonlinearity and so on.
*/
private:
  double SCCC, SCCT, EPP, EPU, EOPP, EOPU, SHEARTOL;
  double E, n;
  // stirrups
  double stirrE, stirrFt, stirrA, stirrTOL, stirrEREF, stirrLAMBDA;
  int IS_PLASTIC_FLOW, IFAD; 
  
  LinearElasticMaterial *linearElasticMaterial;
 
public:
 
 Concrete2  (int n,Domain* d);
 ~Concrete2 () ;
 
 virtual void  giveRealStressVector (FloatArray& answer, MatResponseForm, GaussPoint*, 
                   const FloatArray&, TimeStep* atTime) ;

 virtual void  give3dMaterialStiffnessMatrix (FloatMatrix& answer,
                        MatResponseForm, MatResponseMode,
                        GaussPoint* gp,
                        TimeStep* atTime);

 MaterialStatus* CreateStatus (GaussPoint* gp) const;

 protected:
   
   void  giveRealStresses3dShellLayer (FloatArray& answer, MatResponseForm form, GaussPoint* gp, 
                    const FloatArray& strain, TimeStep* atTime) ;
   void  dtp3 (GaussPoint* gp, FloatArray* e, FloatArray* s, FloatArray* ep, 
               double SCC, double SCT, int* ifplas);
   void  dtp2 (GaussPoint* gp, FloatArray* e, FloatArray* s, FloatArray* ep, 
               double SCC, double SCT, int* ifplas);
   void  stirr (double dez, double srf);
   void  strsoft (GaussPoint *gp, double epsult, FloatArray *ep, double &ep1, 
      double &ep2, double &ep3, double SCC, double SCT ,int& ifupd);

 // two functions used to initialize and updating temporary variables in
 // gp's status. These variables are used to controll process, when
 // we try to find equlibrium state.

 void updateStirrups (GaussPoint* gp, FloatArray* strainIncrement);
public:

 // non-standart
 double   give (int) ; 

   // identification and auxiliary functions
 
   int       hasNonLinearBehaviour ()     { return 1 ;}
   const char* giveClassName ()   const   { return "Concrete2" ;}
   classType giveClassID ()            const   {return Concrete2Class;}
   IRResultType initializeFrom (InputRecord* ir);
   //      void     printYourself () ;
   
   
} ;


#define concrete2_h
#endif






