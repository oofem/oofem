/* $Header: /home/cvs/bp/oofem/sm/src/maxwellChM.h,v 1.6 2003/04/06 14:08:31 bp Exp $ */
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

//   *********************************
//   *** CLASS Maxwelll Chain Model ***
//   *********************************

#ifndef maxwellchm_h 
#define maxwellchm_h 

#include "femcmpnn.h"
#include "structuralmaterial.h"
#include "linearelasticmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "matconst.h"
#include "structuralelement.h"
#include "structuralms.h"

#define MNC_NPOINTS 30
#define TIME_DIFF   1.e-10


class MaxwellChainMaterialStatus : public StructuralMaterialStatus
{
/*
 This class implements associated Material Status to MaxwellChainMaterial.
 It is atribute of matStatusDictionary at every GaussPoint, for which this material 
 is active. Isotropic linear elastic material is been assumed.
 DESCRIPTION:
   Idea used there is that we have variables
   describing:
     1) state at previous equilibrium state (variables without temp)
   2) state during searching new equilibrium (variables with temp)
   when we start search new state from previous equilibrium one we copy 
   non-tem variables into temp ones. And after we reach new equilibrium 
   (now decribed by temp variables) we copy tem-var into non-tepm ones
   (see function updateYourself).
   
   variables description:


  TASK:

*/

protected:
 FloatArray** hiddenStreses;
 int nMaxwellUnits;

  /// total values of shrinkage strain (nneded when only incremen tal shrinkage formulation exist)
  FloatArray shrinkageStrain;
public:
 MaxwellChainMaterialStatus (int n, Domain *d, GaussPoint* g, int nUnits);
 ~MaxwellChainMaterialStatus ();
 void printOutputAt (FILE *file, TimeStep* tStep) ;

 FloatArray* giveHiddenStressVector (int i) {return hiddenStreses[i-1];}
 FloatArray* letHiddenStressVectorBe (int i, FloatArray*);

  FloatArray* giveShrinkageStrainVector() {return &shrinkageStrain;}
  void        setShrinkageStrainVector(const FloatArray& src) {shrinkageStrain=src;}

 virtual void initTempStatus () ;
 virtual void updateYourself(TimeStep*);    // update after new equilibrium state reached

 // saves current context(state) into stream
 contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);
 contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);
  
 // definition
 const char* giveClassName () const { return "MaxwellChainMaterialStatus" ;}
   classType             giveClassID () const 
   { return MaxwellChainMaterialStatusClass; }
};




class MaxwellChainMaterial : public StructuralMaterial
{
/*
   This class implements a rheologic Maxwelll chain model in a finite
 element problem. 
 
 DESCRIPTION
 TASK
*/

protected:
 
 int nChainUnits;
 double relMatAge;
 double nu;
 double Eval;
 double EmjuValTime;

  double endOfTimeOfInterest; // local one or taken from e-model
 // linearElasticMaterial should have E = 1;
 LinearElasticMaterial *linearElasticMaterial;
 FloatArray EmjuVal;
 FloatArray relaxationTimes;
 FloatArray discreteTimeScale;

  /**time coefficient to transform solotion time 
  (for example soltime/timeFactor = time in days) required by this model*/
  double timeFactor;  



   public :
 MaxwellChainMaterial (int n,Domain* d);
 ~MaxwellChainMaterial ();
 
 // standart matrial stiffness matrices
   virtual void  giveCharacteristicMatrix (FloatMatrix& answer,
                      MatResponseForm form,
                      MatResponseMode mode,
                      GaussPoint* gp,
                      TimeStep* atTime);
 
 // virtual FloatArray*  GiveRealStressVector3d ( GaussPoint*, FloatArray*) ;
 virtual void giveRealStressVector (FloatArray& answer,  MatResponseForm, GaussPoint*, 
                   const FloatArray&, TimeStep*) ;
 
 // returns a FloatArray(3) of coefficients of thermal dillatation in direction
 // of each (local) axisgiven by principal axis of material
 // 
 virtual void giveThermalDilatationVector (FloatArray& answer, GaussPoint*, TimeStep*) 
  {answer.resize (0);}
 
 // uptates MatStatus to newly reched (equilibrium) state
 virtual void updateYourself (GaussPoint* gp, TimeStep*);
 
 // identification and auxiliary functions
 virtual int hasNonLinearBehaviour ()   { return 0 ;}
  virtual int hasMaterialModeCapability (MaterialMode mode);
 const char*    giveClassName () const      { return "MaxwelllChainMaterial" ;}
 classType giveClassID ()         const      {return MaxwelllChainMaterialClass;}
  IRResultType initializeFrom (InputRecord* ir);
 void     printYourself () ;

 // store & restore context functions
 contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);
 contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);


 virtual void  give3dMaterialStiffnessMatrix (FloatMatrix& answer,
                        MatResponseForm, MatResponseMode, 
                        GaussPoint* gp,
                        TimeStep *atTime);

/**
 Computes strain vector in given integration point, generated by internal processes in
 material, which are independent on loading in particular integration point.
 Default implementation takes only into account temperature induced strains.
 Overloaded to take into account shrinkage and eigen strains.
 @param answer returned strain vector
 @param gp integration point
 @param atTime time step (most models are able to respond only when atTime is current time step)
 @param determines response mode (Total or incremental)
 */
  virtual void computeStressIndependentStrainVector (FloatArray& answer,
                           GaussPoint *gp, TimeStep *stepN, ValueModeType mode);
/**
 Computes strain vector in given integration point, generated by internal processes in
 material, which are independent on loading in particular integration point.
 @param answer returned strain vector
 @param form material response form
 @param gp integration point
 @param atTime time step (most models are able to respond only when atTime is current time step)
 @param determines response mode (Total or incremental)
 */
 virtual void  giveShrinkageStrainVector (FloatArray& answer,
                      MatResponseForm form, 
                      GaussPoint* gp, 
                      TimeStep* atTime,
                      ValueModeType mode) 
  {answer.resize(0);}
 
 // Note: must take LoadResponseMode into account
/**
 Computes strain vector for given integration point, which is induced by stress history 
 in given integration point (typycally creep strain)
 @param answer computed strains
 @param form material response form
 @param gp integration point
 @param atTime time step (most models are able to respond only when atTime is current time step)
 @param mode determines response mode
 */
  virtual void  giveEigenStrainVector (FloatArray& answer, MatResponseForm form, 
                    GaussPoint* gp, TimeStep* atTime, ValueModeType mode);

#ifdef __OOFEG
#endif

 virtual MaterialStatus* CreateStatus (GaussPoint* gp) const;



protected:
 
 /** if only incremental shrinkage strain formulation is provided, then total shrinkage strain must be tracked
  in status in order to be able to compute total value. */
 virtual int  hasIncrementalShrinkageFromulation() {return 0;}

 void         generateLogTimeScale (FloatArray& answer, double from, double to, int nsteps,
                   int fromIncluded = 0);
 const FloatArray&  giveDiscreteTimes ();
 virtual double  computeCreepFunction (GaussPoint* gp, double ofAge, double atTime) = 0;
 void         computeCharCoeficients (FloatArray& answer, GaussPoint* gp, double);
 void         computeDiscreteRelaxationFunction (FloatArray& answer, GaussPoint* gp, 
                         const FloatArray& atTimes,
                         double t0, double tr);
 void         giveGeneralizationBMatrix (FloatMatrix& answer, MatResponseForm, 
                     GaussPoint* gp, TimeStep* tStep);
 void         giveGeneralizationBInvMatrix (FloatMatrix& answer,
                       MatResponseForm, GaussPoint* gp, TimeStep* tStep);
 double       giveEModulus (GaussPoint* gp, TimeStep* atTime);
 void         updateEmjuModuluses (GaussPoint* gp, double atTime);
 double       giveEmjuModulus (int iChain);
 void         computeRelaxationTimes();
 double       giveRelaxationTime (int) ;
 virtual double giveRelaxationTimeExponent (int i) {return 1.0;}
 LinearElasticMaterial* giveLinearElasticMaterial ();

  double giveEndOfTimeOfInterest();

 virtual void  givePlaneStressStiffMtrx (FloatMatrix& answer, 
                     MatResponseForm,MatResponseMode,
                     GaussPoint* gp, 
                     TimeStep* atTime);
 virtual void  givePlaneStrainStiffMtrx (FloatMatrix& answer,
                     MatResponseForm,MatResponseMode,
                     GaussPoint* gp,
                     TimeStep* atTime);
 virtual void  give1dStressStiffMtrx (FloatMatrix& answer,
                    MatResponseForm,MatResponseMode,
                    GaussPoint* gp,
                    TimeStep* atTime);
 virtual void  give2dBeamLayerStiffMtrx (FloatMatrix& answer,
                     MatResponseForm,MatResponseMode,
                     GaussPoint* gp,
                     TimeStep* atTime);
 virtual void give2dPlateLayerStiffMtrx(FloatMatrix& answer,
                     MatResponseForm,MatResponseMode, 
                     GaussPoint* gp,
                     TimeStep* atTime);
 virtual void give3dShellLayerStiffMtrx(FloatMatrix& answer,
                     MatResponseForm,MatResponseMode,
                     GaussPoint* gp,
                     TimeStep* atTime);
/**
 Computes strain vector in given integration point, generated by internal processes in
 material, which are independent on loading in particular integration point.
 Takes only into account temperature and shrinkage induced strains.
 @param answer returned strain vector
 @param gp integration point
 @param atTime time step (most models are able to respond only when atTime is current time step)
 @param determines response mode (Total or incremental)
 */
  void computeTrueStressIndependentStrainVector (FloatArray& answer,GaussPoint *gp, 
                         TimeStep *stepN, ValueModeType mode);
                           
} ;


#endif // maxwellchm_h 
