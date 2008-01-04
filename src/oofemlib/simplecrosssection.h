/* $Header: /home/cvs/bp/oofem/oofemlib/src/simplecrosssection.h,v 1.11 2003/04/14 16:00:47 bp Exp $ */
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


//   **********************************
//   *** CLASS SIMPLE CROSSSECTOION ***
//   **********************************

#ifndef simplecrosssection_h 

#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"

class GaussPoint ;

/**
 Class implementing "simple" cross section model in finite element problem.
 A cross section  is an attribute of a domain. It is usually also attribute of many 
 elements.

 The simple cross section implementation does not perform any integration over cross-section volume,
 it repsesents a cross section model, where the whole cross section model is represented by single integration point.
 and therefore all requests for characteristic contibutions (stiffness) and for real stress computations are simply
 passed to parent StructuralCrossSection class, which invokes corresponding material mode services.
 Please note, that it is assumed that material model will support these material modes and provide
 corresponding services for characteristic components and stress evaluation.
 For description, how to incorporate more elaborate models of cross section, please read 
 base CrossSection docomentation.

 The overloaded methods giveFullCharacteristicVector and giveFullCharacteristicVector add some additional support
 for integrated cross section models - _3dShell, _3dBeam, _2dPlate and _2dBeam. 

 This class also reads into its property dictionary necessary geometric cross section characteristics,
 which are requested by particular material models. These parameters can be requested using get service and 
 include THICKNESS,WIDTH,AREA,INERTIA_MOMENT_Y,INERTIA_MOMENT_Z,TORSION_MOMENT_X and BEAM_SHEAR_COEFF values.
*/
class SimpleCrossSection : public StructuralCrossSection
{
/*
   This class implements a simple cross section in a finite element problem. A cross
 section  is an attribute of a domain. It is usually also attribute of many 
 elements. 

 DESCRIPTION
   
   Simple cross section represent a cross section in which stress strain state is
 described by using only one material point (GaussPoint). It includes also 
 methods for analysis plates, shels and beams (integral material point models).
 It does not include layered, fibred and others integral models. 

 TASK
 - Returning standard material stiffness marices (like 3dstress-strain, 2d plane , 
   plate, 3dbeam, 2d beam ..) according to current state determined by parametr
 StressMode by calling gp->material->GiveMaterialStiffnessMatrix (....) and by
 possible modifiing returned matrix. (for example in plate bending problems 
 by integrating over thickness ).
 - Returning RealStress state in gauss point and for given Stress mode.
 - Returning a properties of cross section like thickness or area.
*/
protected:

public:
 /// Constructor - creates SimpleCrossSection instance with number n and associates it with domain d.
 SimpleCrossSection (int n,Domain* d) : StructuralCrossSection (n,d)
  { }

 /**
  Computes full form of stress/strain from its reduced form, based on stress/strainn mode 
  stored in given integration point.
  Adds support for _3dShell, _3dBeam, _2dBeam and _2dPlate modes, for other modes parent 
  StructuralCrossSection service is invoked.
  @param answer full form of stress/strain vector
  @param gp integration point
  @param strainVector reduced vector
  */
 void  giveFullCharacteristicVector (FloatArray& answer,  GaussPoint*, const FloatArray&) ;
 /**
  Computes reduced stress/strain vector from full stress/strain vector.
  The stress/strain mode is determined form given integration point.
  Adds support for _3dShell, _3dBeam, _2dBeam and _2dPlate modes, for other modes parent 
  StructuralCrossSection service is invoked.
  @param answer charVector3d reduced
  @param gp integration point
  @param charVector3d full 3d stress/strain  vector
  */
 void  giveReducedCharacteristicVector (FloatArray& answer, GaussPoint* gp, 
                     const FloatArray& charVector3d);
/**  
 Computes the real stress vector for given strain and integration point. 
 The service should use previously reached equilibrium history variables. Also
 it should update temporary history variables in status according to newly reached state.
 The temporary history variables are moved into equilibrium ones after global structure
 equlibrium has been reached by iteration process.
 Elements should always pass their requests to their cross section model, which 
 performs necessary integration over its volume and invokes necessary material
 services for corresponding material model defined for given integration point.
 Implementation tests if material mode is supported by corresponding material model
 and pass the request to StructuralCrossSection giveRealStresses method, otherwise
 error is generated.
 @param answer contains result
 @param form material response form
 @param gp integration point 
 @param reducedStrainIncrement strain increment vector in reduced form
 @param tStep current time step (most models are able to respond only when atTime is current time step)
 */ 
 void  giveRealStresses(FloatArray &answer, MatResponseForm, GaussPoint*, 
             const FloatArray&, TimeStep*);

 // next function is intended to be used if we would like to obtain 
 // char matrix form different material which is not associated with gp and its element.
 // (mainly for obtaining linear elastic matrix)
 // stress-strain mode is taken from gp.
 // NORMALLY - PLEASE USE GiveCharMaterialStiffnessMatrix function
 //
/**
 Computes the stiffness matrix of receiver in given integration point, respecting its history.
 The algorithm should use temporary or equlibrium  history variables stored in integration point status
 to compute and return required result. 
 Elements should always pass their requests to their cross section model, which 
 performs necessary integration over its volume and invokes necessary material
 services for corresponding material model defined for given integration point.
 Current implementaion calls directly giveMaterialStiffnessMatrixOf service of receiver.
 @param answer contains result
 @param form material response form
 @param mode  material response mode
 @param gp integration point
 @param atTime time step (most models are able to respond only when atTime is current time step)
 */
   virtual void
  giveCharMaterialStiffnessMatrixOf  (FloatMatrix& answer,
                    MatResponseForm form, MatResponseMode rMode, 
                    GaussPoint*, StructuralMaterial*,
                    TimeStep*tStep) ;


/**
 Computes strain vector not dependent on sresses in given integration point. Returned vector is
 generated by temperature or shrinkage effects, for example.
 The load mode (Incremental or Total Load form) passed as parameter is taken into account.
 Depend on load form, tre resulting strain is total strain or its increment from previous 
 step. Overloaded to support beams, paltes and shells.
 @param answer stress independent strain vector
 @param gp integration point
 @param mode determines load mode
 @param stepN time step (most models are able to respond only when atTime is current time step)
 @param mode determines the response mode.
 */
 virtual void computeStressIndependentStrainVector (FloatArray& answer,
                           GaussPoint *gp, TimeStep *stepN, ValueModeType mode);


 // updates gp - record
 // stressMode is stored in gp

    /*FloatMatrix*    
   GiveCharMaterialStiffnessMatrix  (GaussPoint*, FloatArray *strainIncrement = NULL) ;*/

 /**
  Returns the value of cross section property 'aProperty'. Property must be identified 
  by unique int id.
  Implementation adds support for THICKNESS,WIDTH,AREA,INERTIA_MOMENT_Y,INERTIA_MOMENT_Z,TORSION_MOMENT_X 
  and BEAM_SHEAR_COEFF properties.
  @aProperty id of peroperty requested
  @return property value
  */
 double   give (int) ;
 
// identification and auxiliary functions
 /// Returns "SimpleCrossSection" - class name of the receiver.
 const char*    giveClassName () const {return "SimpleCrossSection" ;}
 /// Returns SimpleCrossSectionClass - classType id of receiver.
 classType giveClassID ()         const {return SimpleCrossSectionClass;}
 /// Returns input record name of the receiver.
 const char*    giveInputRecordName() const {return "SimpleCS";}
 /**
  Initializes receiver acording to object description stored in input record.
  Calls CrossSection initializeFrom service and reads the values of 
  THICKNESS,WIDTH,AREA,INERTIA_MOMENT_Y,INERTIA_MOMENT_Z,TORSION_MOMENT_X 
  and BEAM_SHEAR_COEFF parameters.
  */
 IRResultType initializeFrom (InputRecord* ir);
 /** Setups the input record string of receiver
  @param str string to be filled by input record
  @param keyword print record keyword (default true)
  */
 virtual int giveInputRecordString(std::string &str, bool keyword = true);

protected:
/**
 Evaluates the stiffness matrix of receiver. It is implemented as simple interface to material class, 
 forcing returned matrix to be in reduced form if material mode is supported by corresponding 
 material model (obtained from given integration point), 
 otherwise error is generated.
 For internal usage by cross section model. 
 It is direct interface to material model service giveCharacteristicMatrix.
 @param answer contains result
 @param form material response form
 @param mode  material response mode
 @param gp integration point
 @param mat pointer to material model
 @param tStep time step (most models are able to respond only when atTime is current time step)
*/
 void giveMaterialStiffnessMatrixOf (FloatMatrix& answer,
                   MatResponseForm form, 
                   MatResponseMode rMode,
                   GaussPoint* gp, 
                   StructuralMaterial* mat,
                   TimeStep* tStep);
 //
} ;


#define simplecrosssection_h
#endif
 
 
