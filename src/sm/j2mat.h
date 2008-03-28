/* $Header: /home/cvs/bp/oofem/sm/src/Attic/j2mat.h,v 1.1.2.1 2004/04/05 15:19:47 bp Exp $ */
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

//   *******************************
//   *** CLASS J2 plastic material  
//   *******************************

#ifndef j2mat_h
#define j2mat_h

#include "mplasticmaterial2.h"

class Domain;

class J2Mat : public MPlasticMaterial2
{
/*
   This class implements a isotropic  plastic linear material (J2 plasticity condition is used)
 in a finite element problem. A material
   is an attribute of a domain. It is usually also attribute of many elements.

 DESCRIPTION
   ISOTROPIC PLASTIC Material with J2 plastic condition
 
 TASK
 - Returning standard material stiffness marix for 3d-case.
   according to current state determined by using data stored 
   in Gausspoint.
 - Returning a material property (method 'give'). Only for non-standard elements.
 - Returning real stress state vector(tensor) at gauss point for 3d - case.
*/

protected:
  int  kinematicHardeningFlag, isotropicHardeningFlag;
  double kinematicModuli, isotropicModuli;
  //double E, nu; // isotropic material constants
  double k;
public:
 
  J2Mat (int n, Domain* d) ;
  ~J2Mat ();
  
  IRResultType initializeFrom (InputRecord* ir);
  const char*    giveClassName () const      {return "J2Mat" ;}
  classType giveClassID ()         const      {return J2MatClass;}

  virtual int         giveSizeOfFullHardeningVarsVector();
  virtual int         giveSizeOfReducedHardeningVarsVector(GaussPoint*);
  virtual bool isCharacteristicMtrxSymmetric (MatResponseMode rMode) {return false;}

  MaterialStatus* CreateStatus (GaussPoint* gp) const ;

protected:
  virtual int giveMaxNumberOfActiveYieldConds (GaussPoint* gp) {return 2;}

  //
  // yield(YC-like functions) and loading(LC-like functions) criteria specific section
  //
  virtual double computeYieldValueAt (GaussPoint *gp, int isurf, const FloatArray& stressVector, 
                                      const FloatArray& strainSpaceHardeningVars);

  virtual void computeStressGradientVector (FloatArray& answer, functType ftype, int isurf, GaussPoint *gp, 
                                            const FloatArray& stressVector, const FloatArray& strainSpaceHardeningVars);

  virtual void computeStrainHardeningVarsIncrement (FloatArray& answer, GaussPoint* gp, 
                                                    const FloatArray& stress, const FloatArray& dlambda,
                                                    const FloatArray& dplasticStrain, const IntArray& activeConditionMap) ;
  virtual void computeKGradientVector (FloatArray& answer, functType ftype, int isurf, GaussPoint *gp, FloatArray& fullStressVector,
                                       const FloatArray& strainSpaceHardeningVariables);
  
  virtual void computeReducedHardeningVarsSigmaGradient(FloatMatrix& answer, GaussPoint* gp, const IntArray& activeConditionMap,
                                                        const FloatArray&fullStressVector, 
                                                        const FloatArray& strainSpaceHardeningVars,
                                                        const FloatArray& gamma) ;
  virtual void computeReducedHardeningVarsLamGradient(FloatMatrix& answer, GaussPoint* gp, int actSurf, 
                                                      const IntArray& activeConditionMap, 
                                                      const FloatArray&fullStressVector, 
                                                      const FloatArray& strainSpaceHardeningVars,
                                                      const FloatArray& gamma) ;
  /**
    Indicates, whether receiver model has hardening/softening behaviour or behaves according to perfect plasticity theory
  */
  virtual int hasHardening () ;
 /* virtual void  computeReducedGradientMatrix (FloatMatrix& answer, int isurf,
                                             GaussPoint *gp,
                                             const FloatArray& stressVector, 
                                             const FloatArray& stressSpaceHardeningVars) = 0;*/
 virtual void  computeReducedSSGradientMatrix (FloatMatrix& gradientMatrix,  int i, GaussPoint* gp , const FloatArray& fullStressVector, 
                                               const FloatArray& strainSpaceHardeningVariables);
 virtual void  computeReducedSKGradientMatrix (FloatMatrix& gradientMatrix,  int i, GaussPoint* gp , const FloatArray& fullStressVector, 
                                               const FloatArray& strainSpaceHardeningVariables);
  
 // auxiliary function
 double      computeJ2InvariantAt(const FloatArray&) ;
 double      giveIsotropicHardeningVar(GaussPoint* gp, const FloatArray& strainSpaceHardeningVars);
 void        giveStressBackVector (FloatArray& answer, GaussPoint* gp, const FloatArray& strainSpaceHardeningVars);
};

#endif // j2mat_h
