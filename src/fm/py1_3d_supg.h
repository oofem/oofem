/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.h,v 1.3 2003/04/23 14:22:15 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2002   Borek Patzak                                       



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

//   ************************************************************************************************
//   *** 3D LINEAR PYRAMID ELEMENT FOR FLUID DYNAMIC PROBLEMS SOLVED WITH SUPG/PSSPG ALGORITHM ***
//   ************************************************************************************************
 
#ifndef py1_3d_supg_h


#include "supgelement2.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "levelsetpcs.h"
#include "fei3dtrlin.h"



class TimeStep ; class Node ; class Material ; class GaussPoint ;
class FloatMatrix ; class FloatArray ; class IntArray ;

/**
  Class representing 3d linear pyramid element
  for solving incompressible fluid with SUPG solver

*/
class PY1_3D_SUPG : public SUPGElement2, public LevelSetPCSElementInterface
{
protected:
  static FEI3dTrLin interpolation;
public:
  // constructor
  PY1_3D_SUPG (int,Domain*) ;
  ~PY1_3D_SUPG () ;                        // destructor
  
  // definition
  const char* giveClassName () const { return "PY1_3D_SUPG" ;}
  classType                giveClassID () const { return SUPGElementClass; }
  Element_Geometry_Type giveGeometryType() const {return EGT_tetra_1;}
  
  virtual void           giveElementDofIDMask  (EquationID, IntArray& answer) const;
  virtual void           giveDofManDofIDMask  (int inode, EquationID ut, IntArray& answer) const ;
  virtual int            computeNumberOfDofs (EquationID ut);
  IRResultType           initializeFrom (InputRecord* ir);

  /** Interface requesting service */
  Interface* giveInterface (InterfaceType);

  double                computeCriticalTimeStep (TimeStep* tStep) ;
  double                computeVolumeAround (GaussPoint*);

  /**
     @name The element interface required by LevelSetPCSElementInterface
  */
  //@{
  /** Evaluetes F in level set equation of the form
      fi_t+F(grad(fi), x)*norm(grad(fi)) = 0
      where for interface position driven by flow with speed u: 
      F=integral of {dotProduct(u,grad(fi))/norm(grad(fi))}
  */
  double LS_PCS_computeF (LevelSetPCS*, TimeStep*);
  /// Returns gradient of shape functions (assumed constatnt <- linear approx)
  void LS_PCS_computedN (FloatMatrix& answer);
  /// Returns receiver's volume
  double LS_PCS_computeVolume() ;
  virtual double LS_PCS_computeS (LevelSetPCS*, TimeStep*);
  void LS_PCS_computeVOFFractions (FloatArray& answer, FloatArray& fi);
  //@}


 protected:
  void                  computeGaussPoints ();
  virtual void computeNuMatrix (FloatMatrix& answer, GaussPoint* gp);
  virtual void computeUDotGradUMatrix (FloatMatrix& answer, GaussPoint* gp, TimeStep* atTime);
  virtual void computeBMatrix (FloatMatrix& anwer, GaussPoint* gp);
  virtual void computeDivUMatrix (FloatMatrix& answer, GaussPoint* gp);
  virtual void computeNpMatrix (FloatMatrix& answer, GaussPoint* gp);
  virtual void computeGradPMatrix (FloatMatrix& answer, GaussPoint* gp);
  virtual int  giveNumberOfSpatialDimensions ();
} ;

#define py1_3d_supg_h
#endif







