/* $Header: /home/cvs/bp/oofem/sm/src/beam2d.h,v 1.5 2003/04/06 14:08:30 bp Exp $ */
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

//   ********************
//   *** CLASS Beam2d ***
//   ********************

#ifndef beam2d_h
#define beam2d_h

#include "structuralelement.h"
#include "layeredcrosssection.h"

class Beam2d : public  StructuralElement, public LayeredCrossSectionInterface
{
/*
   This class implements a 2-dimensional beam element
 with cubic lateral displace,ent interpolation (rotations are quadratic)
 and longitudial displacements are linear.
 This is an exact displacement approximation for beam with no
 nonnodal loading.

 This class is not derived from liBeam2d or truss element, because it does not support
 any material nonlinearities (if shoul, stiffness must be integrated)
*/
protected:
  double kappa,pitch, length;
  IntArray *dofsToCondense;
  
public:
  Beam2d (int,Domain*) ;                       // constructor
  ~Beam2d () ;                                 // destructor

  void          computeConsistentMassMatrix (FloatMatrix& answer, TimeStep* tStep, double &mass);
  void          computeInitialStressMatrix (FloatMatrix& answer, TimeStep* tStep) ;
  /*
  Returns mask indicating, which unknowns (their type and ordering is the same as in vector 
  of receivers unknowns, interpolated by N matrix. Nonzero value at i-th position 
  indicates that corresponding row in interpolation matrix N will participate in
  mass matrix integration (typically only displacements are taken into account).
  */
  //virtual void          giveMassMtrxIntegrationgMask (IntArray& answer) ;

  void          computeStiffnessMatrix (FloatMatrix& answer,
                    MatResponseMode rMode, TimeStep* tStep) ;
  int           giveLocalCoordinateSystem (FloatMatrix& answer);
  void          computeLocalForceLoadVector (FloatArray& answer, TimeStep* stepN, ValueModeType mode);
  void          giveInternalForcesVector (FloatArray& answer, 
                     TimeStep*, int useUpdatedGpRecord = 0) ;
  void          giveEndForcesVector  (FloatArray& answer, TimeStep* tStep);
 /**
  Computes the global coordinates from given element's local coordinates.
  @returns retutns nonzero if successful
  */
 virtual int computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) ;


  virtual int testElementExtension (ElementExtension ext) {return ((ext==Element_EdgeLoadSupport)?1:0);}
  //int hasLayeredSupport () {return 1;}

  /** Interface requesting service */
  Interface* giveInterface (InterfaceType);

 virtual int            computeNumberOfDofs (EquationID ut) {return 6;}
 virtual void           giveDofManDofIDMask  (int inode, EquationID, IntArray& ) const;
  double        computeVolumeAround (GaussPoint*) ;

   void  printOutputAt (FILE* , TimeStep* );
// 
// definition & identification
//
  const char* giveClassName () const { return "Beam2d" ;}
  classType            giveClassID () const { return Beam2dClass; } 
  IRResultType initializeFrom (InputRecord* ir);

#ifdef __OOFEG
      void          drawRawGeometry (oofegGraphicContext&);
      void          drawDeformedGeometry(oofegGraphicContext&, UnknownType);
#endif


/**
 @name The element interface required by LayeredCrossSectionInterface
*/
//@{
/**
 Computes full 3d strain vector in element layer. This function is necesary 
 if layered cross section is specified. If it is implemented, the testElementExtension
 servise should return nonzero for Element_LayeredSupport parameter. This service is used by
 layered cross section models.
 @param answer full layer starin vector
 @param masterGp element integration point
 @param slaveGp slave integration point representing particular layer
 @tStep time step
*/
 void  computeStrainVectorInLayer (FloatArray& answer, GaussPoint* masterGp, 
                  GaussPoint* slaveGp, TimeStep* tStep);
//@}

protected:
  void          computeEdgeLoadVectorAt (FloatArray &answer, Load*, int, TimeStep*, ValueModeType mode) ;
  void          computePrescribedStrainLocalLoadVectorAt (FloatArray& answer, TimeStep* tStep, ValueModeType mode) ;
  //void          computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);
  void          computeBmatrixAt (GaussPoint*, FloatMatrix&, int=1, int=ALL_STRAINS) ;
  void          computeNmatrixAt (GaussPoint*, FloatMatrix &)  ;
  int           computeGtoLRotationMatrix (FloatMatrix&);  // giveRotationMatrix () ;
//  int           computeGtoNRotationMatrix (FloatMatrix&);

  double giveKappaCoeff () ;
  double        giveLength () ;
  double        givePitch () ;
  virtual void computeClampedStiffnessMatrix (FloatMatrix& answer, 
                       MatResponseMode rMode, TimeStep* tStep) ;
  virtual void computeLocalStiffnessMatrix (FloatMatrix& answer, 
                       MatResponseMode rMode, TimeStep* tStep) ;
  void          computeGaussPoints () ;
 integrationDomain  giveIntegrationDomain () {return _Line;}
  // return desired number of integration points for consistent mass matrix 
  // computation, if required.
  virtual int  giveNumberOfIPForMassMtrxIntegration () {return 4;}
  
} ;

#endif
