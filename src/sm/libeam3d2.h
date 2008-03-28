/* $Header: /home/cvs/bp/oofem/sm/src/libeam3d2.h,v 1.6 2003/04/06 14:08:30 bp Exp $ */
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

//   **********************
//   *** CLASS LIBeam3d2 ***
//   **********************

#ifndef libeam3d2_h
#define libeam3d2_h

#include "nlstructuralelement.h"
#include "fiberedcs.h"

class LIBeam3d2 : public  NLStructuralElement, public FiberedCrossSectionInterface
{
/*
   This class implements a 3-dimensional Linear Isoparametric
   Mindlin theory beam element, with reduced integration.
  VERY SIMPLE support for geometric nonlinearity is taken into account.
  The central triad (which defines local coordinate system - corotational system)
  is updated using average increment of rotational pseudovector (same as in Simo -Qu Voc element).
  Nonlinearity in geometric equations is not taken into account, 
  the only nonlinearity accounted is due to rotation of local (corotational) frame.
  Increments of eleemnt's strains are computed using increment of displacements in
  updated local system and total strain by adding this increnet to stored total deformation.
  The local stiffness or nodal forces are transformed using updated centre triad to global system.
*/
  
private:
  double length;
  int referenceNode;


  // parameters related to nonlinearity
  //
  /// last equlibrium triad at the centre
  FloatMatrix tc;
  /// temporary triad at the centre
  FloatMatrix tempTc;
    /// time stamp of temporary centre triad
  StateCounterType tempTcCounter;


public :
    LIBeam3d2 (int,Domain*) ;                       // constructor
    ~LIBeam3d2 ()  {}                               // destructor

  // FloatMatrix*  ComputeConstitutiveMatrixAt (GaussPoint*) ;
  // FloatArray*   ComputeResultingBodyForceAt (TimeStep*) ;
  void          computeLumpedMassMatrix (FloatMatrix& answer, TimeStep* tStep) ;
  void          computeMassMatrix (FloatMatrix& answer, TimeStep* tStep) 
  {computeLumpedMassMatrix (answer, tStep);}
  void          computeStiffnessMatrix (FloatMatrix& answer, 
                    MatResponseMode rMode, TimeStep* tStep) ;
  int           computeGtoLRotationMatrix (FloatMatrix&);  // giveRotationMatrix () ;

  virtual int testElementExtension       (ElementExtension ext);

 virtual int            computeNumberOfDofs (EquationID ut) {return 12;}
 virtual void           giveDofManDofIDMask  (int inode, EquationID, IntArray& ) const;
  double        computeVolumeAround (GaussPoint*) ;
  int           giveLocalCoordinateSystem (FloatMatrix& answer);
  void          computeStrainVector (FloatArray& answer, GaussPoint* gp,TimeStep* tStep);
 /**
  Computes the global coordinates from given element's local coordinates.
  @returns nonzero if successful
  */
 virtual int computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) ;

//
// fibered cross section support functions
//
  void FiberedCrossSectionInterface_computeStrainVectorInFiber (FloatArray& answer, GaussPoint* masterGp, 
                                GaussPoint* slaveGp, TimeStep* tStep);


  /** Interface requesting service */
  Interface* giveInterface (InterfaceType);



// 
// definition & identification
//
  const char* giveClassName () const { return "LIBeam3d" ;}
  classType             giveClassID ()          const { return LIBeam3dClass; } 
  IRResultType initializeFrom (InputRecord* ir);

  // due to nonlinearity
  contextIOResultType           saveContext (DataStream* stream, ContextMode mode, void *obj);
  contextIOResultType           restoreContext (DataStream* stream, ContextMode mode, void *obj);

#ifdef __OOFEG
      void          drawRawGeometry (oofegGraphicContext&);
      void          drawDeformedGeometry(oofegGraphicContext&, UnknownType);
      void          drawScalar   (oofegGraphicContext& context);
      
#endif

protected:
  // edge load support
  void  computeEgdeNMatrixAt (FloatMatrix& answer, GaussPoint*) ;
  void  giveEdgeDofMapping (IntArray& answer, int) const;
  double        computeEdgeVolumeAround (GaussPoint*, int);
  void          computeEdgeIpGlobalCoords (FloatArray& answer, GaussPoint* gp, int iEdge) 
  {computeGlobalCoordinates(answer, *(gp->giveCoordinates()));}
  int   computeLoadLEToLRotationMatrix (FloatMatrix&, int, GaussPoint*) ;
 int  computeLoadGToLRotationMtrx (FloatMatrix& answer) ;
  //void          computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);
  void          computeBmatrixAt (GaussPoint*, FloatMatrix&, int=1, int=ALL_STRAINS) ;
  void          computeNmatrixAt (GaussPoint*, FloatMatrix &) ;
  void          computeGaussPoints () ;
 integrationDomain  giveIntegrationDomain () {return _Line;}
  double        giveLength () ;

  // nonlinerity
  void          updateTempTriad (TimeStep* tStep);
  void          computeSMtrx (FloatMatrix& answer, FloatArray& vec);
  void          computeRotMtrx (FloatMatrix& answer, FloatArray& psi);
  double        giveCurrentLength (TimeStep* tStep);
  void          updateYourself (TimeStep* tStep);
  void          initForNewStep ();
  int           updateRotationMatrix () ;
} ;

#endif // libeam3d2_h



