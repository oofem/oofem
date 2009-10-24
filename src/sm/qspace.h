/* $Header: /home/cvs/bp/oofem/sm/src/qspace.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ************************************
//   *** CLASS QUADRATIC 3D Element   ***
//   ************************************

//
//  Recasted by L. Svoboda
//
#ifndef qspace_h
#define qspace_h

#include "structuralelement.h"
#include "fei3dhexaquad.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "eleminterpmapperinterface.h"
#include "huertaerrorestimator.h"

class QSpace : public StructuralElement, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface

{
  /*
    This class implements an Quadratic 3d  20 - node
    elasticity finite element. Each node has 3 degrees of freedom.
    DESCRIPTION :
    One single additional attribute is needed for Gauss integration purpose :
    'jacobianMatrix'. This 3x3 matrix contains polynomials.
    TASKS :
    - calculating its Gauss points ;
    - calculating its B,D,N matrices and dV.
  */

 protected:
  int numberOfGaussPoints;
  static FEI3dHexaQuad interpolation;

 public:
  QSpace (int,Domain*) ;     // constructor
  ~QSpace ()  {}             // destructor

  IRResultType initializeFrom (InputRecord* ir);
  virtual void giveDofManDofIDMask (int inode, EquationID ut, IntArray& answer) const;
  double       computeVolumeAround (GaussPoint*) ;

  /**
   Computes the global coordinates from given element's local coordinates.
   Required by nonlocal material models.
   @returns nonzero if successful
  */
  virtual int  computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) ;

  /// characteristic length in gp (for some material models)
  double       giveCharacteristicLenght (GaussPoint* gp, const FloatArray &normalToCrackPlane);

  /// Interface requesting service
  Interface*  giveInterface (InterfaceType);
  virtual int testElementExtension (ElementExtension ext) {return ((ext==Element_SurfaceLoadSupport)?1:0);}
  /**
  @name The element interface required by ZZNodalRecoveryModel
  */
  //@{
  /**
   Returns the size of DofManger record required to hold recovered values for given mode.
   @param type determines the type of internal variable to be recovered
   @return size of DofManger record required to hold recovered values
  */
  int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);

  /// Returns the corresponding element to interface
  Element *ZZNodalRecoveryMI_giveElement() { return this; }

  /// Evaluates N matrix (interpolation estimated stress matrix).
  void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx (FloatMatrix& answer, GaussPoint* aGaussPoint, InternalStateType type);
  //@}
  /**
  @name The element interface required by NodalAveragingRecoveryModel
  */
  //@{
  /**
   Computes the element value in given node.
   @param answer contains the result
   @param node element node number
   @param type determines the type of internal variable to be recovered
   @param tStep time step
  */
  void NodalAveragingRecoveryMI_computeNodalValue (FloatArray& answer, int node, InternalStateType type, TimeStep* tStep);

  /**
   Computes the element value in given side.
   @param answer contains the result
   @param side element side number
   @param type determines the type of internal variable to be recovered
   @param tStep time step
  */
  void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);

  /**
   Returns the size of DofManger record required to hold recovered values for given mode.
   @param type determines the type of internal variable to be recovered
   @return size of DofManger record required to hold recovered values
  */
  virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
  //@}

  //
  // definition & identification
  //
  const char* giveClassName () const { return "QSpace"; }
  classType   giveClassID   () const { return QSpaceClass; }
  Element_Geometry_Type giveGeometryType() const {return EGT_hexa_2;}
  virtual int  computeNumberOfDofs (EquationID ut) {return 60;}

 protected:
  void computeGaussPoints ();
  void computeNmatrixAt (GaussPoint* ,FloatMatrix& );
  void computeBmatrixAt (GaussPoint* ,FloatMatrix& ,int=1,int=ALL_STRAINS);

  integrationDomain  giveIntegrationDomain () {return _Cube;}
  int giveApproxOrder () {return 2;}
  int giveNumberOfIPForMassMtrxIntegration () {return 27;}

  //@}
  /**
     @name Surface load support
  */
  //@{
  virtual IntegrationRule* GetSurfaceIntegrationRule (int) ;
  virtual void   computeSurfaceNMatrixAt (FloatMatrix& answer, GaussPoint*) ;
  virtual void   giveSurfaceDofMapping (IntArray& answer, int) const;
  virtual double computeSurfaceVolumeAround (GaussPoint*, int) ;
  virtual void   computeSurfIpGlobalCoords (FloatArray& answer, GaussPoint*, int) ;
  virtual int    computeLoadLSToLRotationMatrix (FloatMatrix& answer, int, GaussPoint*) ;
  //@}
} ;

#endif
