/* $Header: /home/cvs/bp/oofem/sm/src/QSpace.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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
#ifndef qspacegrad_h
#define qspacegrad_h
#include "qspace.h"
#include "fei3dhexalin.h"
#include "graddpelement.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "eleminterpmapperinterface.h"
#include "huertaerrorestimator.h"
#include "sprnodalrecoverymodel.h"


namespace oofem {

class QSpaceGrad : public QSpace,public GradDpElement {
 protected:
  int numberOfGaussPoints;
  static FEI3dHexaLin interpolation;

 public:
  QSpaceGrad (int,Domain*) ;     // constructor
  ~QSpaceGrad ()  {}             // destructor

IRResultType initializeFrom (InputRecord* ir);
  virtual void giveDofManDofIDMask (int inode, EquationID ut, IntArray& answer) const;
  //
  // definition & identification
  //
  const char* giveClassName () const { return "QSpaceGrad"; }
  classType   giveClassID   () const { return QSpaceGradClass; }
  virtual int  computeNumberOfDofs (EquationID ut) {return 68;}
  
 protected:
 ///////////////////////////////////////////////////////////////////////////////
   void computeGaussPoints ();
   void computeNkappaMatrixAt(GaussPoint*,FloatMatrix&);
   void computeBkappaMatrixAt(GaussPoint*,FloatMatrix&);
   StructuralElement* giveStructuralElement(){return this;}
};
 
}
#endif // end namespace oofem

