/* $Header: $ */
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

#ifndef randomfieldgenerator_h
#define randomfieldgenerator_h
#include "gausspnt.h"
#include "flotarry.h"

namespace oofem {

class RandomFieldGenerator
{
protected:
  Domain* domain;
public:
  /// Constructor. Creates empty RandomFieldGenerator
  RandomFieldGenerator (int n, Domain *d) {this->domain = d;}
  /// Destructor
  virtual ~RandomFieldGenerator () {}
  /**
     Generates random value
  */
  virtual void generateRandomValue(double& value, FloatArray *position){;}
  virtual void generateRandomValueAt(double& value, GaussPoint *gp) {
    FloatArray globalCoordinates;
    if (gp->giveElement()->computeGlobalCoordinates(globalCoordinates,*(gp->giveLocalCoordinates())))
      this->generateRandomValue (value, &globalCoordinates);
    else
      OOFEM_ERROR ("RandomFieldGenerator::generateRandomValue computeGlobalCoordinates failed");
  }
  
  virtual IRResultType initializeFrom(InputRecord *ir) {return IRRT_OK;}
  /// Returns class name of the receiver.
  virtual const char *giveClassName() const { return "RandomFieldGenerator"; }
 
};

} // end namespace oofem
#endif
