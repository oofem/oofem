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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef latticestructuralelement_h
#define latticestructuralelement_h

#include "structuralelement.h"

namespace oofem {

class LatticeStructuralElement : public StructuralElement
{
  /**
  *  This class implements a base lattice element    
  *
  */

 public:
  LatticeStructuralElement (int,Domain*) ;                       // constructor
  ~LatticeStructuralElement () ;                                 // destructor
  
  IRResultType initializeFrom (InputRecord* ir);

  virtual  double giveArea() {return 0;}
  virtual  double giveLength() {return 0;}
  virtual  int giveCrackFlag() {return 0;}
  virtual  int giveNumberOfCrossSectionNodes() {return 0;}
  virtual void giveCrossSectionCoordinates(FloatArray & coords){};
  virtual  double giveCrackWidth() {return 0;}
  virtual  void giveGPCoords(FloatArray & gpcoords) { }
  virtual  double giveDissipation() {return 0;}
  virtual  double giveDeltaDissipation() {return 0;}
 

};
} // end namespace oofem
#endif
