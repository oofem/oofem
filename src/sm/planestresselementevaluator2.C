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

#include "planestresselementevaluator2.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "domain.h"
#include "node.h"
#include "elementgeometry.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "iga.h"

namespace oofem {


  PlaneStressStructuralElementEvaluator2 :: PlaneStressStructuralElementEvaluator2(int dispInterpNumber):StructuralElementEvaluator2(dispInterpNumber)
  {
    displacementInterpolationNumber = dispInterpNumber;
    numberOfDofs = 2;
    
  }
 PlaneStressStructuralElementEvaluator2 :: PlaneStressStructuralElementEvaluator2():StructuralElementEvaluator2()
  {
        
  }

  void PlaneStressStructuralElementEvaluator2 :: computeBmatrixAt( GaussPoint *gp,FloatMatrix &answer, int lowerIndx, int upperIndx) {
    int i, nDofMan;
    //IntArray dofmanSubElementMask;
    FloatMatrix d;
    nDofMan  = this->giveElement()->giveNumberOfDofManagers();
    FEInterpolation *interp = gp->giveElement()->giveInterpolation(displacementInterpolationNumber);
    // this uses FEIInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx(d, * gp->giveCoordinates(), FEIElementGeometryWrapper( gp->giveElement()));

    answer.resize(3, nDofMan);
    answer.zero();

    for ( i = 1; i <= nDofMan/2.; i++ )
      {
	answer.at(1, i * 2 - 1) = d.at(i, 1);
	answer.at(2, i * 2 - 0)   = d.at(i, 2);
	
	answer.at(3, 2 * i - 1) = d.at(i, 2);
	answer.at(3, 2 * i - 0) = d.at(i, 1);
      }
  }
  

  
  double PlaneStressStructuralElementEvaluator2 :: computeVolumeAround(GaussPoint *gp,ElementGeometry* elem) {
    double determinant, weight, thickness, volume;
    determinant = fabs( elem->giveInterpolation(displacementInterpolationNumber)->giveTransformationJacobian(* gp->giveCoordinates(),FEIElementGeometryWrapper( elem)) );
    weight      = gp->giveWeight();
    thickness   = elem->giveCrossSection()->give(CS_Thickness);
    volume      = determinant * weight * thickness;

    return volume;
}








} // end namespace oofem
