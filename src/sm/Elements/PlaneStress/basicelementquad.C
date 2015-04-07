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
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "Elements/PlaneStress/basicelementquad.h"
#include "fei2dquadlin.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(BasicElementQuad);

// Interpolator describing shape functions for the approximated unknowns.
// 1 -> first spatial index, 2 -> second spatial index
FEI2dQuadLin BasicElementQuad :: interp(1, 2);    

BasicElementQuad :: BasicElementQuad(int n, Domain *aDomain) : PlaneStressElement(n, aDomain)
{
    this->numberOfGaussPoints = 4;
}


FEInterpolation *BasicElementQuad :: giveInterpolation() const 
{ 
    /* Returns the interpolator used for the element which provide
     * shape functions, their derivatives, area of the element etc.
     */
    return & interp; 
}


} // end namespace oofem
