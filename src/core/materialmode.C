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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "materialmode.h"

namespace oofem {

bool mmodeIs1D(MaterialMode mmode){
    return (mmode == _1dMat) || (mmode == _1dMatGrad) || (mmode == _1dInterface) || (mmode == _1dHeat) || (mmode == _1dHeMo) || (mmode == _1dUP) || (mmode == _1dLattice);
}
bool mmodeIs2D(MaterialMode mmode){
    return (mmode == _PlaneStress) || (mmode == _PlaneStrain) || (mmode == _PlaneStressGrad) || (mmode == _PlaneStrainGrad) || (mmode == _2dPlate) || (mmode == _2dPlateSubSoil) || (mmode == _2dBeam) || (mmode == _2dInterface) || (mmode == _2dHeat) || (mmode == _2dHeMo) || (mmode == _2dFlow) || (mmode == _2dAxiFlow) || (mmode == _2dUP) || (mmode == _2dUPV) || (mmode == _2dLattice) || (mmode == _2dMTLattice);  
}
bool mmodeIs3D(MaterialMode mmode){
    return (mmode == _3dMat) || (mmode == _3dMatGrad) || (mmode == _3dBeam) || (mmode == _3dShell) || (mmode == _3dShellRot) || (mmode == _3dDegeneratedShell) || (mmode == _3dInterface) || (mmode == _3dHeat) || (mmode == _3dHeMo) || (mmode == _3dFlow) || (mmode == _3dUP) || (mmode == _3dUPV) || (mmode == _3dLattice) || (mmode == _3dMTLattice) || (mmode == _Warping);
}

} // end namespace oofem
