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



#include "domain.h"
#include "mixedpressurematerialextensioninterface.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "../sm/Materials/structuralmaterial.h"

#include <list>


namespace oofem {
// constructor
MixedPressureMaterialExtensionInterface :: MixedPressureMaterialExtensionInterface(Domain *d)  : Interface()
{
    dom = d;
}




void
MixedPressureMaterialExtensionInterface :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, double pressure, TimeStep *tStep)
{
    ///@todo Move this to StructuralCrossSection ?
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        this->giveRealStressVector_3d(answer, gp, reducedStrain, pressure, tStep);
    } else if ( mode == _PlaneStrain ) {
        this->giveRealStressVector_PlaneStrain(answer, gp, reducedStrain, pressure, tStep);
    } else {
        OOFEM_ERROR("Unknown material mode for the mixed u-p formulation");
    }
}



void
MixedPressureMaterialExtensionInterface ::  giveDeviatoricConstitutiveMatrix(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep)
{
    ///@todo Move this to StructuralCrossSection ?
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        this->giveDeviatoric3dMaterialStiffnessMatrix(answer, rmode, gp, tStep);
    } else if ( mode == _PlaneStrain ) {
        this->giveDeviatoricPlaneStrainStiffMtrx(answer, rmode, gp, tStep);
    } else {
        OOFEM_ERROR("Unknown material mode for the mixed u-p formulation");
    }
}




void
MixedPressureMaterialExtensionInterface :: giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, double pressure, TimeStep *tStep)
{
    FloatArray vE, vS;
    StructuralMaterial :: giveFullSymVectorForm(vE, reducedStrain, _PlaneStrain);
    this->giveRealStressVector_3d(vS, gp, vE, pressure, tStep);
    StructuralMaterial :: giveReducedSymVectorForm(answer, vS, _PlaneStrain);
}
}
