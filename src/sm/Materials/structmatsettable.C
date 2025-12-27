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

#include "structmatsettable.h"

#include "floatarray.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Material( StructuralMaterialSettable );

StructuralMaterialSettable :: StructuralMaterialSettable(int n, Domain *d) :
    StructuralMaterial(n, d),
    isoLE(n,d)
{}


void
StructuralMaterialSettable :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);
    isoLE.initializeFrom(ir);
}


FloatArrayF<6>
StructuralMaterialSettable :: giveRealStressVector_3d(const FloatArrayF<6> &strain,GaussPoint *gp, 
                                                      TimeStep *atTime) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    const auto &stressVector = status->giveStressVector();

    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(stressVector);
    return stressVector;
}


// TODO
FloatMatrixF<6,6>
StructuralMaterialSettable :: give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp,
                                                            TimeStep *atTime) const
{
    return isoLE.give3dMaterialStiffnessMatrix(mode, gp, atTime);
}

std::unique_ptr<MaterialStatus> 
StructuralMaterialSettable :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<StructuralMaterialStatus>(gp);
}


} // end namespace oofem
