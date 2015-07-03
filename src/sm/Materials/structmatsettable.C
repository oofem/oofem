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
#include "../sm/CrossSections/structuralcrosssection.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Material( StructuralMaterialSettable );

StructuralMaterialSettable :: StructuralMaterialSettable(int n, Domain *d) :
    StructuralMaterial(n, d)
{
    isoLE = new IsotropicLinearElasticMaterial(n,d);
}

StructuralMaterialSettable :: ~StructuralMaterialSettable()
{
    delete isoLE;
}

IRResultType
StructuralMaterialSettable :: initializeFrom(InputRecord *ir)
{
    //IRResultType result;                // Required by IR_GIVE_FIELD macro
    StructuralMaterial :: initializeFrom(ir);
    return isoLE->initializeFrom(ir);
}

void
StructuralMaterialSettable :: giveRealStressVector_3d(FloatArray &answer,
                                                  GaussPoint *gp,
                                                  const FloatArray &totalStrain,
                                                  TimeStep *atTime)
{

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    const FloatArray& stressVector = status->giveStressVector();

    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stressVector);
    answer = stressVector;
    return;
}


// TODO
void
StructuralMaterialSettable :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *atTime)
{
    isoLE->give3dMaterialStiffnessMatrix(answer,mode,gp,atTime);
}

MaterialStatus *
StructuralMaterialSettable :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, StructuralMaterial :: giveDomain(), gp);
}


} // end namespace oofem
