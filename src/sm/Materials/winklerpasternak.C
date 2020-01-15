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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#include "winklerpasternak.h"
#include "sm/Materials/structuralms.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(WinklerPasternakMaterial);

WinklerPasternakMaterial:: WinklerPasternakMaterial (int n, Domain* d): StructuralMaterial(n, d) 
{ }


void
WinklerPasternakMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, c1, _IFT_WinklerPasternakMaterial_C1);
    if ( ir.hasField(_IFT_WinklerPasternakMaterial_C2)) {
      // isotropic case
      IR_GIVE_FIELD(ir, c2x, _IFT_WinklerPasternakMaterial_C2);
      c2y = c2x;
    } else {
      IR_GIVE_FIELD(ir, c2x, _IFT_WinklerPasternakMaterial_C2X);
      IR_GIVE_FIELD(ir, c2y, _IFT_WinklerPasternakMaterial_C2Y);
    }
}


void
WinklerPasternakMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);

    input.setField(this->c1, _IFT_WinklerPasternakMaterial_C1);
    input.setField(this->c2x, _IFT_WinklerPasternakMaterial_C2X);
    input.setField(this->c2y, _IFT_WinklerPasternakMaterial_C2Y);
}


FloatArrayF<3>
WinklerPasternakMaterial::giveRealStressVector_2dPlateSubSoil(const FloatArrayF<3> &reducedE, GaussPoint *gp, TimeStep *tStep) const
{
    auto tangent = this->give2dPlateSubSoilStiffMtrx(ElasticStiffness, gp, tStep);
    auto answer = dot(tangent, reducedE);

    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    status->letTempStrainVectorBe(reducedE);
    status->letTempStressVectorBe(answer);
    return answer;
}


FloatMatrixF<3,3>
WinklerPasternakMaterial :: give2dPlateSubSoilStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    return diag<3>({c1, c2x, c2y});
}


MaterialStatus *
WinklerPasternakMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(gp);
}


bool
WinklerPasternakMaterial :: hasMaterialModeCapability(MaterialMode mode) const
//
// returns whether receiver supports given mode
//
{
    return mode == _2dPlateSubSoil;
}

} // end namespace oofem
