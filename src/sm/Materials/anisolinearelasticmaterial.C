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

#include "linearelasticmaterial.h"
#include "anisolinearelasticmaterial.h"
#include "structuralms.h"
#include "floatmatrix.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(AnisotropicLinearElasticMaterial);

void
AnisotropicLinearElasticMaterial :: initializeFrom(InputRecord &ir)
{
    LinearElasticMaterial :: initializeFrom(ir);

    // read the stiffness coefficients arranged by rows from the diagonal to the right (21 values)
    FloatArray stiffness;
    IR_GIVE_FIELD(ir, stiffness, _IFT_AnisotropicLinearElasticMaterial_stiff);
    if ( stiffness.giveSize() != 21 ) {
        OOFEM_ERROR( "Incorrect size of stiff - should be 21, is %d\n", stiffness.giveSize() );
    }

    // put the stiffness coefficients into a 6x6 matrix
    for ( int k = 1, i = 1; i <= 6; i++ ) {
        tangent.at(i, i) = stiffness.at(k++);
        for ( int j = i + 1; j <= 6; j++ ) {
            tangent.at(i, j) = tangent.at(j, i) = stiffness.at(k++);
        }
    }
    this->computesSubTangents();

    FloatArray alpha_input(6);
    IR_GIVE_FIELD(ir, alpha_input, _IFT_AnisotropicLinearElasticMaterial_talpha);
    if ( alpha_input.giveSize() != 6 ) {
        OOFEM_ERROR( "Incorrect size of talpha - should be 6, is %d\n", alpha.giveSize() );
    }
    alpha = alpha_input;
}


void
AnisotropicLinearElasticMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    Material :: giveInputRecord(input);
    FloatArray stiffness(21);
    for ( int k = 1, i = 1; i <= 6; i++ ) {
        for ( int j = i; j <= 6; j++ ) {
            stiffness.at(k++) = tangent.at(i, j);
        }
    }
    input.setField(stiffness, _IFT_AnisotropicLinearElasticMaterial_stiff);

    input.setField(alpha, _IFT_AnisotropicLinearElasticMaterial_talpha);
}

} // end namespace oofem
