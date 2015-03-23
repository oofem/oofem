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

IRResultType
AnisotropicLinearElasticMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // read the stiffness coefficients arranged by rows from the diagonal to the right (21 values)
    FloatArray stiffness;
    IR_GIVE_FIELD(ir, stiffness, _IFT_AnisotropicLinearElasticMaterial_stiff);
    if ( stiffness.giveSize() != 21 ) {
        OOFEM_ERROR( "Incorrect size of stiff - should be 21, is %d\n", stiffness.giveSize() );
    }

    // put the stiffness coefficients into a 6x6 matrix
    int k = 1;
    for ( int i = 1; i <= 6; i++ ) {
        stiffmat.at(i, i) = stiffness.at(k++);
        for ( int j = i + 1; j <= 6; j++ ) {
            stiffmat.at(i, j) = stiffmat.at(j, i) = stiffness.at(k++);
        }
    }

    // read the thermal expansion coefficients (3 values)
    IR_GIVE_FIELD(ir, alpha, _IFT_AnisotropicLinearElasticMaterial_talpha);
    if ( alpha.giveSize() == 0 ) {
        alpha.resize(6);
        alpha.zero();
    } else if ( alpha.giveSize() != 6 )     {
        OOFEM_ERROR( "Incorrect size of talpha - should be 0 or 6, is %d\n", alpha.giveSize() );
    }

    return LinearElasticMaterial :: initializeFrom(ir);
}


void
AnisotropicLinearElasticMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    Material :: giveInputRecord(input);
    FloatArray stiffness(21);
    int k = 1;
    for ( int i = 1; i <= 6; i++ ) {
        for ( int j = i; j <= 6; j++ ) {
            stiffness.at(k++) = stiffmat.at(i, j);
        }
    }
    input.setField(stiffness, _IFT_AnisotropicLinearElasticMaterial_stiff);

    input.setField(alpha, _IFT_AnisotropicLinearElasticMaterial_talpha);
}


void
AnisotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                                  MatResponseMode mode,
                                                                  GaussPoint *gp,
                                                                  TimeStep *tStep)
{
    answer = stiffmat;
}


void
AnisotropicLinearElasticMaterial :: giveThermalDilatationVector(FloatArray &answer,
                                                                GaussPoint *gp, TimeStep *tStep)
{
    answer = alpha;
}


MaterialStatus *
AnisotropicLinearElasticMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 */
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}
} // end namespace oofem
