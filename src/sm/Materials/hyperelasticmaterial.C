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

#include "hyperelasticmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(HyperElasticMaterial);

HyperElasticMaterial :: HyperElasticMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{ }


FloatMatrixF<6,6>
HyperElasticMaterial :: give3dMaterialStiffnessMatrix(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    FloatMatrixF<3,3> C;
    C.at(1, 1) = 1. + 2. * status->giveTempStrainVector().at(1);
    C.at(2, 2) = 1. + 2. * status->giveTempStrainVector().at(2);
    C.at(3, 3) = 1. + 2. * status->giveTempStrainVector().at(3);
    C.at(2, 3) = C.at(3, 2) = status->giveTempStrainVector().at(4);
    C.at(1, 3) = C.at(3, 1) = status->giveTempStrainVector().at(5);
    C.at(1, 2) = C.at(2, 1) = status->giveTempStrainVector().at(6);

    auto invC = inv(C);
    auto J2 = det(C);

    double c11 = invC.at(1, 1);
    double c22 = invC.at(2, 2);
    double c33 = invC.at(3, 3);
    double c12 = invC.at(1, 2);
    double c13 = invC.at(1, 3);
    double c23 = invC.at(2, 3);

    double A = ( K - 2. / 3. * G ) * J2;
    double B = -( K - 2. / 3. * G ) * ( J2 - 1. ) + 2. * G;

    FloatMatrixF<6,6> answer;
    answer.at(1, 1) = ( A + B ) * c11 * c11;
    answer.at(2, 2) = ( A + B ) * c22 * c22;
    answer.at(3, 3) = ( A + B ) * c33 * c33;
    answer.at(4, 4) = A * c23 * c23 + B / 2. * ( c22 * c33 + c23 * c23 );
    answer.at(5, 5) = A * c13 * c13 + B / 2. * ( c11 * c33 + c13 * c13 );
    answer.at(6, 6) = A * c12 * c12 + B / 2. * ( c11 * c22 + c12 * c12 );
    answer.at(1, 2) = answer.at(2, 1) = A * c11 * c22 + B * c12 * c12;
    answer.at(1, 3) = answer.at(3, 1) = A * c11 * c33 + B * c13 * c13;
    answer.at(1, 4) = answer.at(4, 1) = A * c11 * c23 + B * c12 * c13;
    answer.at(1, 5) = answer.at(5, 1) = A * c11 * c13 + B * c11 * c13;
    answer.at(1, 6) = answer.at(6, 1) = A * c11 * c12 + B * c11 * c12;
    answer.at(2, 3) = answer.at(3, 2) = A * c22 * c33 + B * c23 * c23;
    answer.at(2, 4) = answer.at(4, 2) = A * c22 * c23 + B * c22 * c23;
    answer.at(2, 5) = answer.at(5, 2) = A * c22 * c13 + B * c12 * c23;
    answer.at(2, 6) = answer.at(6, 2) = A * c22 * c12 + B * c22 * c12;
    answer.at(3, 4) = answer.at(4, 3) = A * c33 * c23 + B * c33 * c23;
    answer.at(3, 5) = answer.at(5, 3) = A * c33 * c13 + B * c33 * c13;
    answer.at(3, 6) = answer.at(6, 3) = A * c33 * c12 + B * c13 * c23;
    answer.at(4, 5) = answer.at(5, 4) = A * c23 * c13 + B / 2. * ( c12 * c33 + c13 * c23 );
    answer.at(4, 6) = answer.at(6, 4) = A * c23 * c12 + B / 2. * ( c12 * c23 + c22 * c13 );
    answer.at(5, 6) = answer.at(6, 5) = A * c13 * c12 + B / 2. * ( c11 * c23 + c12 * c13 );
    return answer;
}


FloatArrayF<6>
HyperElasticMaterial :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    
    auto thermalStrain = computeStressIndependentStrainVector_3d(gp, tStep, VM_Total);
    auto elasticStrain = strain - thermalStrain;

    FloatMatrixF<3,3> C;
    C.at(1, 1) = 1. + 2. * elasticStrain.at(1);
    C.at(2, 2) = 1. + 2. * elasticStrain.at(2);
    C.at(3, 3) = 1. + 2. * elasticStrain.at(3);
    C.at(1, 2) = C.at(2, 1) = elasticStrain.at(6);
    C.at(1, 3) = C.at(3, 1) = elasticStrain.at(5);
    C.at(2, 3) = C.at(3, 2) = elasticStrain.at(4);
    auto invC = inv(C);
    double J2 = det(C);

    double aux = ( K - 2. / 3. * G ) * ( J2 - 1. ) / 2. - G;
    FloatArrayF<6> stress = {
        aux * invC.at(1, 1) + G,
        aux * invC.at(2, 2) + G,
        aux * invC.at(3, 3) + G,
        aux * invC.at(2, 3),
        aux * invC.at(1, 3),
        aux * invC.at(1, 2),
    };

    // update gp
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(stress);
    
    return stress;
}


MaterialStatus *
HyperElasticMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(gp);
}


void
HyperElasticMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, K, _IFT_HyperElasticMaterial_k);
    IR_GIVE_FIELD(ir, G, _IFT_HyperElasticMaterial_g);
}

} // end namespace oofem
