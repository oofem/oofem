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
#include "classfactory.h"

namespace oofem {
REGISTER_Material(HyperElasticMaterial);

HyperElasticMaterial :: HyperElasticMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{ }


void
HyperElasticMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode, GaussPoint *gp, TimeStep *tStep)

// returns the 6x6 tangent stiffness matrix - dS/dE

{
    double J2, c11, c22, c33, c12, c13, c23, A, B;
    FloatMatrix C(3, 3);
    FloatMatrix invC(3, 3);

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    C.at(1, 1) = 1. + 2. * status->giveTempStrainVector().at(1);
    C.at(2, 2) = 1. + 2. * status->giveTempStrainVector().at(2);
    C.at(3, 3) = 1. + 2. * status->giveTempStrainVector().at(3);
    C.at(2, 3) = C.at(3, 2) = status->giveTempStrainVector().at(4);
    C.at(1, 3) = C.at(3, 1) = status->giveTempStrainVector().at(5);
    C.at(1, 2) = C.at(2, 1) = status->giveTempStrainVector().at(6);

    invC.beInverseOf(C);
    J2 = C.giveDeterminant();

    c11 = invC.at(1, 1);
    c22 = invC.at(2, 2);
    c33 = invC.at(3, 3);
    c12 = invC.at(1, 2);
    c13 = invC.at(1, 3);
    c23 = invC.at(2, 3);

    A = ( K - 2. / 3. * G ) * J2;
    B = -( K - 2. / 3. * G ) * ( J2 - 1. ) + 2. * G;

    answer.resize(6, 6);

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
}


void
HyperElasticMaterial :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep)

// returns 6 components of the stress corresponding to the given total strain

{
    double J2;
    FloatMatrix C(3, 3);
    FloatMatrix invC(3, 3);
    FloatArray strainVector;

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    this->giveStressDependentPartOfStrainVector(strainVector, gp,
                                                totalStrain,
                                                tStep, VM_Total);

    C.at(1, 1) = 1. + 2. * strainVector.at(1);
    C.at(2, 2) = 1. + 2. * strainVector.at(2);
    C.at(3, 3) = 1. + 2. * strainVector.at(3);
    C.at(1, 2) = C.at(2, 1) = strainVector.at(6);
    C.at(1, 3) = C.at(3, 1) = strainVector.at(5);
    C.at(2, 3) = C.at(3, 2) = strainVector.at(4);
    invC.beInverseOf(C);
    J2 = C.giveDeterminant();

    answer.resize(6);
    double aux = ( K - 2. / 3. * G ) * ( J2 - 1. ) / 2. - G;
    answer.at(1) = aux * invC.at(1, 1) + G;
    answer.at(2) = aux * invC.at(2, 2) + G;
    answer.at(3) = aux * invC.at(3, 3) + G;
    answer.at(4) = aux * invC.at(2, 3);
    answer.at(5) = aux * invC.at(1, 3);
    answer.at(6) = aux * invC.at(1, 2);

    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}


MaterialStatus *
HyperElasticMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}


IRResultType
HyperElasticMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, K, _IFT_HyperElasticMaterial_k);
    IR_GIVE_FIELD(ir, G, _IFT_HyperElasticMaterial_g);

    return StructuralMaterial :: initializeFrom(ir);
}

} // end namespace oofem
