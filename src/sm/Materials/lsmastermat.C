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

#include "lsmastermat.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "intarray.h"
#include "stressvector.h"
#include "strainvector.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LargeStrainMasterMaterial);

LargeStrainMasterMaterial :: LargeStrainMasterMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{
}

void
LargeStrainMasterMaterial :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_OPTIONAL_FIELD(ir, slaveMat, _IFT_LargeStrainMasterMaterial_slaveMat); // number of slave material
    IR_GIVE_OPTIONAL_FIELD(ir, m, _IFT_LargeStrainMasterMaterial_m); // type of Set-Hill strain tensor
}

MaterialStatus *
LargeStrainMasterMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new LargeStrainMasterMaterialStatus(gp, domain, slaveMat);
}


FloatArrayF<9>
LargeStrainMasterMaterial :: giveFirstPKStressVector_3d(const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);

    auto sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );

    //store of deformation gradient into 3x3 matrix
    auto F = from_voigt_form(vF);
    //compute right Cauchy-Green tensor(C), its eigenvalues and eigenvectors
    auto C = Tdot(F, F);
    // compute eigen values and eigen vectors of C
    FloatArray eVals;
    FloatMatrix eVecs;
    FloatMatrix(C).jaco_(eVals, eVecs, 15);

    // compute Seth - Hill's strain measure, it depends on mParameter
    double lambda1 = eVals.at(1);
    double lambda2 = eVals.at(2);
    double lambda3 = eVals.at(3);
    double E1, E2, E3;
    if ( m == 0 ) {
        E1 = 1. / 2. * log(lambda1);
        E2 = 1. / 2. * log(lambda2);
        E3 = 1. / 2. * log(lambda3);
    } else {
        E1 = 1. / ( 2. * m ) * ( pow(lambda1, m) - 1. );
        E2 = 1. / ( 2. * m ) * ( pow(lambda2, m) - 1. );
        E3 = 1. / ( 2. * m ) * ( pow(lambda3, m) - 1. );
    }

    FloatMatrixF<3,3> SethHillStrain;
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
            SethHillStrain.at(i, j) = E1 * eVecs.at(i, 1) * eVecs.at(j, 1) + E2 *eVecs.at(i, 2) * eVecs.at(j, 2) + E3 *eVecs.at(i, 3) * eVecs.at(j, 3);
        }
    }

    auto SethHillStrainVector = to_voigt_strain(SethHillStrain);
    auto stressVector = sMat->giveRealStressVector_3d(SethHillStrainVector, gp, tStep);

    auto T = this->constructTransformationMatrix(eVecs);
    stressVector.at(4) = 2;
    stressVector.at(5) = 2;
    stressVector.at(6) = 2;

    auto stressM = dot(T, FloatArrayF<6>(stressVector));
    stressM.at(4) *= 1. / 2.;
    stressM.at(5) *= 1. / 2.;
    stressM.at(6) *= 1. / 2.;

    //auto [L1, L2] = this->constructL1L2TransformationMatrices(eVals, stressM, E1, E2, E3); // c++17
    auto tmp = this->constructL1L2TransformationMatrices(eVals, stressM, E1, E2, E3);
    auto L1 = tmp.first;
    auto L2 = tmp.second;

    auto P = Tdot(T, dot(L1, T));

    //transformation of the stress to the 2PK stress and then to 1PK
    stressVector.at(4) *= 0.5;
    stressVector.at(5) *= 0.5;
    stressVector.at(6) *= 0.5;
    auto secondPK = dot(P, FloatArrayF<6>(stressVector));

    auto firstPK = dot(F, from_voigt_stress(secondPK)); // P = F*S
    auto firstPKv = to_voigt_form(firstPK);
    auto TL = Tdot(T, dot(L2, T));
    status->setPmatrix(P);
    status->setTLmatrix(TL);
    status->letTempStressVectorBe(firstPKv);
    return firstPKv;
}


FloatMatrixF<6,6>
LargeStrainMasterMaterial :: constructTransformationMatrix(const FloatMatrixF<3,3> &eVecs) const
{
    FloatMatrixF<6,6> answer;
    answer.at(1, 1) = eVecs.at(1, 1) * eVecs.at(1, 1);
    answer.at(1, 2) = eVecs.at(2, 1) * eVecs.at(2, 1);
    answer.at(1, 3) = eVecs.at(3, 1) * eVecs.at(3, 1);
    answer.at(1, 4) = eVecs.at(2, 1) * eVecs.at(3, 1);
    answer.at(1, 5) = eVecs.at(1, 1) * eVecs.at(3, 1);
    answer.at(1, 6) = eVecs.at(1, 1) * eVecs.at(2, 1);

    answer.at(2, 1) = eVecs.at(1, 2) * eVecs.at(1, 2);
    answer.at(2, 2) = eVecs.at(2, 2) * eVecs.at(2, 2);
    answer.at(2, 3) = eVecs.at(3, 2) * eVecs.at(3, 2);
    answer.at(2, 4) = eVecs.at(2, 2) * eVecs.at(3, 2);
    answer.at(2, 5) = eVecs.at(1, 2) * eVecs.at(3, 2);
    answer.at(2, 6) = eVecs.at(1, 2) * eVecs.at(2, 2);

    answer.at(3, 1) = eVecs.at(1, 3) * eVecs.at(1, 3);
    answer.at(3, 2) = eVecs.at(2, 3) * eVecs.at(2, 3);
    answer.at(3, 3) = eVecs.at(3, 3) * eVecs.at(3, 3);
    answer.at(3, 4) = eVecs.at(2, 3) * eVecs.at(3, 3);
    answer.at(3, 5) = eVecs.at(1, 3) * eVecs.at(3, 3);
    answer.at(3, 6) = eVecs.at(1, 3) * eVecs.at(2, 3);

    answer.at(4, 1) = 2 * eVecs.at(1, 2) * eVecs.at(1, 3);
    answer.at(4, 2) = 2 * eVecs.at(2, 2) * eVecs.at(2, 3);
    answer.at(4, 3) = 2 * eVecs.at(3, 2) * eVecs.at(3, 3);
    answer.at(4, 4) = eVecs.at(2, 2) * eVecs.at(3, 3) + eVecs.at(3, 2) * eVecs.at(2, 3);
    answer.at(4, 5) = eVecs.at(1, 2) * eVecs.at(3, 3) + eVecs.at(3, 2) * eVecs.at(1, 3);
    answer.at(4, 6) = eVecs.at(1, 2) * eVecs.at(2, 3) + eVecs.at(2, 2) * eVecs.at(1, 3);

    answer.at(5, 1) = 2 * eVecs.at(1, 1) * eVecs.at(1, 3);
    answer.at(5, 2) = 2 * eVecs.at(2, 1) * eVecs.at(2, 3);
    answer.at(5, 3) = 2 * eVecs.at(3, 1) * eVecs.at(3, 3);
    answer.at(5, 4) = eVecs.at(2, 1) * eVecs.at(3, 3) + eVecs.at(3, 1) * eVecs.at(2, 3);
    answer.at(5, 5) = eVecs.at(1, 1) * eVecs.at(3, 3) + eVecs.at(3, 1) * eVecs.at(1, 3);
    answer.at(5, 6) = eVecs.at(1, 1) * eVecs.at(2, 3) + eVecs.at(2, 1) * eVecs.at(1, 3);

    answer.at(6, 1) = 2 * eVecs.at(1, 1) * eVecs.at(1, 2);
    answer.at(6, 2) = 2 * eVecs.at(2, 1) * eVecs.at(2, 2);
    answer.at(6, 3) = 2 * eVecs.at(3, 1) * eVecs.at(3, 2);
    answer.at(6, 4) = eVecs.at(2, 1) * eVecs.at(3, 2) + eVecs.at(3, 1) * eVecs.at(2, 2);
    answer.at(6, 5) = eVecs.at(1, 1) * eVecs.at(3, 2) + eVecs.at(3, 1) * eVecs.at(1, 2);
    answer.at(6, 6) = eVecs.at(1, 1) * eVecs.at(2, 2) + eVecs.at(2, 1) * eVecs.at(1, 2);
    return answer;
}


std::pair<FloatMatrixF<6,6>, FloatMatrixF<6,6>>
LargeStrainMasterMaterial :: constructL1L2TransformationMatrices(const FloatArrayF<3> &eigenValues, const FloatArrayF<6> &stressM, double E1, double E2, double E3) const
{
    FloatMatrixF<6,6> answer1, answer2;
    double gamma12, gamma13, gamma23, gamma112, gamma221, gamma113, gamma331, gamma223, gamma332, gamma;
    double lambda1 = eigenValues.at(1);
    double lambda2 = eigenValues.at(2);
    double lambda3 = eigenValues.at(3);
    double lambda1P =  pow(lambda1, m - 1);
    double lambda2P =  pow(lambda2, m - 1);
    double lambda3P =  pow(lambda3, m - 1);
    if ( ( lambda1 == lambda2 ) && ( lambda2 == lambda3 ) ) {     // three equal eigenvalues
        gamma12 = gamma13 = gamma23  = 1. / 2. * lambda1P;

        answer2.at(1, 1) = 2 * stressM.at(1) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(2, 2) = 2 * stressM.at(2) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(3, 3) = 2 * stressM.at(3) * ( m - 1 ) * pow(lambda3, m - 2);
        answer2.at(4, 4) = 1. / 2. * ( stressM.at(2) + stressM.at(3) ) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(5, 5) = 1. / 2. * ( stressM.at(1) + stressM.at(3) ) * ( m - 1 ) * pow(lambda3, m - 2);
        answer2.at(6, 6) = 1. / 2. * ( stressM.at(1) + stressM.at(2) ) * ( m - 1 ) * pow(lambda1, m - 2);

        answer2.at(1, 5) = answer2.at(5, 1) = stressM.at(5) * ( m - 1 ) * pow(lambda3, m - 2);
        answer2.at(1, 6) = answer2.at(6, 1) = stressM.at(6) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(2, 4) = answer2.at(4, 2) = stressM.at(4) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(2, 6) = answer2.at(6, 2) = stressM.at(6) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(3, 4) = answer2.at(4, 3) = stressM.at(4) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(3, 5) = answer2.at(5, 3) = stressM.at(5) * ( m - 1 ) * pow(lambda3, m - 2);

        answer2.at(4, 5) = answer2.at(5, 4) =  1. / 2. * stressM.at(6) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(4, 6) = answer2.at(6, 4) = 1. / 2. * stressM.at(5) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(5, 6) = answer2.at(6, 5) = 1. / 2. * stressM.at(4) * ( m - 1 ) * pow(lambda1, m - 2);
    } else if ( lambda1 == lambda2 ) {     //two equal eigenvalues
        gamma12  = 1. / 2. * lambda1P;
        gamma13 = ( E1 - E3 ) / ( lambda1 - lambda3 );
        gamma23 = ( E2 - E3 ) / ( lambda2 - lambda3 );
        gamma113 = ( pow(lambda1, m - 1) * ( lambda1 - lambda3 ) - 2 * ( E1 - E3 ) ) / ( ( lambda1 - lambda3 ) * ( lambda1 - lambda3 ) );
        gamma331 = ( pow(lambda3, m - 1) * ( lambda3 - lambda1 ) - 2 * ( E3 - E1 ) ) / ( ( lambda3 - lambda1 ) * ( lambda3 - lambda1 ) );


        answer2.at(1, 1) = 2 * stressM.at(1) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(2, 2) = 2 * stressM.at(2) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(3, 3) = 2 * stressM.at(3) * ( m - 1 ) * pow(lambda3, m - 2);

        answer2.at(4, 4) = stressM.at(2) * gamma113 + stressM.at(3) * gamma331;
        answer2.at(5, 5) = stressM.at(1) * gamma113 + stressM.at(3) * gamma331;
        answer2.at(6, 6) = 1. / 2. * ( stressM.at(1) + stressM.at(2) ) * ( m - 1 ) * pow(lambda1, m - 2);

        answer2.at(1, 5) = answer2.at(5, 1) = 2. * stressM.at(5) * gamma113;
        answer2.at(1, 6) = answer2.at(6, 1) = stressM.at(6) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(2, 4) = answer2.at(4, 2) = 2. * stressM.at(4) * gamma113;
        answer2.at(2, 6) = answer2.at(6, 2) = stressM.at(6) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(3, 4) = answer2.at(4, 3) = 2. * stressM.at(4) * gamma331;
        answer2.at(3, 5) = answer2.at(5, 3) = 2. * stressM.at(5) * gamma331;
        answer2.at(4, 5) = answer2.at(5, 4) = stressM.at(6) * gamma113;
        answer2.at(4, 6) = answer2.at(6, 4) = stressM.at(5) * gamma113;
        answer2.at(5, 6) = answer2.at(6, 5) = stressM.at(4) * gamma113;
    } else if ( lambda2 == lambda3 ) {
        gamma23  = 1. / 2. * lambda2P;
        gamma12 = ( E1 - E2 ) / ( lambda1 - lambda2 );
        gamma13 = ( E1 - E3 ) / ( lambda1 - lambda3 );
        gamma112 = ( pow(lambda1, m - 1) * ( lambda1 - lambda2 ) - 2 * ( E1 - E2 ) ) / ( ( lambda1 - lambda2 ) * ( lambda1 - lambda2 ) );
        gamma221 = ( pow(lambda2, m - 1) * ( lambda2 - lambda1 ) - 2 * ( E2 - E1 ) ) / ( ( lambda2 - lambda1 ) * ( lambda2 - lambda1 ) );

        answer2.at(1, 1) = 2 * stressM.at(1) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(2, 2) = 2 * stressM.at(2) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(3, 3) = 2 * stressM.at(3) * ( m - 1 ) * pow(lambda3, m - 2);

        answer2.at(4, 4) = 1. / 2. * ( stressM.at(2) + stressM.at(3) ) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(5, 5) = stressM.at(1) * gamma112 + stressM.at(3) * gamma221;
        answer2.at(6, 6) = stressM.at(1) * gamma112 + stressM.at(2) * gamma221;

        answer2.at(1, 5) = answer2.at(5, 1) = 2. *  stressM.at(5) * gamma112;
        answer2.at(1, 6) = answer2.at(6, 1) = 2. *  stressM.at(6) * gamma112;
        answer2.at(2, 4) = answer2.at(4, 2) = stressM.at(4) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(2, 6) = answer2.at(6, 2) = 2. * stressM.at(6) * gamma221;
        answer2.at(3, 4) = answer2.at(4, 3) = stressM.at(4) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(3, 5) = answer2.at(5, 3) = 2. * stressM.at(5) * gamma221;
        answer2.at(4, 5) = answer2.at(5, 4) = stressM.at(6) * gamma221;
        answer2.at(4, 6) = answer2.at(6, 4) = stressM.at(5) * gamma221;
        answer2.at(5, 6) = answer2.at(6, 5) = stressM.at(4) * gamma221;
    } else if ( lambda1 == lambda3 ) {
        gamma13 = 1. / 2. * lambda1P;
        gamma12 = ( E1 - E2 ) / ( lambda1 - lambda2 );
        gamma23 = ( E2 - E3 ) / ( lambda2 - lambda3 );
        gamma223 = ( pow(lambda2, m - 1) * ( lambda2 - lambda3 ) - 2 * ( E2 - E3 ) ) / ( ( lambda2 - lambda3 ) * ( lambda2 - lambda3 ) );
        gamma332 = ( pow(lambda3, m - 1) * ( lambda3 - lambda2 ) - 2 * ( E3 - E2 ) ) / ( ( lambda3 - lambda2 ) * ( lambda3 - lambda2 ) );

        answer2.at(1, 1) = 2 * stressM.at(1) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(2, 2) = 2 * stressM.at(2) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(3, 3) = 2 * stressM.at(3) * ( m - 1 ) * pow(lambda3, m - 2);

        answer2.at(4, 4) = stressM.at(2) * gamma223 + stressM.at(3) * gamma332;
        answer2.at(5, 5) = 1. / 2. * ( stressM.at(1) + stressM.at(3) ) * ( m - 1 ) * pow(lambda3, m - 2);
        answer2.at(6, 6) = stressM.at(1) * gamma332 + stressM.at(2) * gamma223;

        answer2.at(1, 5) = answer2.at(5, 1) = stressM.at(5) * ( m - 1 ) * pow(lambda3, m - 2);
        answer2.at(1, 6) = answer2.at(6, 1) = 2. * stressM.at(6) * gamma332;
        answer2.at(2, 4) = answer2.at(4, 2) = 2. * stressM.at(4) * gamma223;
        answer2.at(2, 6) = answer2.at(6, 2) = 2. * stressM.at(4) * gamma223;
        answer2.at(3, 4) = answer2.at(4, 3) = 2. * stressM.at(4) * gamma332;
        answer2.at(3, 5) = answer2.at(5, 3) = stressM.at(5) * ( m - 1 ) * pow(lambda3, m - 2);
        answer2.at(4, 5) = answer2.at(5, 4) = stressM.at(6) * gamma332;
        answer2.at(4, 6) = answer2.at(6, 4) = stressM.at(5) * gamma332;
        answer2.at(5, 6) = answer2.at(6, 5) = stressM.at(4) * gamma332;
    } else {             //three different eigenvalues
        gamma12 = ( E1 - E2 ) / ( lambda1 - lambda2 );
        gamma13 = ( E1 - E3 ) / ( lambda1 - lambda3 );
        gamma23 = ( E2 - E3 ) / ( lambda2 - lambda3 );

        gamma112 = ( pow(lambda1, m - 1) * ( lambda1 - lambda2 ) - 2 * ( E1 - E2 ) ) / ( ( lambda1 - lambda2 ) * ( lambda1 - lambda2 ) );
        gamma221 = ( pow(lambda2, m - 1) * ( lambda2 - lambda1 ) - 2 * ( E2 - E1 ) ) / ( ( lambda2 - lambda1 ) * ( lambda2 - lambda1 ) );
        gamma113 = ( pow(lambda1, m - 1) * ( lambda1 - lambda3 ) - 2 * ( E1 - E3 ) ) / ( ( lambda1 - lambda3 ) * ( lambda1 - lambda3 ) );
        gamma331 = ( pow(lambda3, m - 1) * ( lambda3 - lambda1 ) - 2 * ( E3 - E1 ) ) / ( ( lambda3 - lambda1 ) * ( lambda3 - lambda1 ) );
        gamma223 = ( pow(lambda2, m - 1) * ( lambda2 - lambda3 ) - 2 * ( E2 - E3 ) ) / ( ( lambda2 - lambda3 ) * ( lambda2 - lambda3 ) );
        gamma332 = ( pow(lambda3, m - 1) * ( lambda3 - lambda2 ) - 2 * ( E3 - E2 ) ) / ( ( lambda3 - lambda2 ) * ( lambda3 - lambda2 ) );

        gamma = ( lambda1 * ( E2 - E3 ) + lambda2 * ( E3 - E1 ) + lambda3 * ( E1 - E2 ) ) / ( ( lambda1 - lambda2 ) * ( lambda2 - lambda3 ) * ( lambda3 - lambda1 ) );

        answer2.at(1, 1) = 2 * stressM.at(1) * ( m - 1 ) * pow(lambda1, m - 2);
        answer2.at(2, 2) = 2 * stressM.at(2) * ( m - 1 ) * pow(lambda2, m - 2);
        answer2.at(3, 3) = 2 * stressM.at(3) * ( m - 1 ) * pow(lambda3, m - 2);

        answer2.at(4, 4) = stressM.at(2) * gamma223 + stressM.at(3) * gamma332;
        answer2.at(5, 5) = stressM.at(1) * gamma113 + stressM.at(3) * gamma331;
        answer2.at(6, 6) = stressM.at(1) * gamma112 + stressM.at(2) * gamma221;

        answer2.at(1, 5) = answer2.at(5, 1) = 2. * stressM.at(5) * gamma113;
        answer2.at(1, 6) = answer2.at(6, 1) = 2. * stressM.at(6) * gamma112;
        answer2.at(2, 4) = answer2.at(4, 2) = 2. * stressM.at(4) * gamma223;
        answer2.at(2, 6) = answer2.at(6, 2) = 2. * stressM.at(6) * gamma221;
        answer2.at(3, 4) = answer2.at(4, 3) = 2. * stressM.at(4) * gamma332;
        answer2.at(3, 5) = answer2.at(5, 3) = 2. * stressM.at(5) * gamma331;
        answer2.at(4, 5) = answer2.at(5, 4) = 2. * stressM.at(6) * gamma;
        answer2.at(4, 6) = answer2.at(6, 4) = 2. * stressM.at(5) * gamma;
        answer2.at(5, 6) = answer2.at(6, 5) = 2. * stressM.at(4) * gamma;
    }


    answer1.at(1, 1) = lambda1P;
    answer1.at(2, 2) = lambda2P;
    answer1.at(3, 3) = lambda3P;
    answer1.at(4, 4) = gamma23;
    answer1.at(5, 5) = gamma13;
    answer1.at(6, 6) = gamma12;
    return {answer1, answer2};
}

FloatMatrixF<9,9>
LargeStrainMasterMaterial :: give3dMaterialStiffnessMatrix_dPdF(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );
    auto sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );

    auto stiffness = sMat->give3dMaterialStiffnessMatrix(mode, gp, tStep);
    ///////////////////////////////////////////////////////////
    stiffness.at(1, 4) = 2. * stiffness.at(1, 4);
    stiffness.at(4, 1) = 2. * stiffness.at(4, 1);
    stiffness.at(1, 5) = 2. * stiffness.at(1, 5);
    stiffness.at(5, 1) = 2. * stiffness.at(5, 1);
    stiffness.at(1, 6) = 2. * stiffness.at(1, 6);
    stiffness.at(6, 1) = 2. * stiffness.at(6, 1);
    stiffness.at(2, 4) = 2. * stiffness.at(2, 4);
    stiffness.at(4, 2) = 2. * stiffness.at(4, 2);
    stiffness.at(2, 5) = 2. * stiffness.at(2, 5);
    stiffness.at(5, 2) = 2. * stiffness.at(5, 2);
    stiffness.at(2, 6) = 2. * stiffness.at(2, 6);
    stiffness.at(6, 2) = 2. * stiffness.at(6, 2);
    stiffness.at(3, 4) = 2. * stiffness.at(3, 4);
    stiffness.at(4, 3) = 2. * stiffness.at(4, 3);
    stiffness.at(3, 5) = 2. * stiffness.at(3, 5);
    stiffness.at(5, 3) = 2. * stiffness.at(5, 3);
    stiffness.at(3, 6) = 2. * stiffness.at(3, 6);
    stiffness.at(6, 3) = 2. * stiffness.at(6, 3);
    stiffness.at(4, 4) = 4. * stiffness.at(4, 4);
    stiffness.at(4, 5) = 4. * stiffness.at(4, 5);
    stiffness.at(5, 4) = 4. * stiffness.at(5, 4);
    stiffness.at(4, 6) = 4. * stiffness.at(4, 6);
    stiffness.at(6, 4) = 4. * stiffness.at(6, 4);
    stiffness.at(5, 5) = 4. * stiffness.at(5, 5);
    stiffness.at(5, 6) = 4. * stiffness.at(5, 6);
    stiffness.at(6, 5) = 4. * stiffness.at(6, 5);
    stiffness.at(6, 6) = 4. * stiffness.at(6, 6);
    /////////////////////////////////////////////////////////////
    auto junk = dot(FloatMatrixF<6,6>(stiffness), status->givePmatrix());
    auto stiffness2 = dot(status->givePmatrix(), junk) + status->giveTLmatrix();

    return convert_dSdE_2_dPdF_3D(stiffness2, status->giveTempStressVector(), status->giveTempFVector());
}


int
LargeStrainMasterMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    auto status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );

    if ( type == IST_StressTensor ) {
        answer = status->giveStressVector();
        return 1;
    } else if ( type == IST_StrainTensor ) {
        answer = status->giveStrainVector();
        return 1;
    } else {
        StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
        return sMat->giveIPValue(answer, gp, type, tStep);
    }
}


//=============================================================================

LargeStrainMasterMaterialStatus :: LargeStrainMasterMaterialStatus(GaussPoint *g, Domain *d, int s) :
    StructuralMaterialStatus(g),
    domain(d),
    slaveMat(s)
{
    Pmatrix = eye<6>();
}


void
LargeStrainMasterMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    auto sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    auto mS = sMat->giveStatus(gp);

    mS->printOutputAt(file, tStep);
    //  StructuralMaterialStatus :: printOutputAt(file, tStep);
}


// initializes temporary variables based on their values at the previous equlibrium state
void LargeStrainMasterMaterialStatus :: initTempStatus()
{
    auto sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    auto mS = sMat->giveStatus(gp);
    mS->initTempStatus();
    //StructuralMaterialStatus :: initTempStatus();
}


// updates internal variables when equilibrium is reached
void
LargeStrainMasterMaterialStatus :: updateYourself(TimeStep *tStep)
{
    auto sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    auto mS = sMat->giveStatus(gp);
    mS->updateYourself(tStep);
    //  StructuralMaterialStatus :: updateYourself(tStep);
}


void
LargeStrainMasterMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    auto sMat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    auto mS = sMat->giveStatus(gp);
    // save parent class status
    mS->saveContext(stream, mode);
}


void
LargeStrainMasterMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: restoreContext(stream, mode);
}
} // end namespace oofem
