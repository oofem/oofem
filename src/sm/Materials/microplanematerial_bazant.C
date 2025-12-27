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

#include "microplanematerial_bazant.h"
#include "floatarray.h"

namespace oofem {
MicroplaneMaterial_Bazant :: MicroplaneMaterial_Bazant(int n, Domain *d) : MicroplaneMaterial(n, d)
{ }


FloatArrayF<6>
MicroplaneMaterial_Bazant :: giveRealStressVector_3d(const FloatArrayF<6> &strain,
                                                     GaussPoint *gp, TimeStep *tStep) const
{
    double SvDash = 0., SvSum = 0.;
    FloatArray mPlaneNormalStress(numberOfMicroplanes), mPlaneShear_L_Stress(numberOfMicroplanes),
    mPlaneShear_M_Stress(numberOfMicroplanes);

    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);

    FloatArrayF<6> answer;

    for ( int mPlaneIndex = 0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++ ) {
        int mPlaneIndex1 = mPlaneIndex + 1;
        // compute strain projections on mPlaneIndex-th microplane
        auto mPlaneStrainCmpns = computeStrainVectorComponents(mPlaneIndex1, strain);
        // compute real stresses on this microplane
        auto mPlaneStressCmpns = giveRealMicroplaneStressVector(gp, mPlaneIndex1, mPlaneStrainCmpns, tStep);

        mPlaneNormalStress.at(mPlaneIndex1) = mPlaneStressCmpns.n;
        mPlaneShear_L_Stress.at(mPlaneIndex1) = mPlaneStressCmpns.l;
        mPlaneShear_M_Stress.at(mPlaneIndex1) = mPlaneStressCmpns.m;
        double mPlaneIntegrationWeight = this->giveMicroplaneIntegrationWeight(mPlaneIndex1);

        SvSum += mPlaneNormalStress.at(mPlaneIndex1) * mPlaneIntegrationWeight;
        double SD = mPlaneNormalStress.at(mPlaneIndex1) - mPlaneStressCmpns.v;
        //SDSum +=  SD* mPlaneIntegrationWeight;

        SvDash = mPlaneStressCmpns.v;
        //volumetric stress is the same for all  mplanes
        //and does not need to be homogenized .
        //Only updating accordinging to mean normal stress must be done.
        //Use  updateVolumetricStressTo() if necessary

        // perform homogenization
        // mPlaneStressCmpns.at(1) je SVdash
        // mPlaneStressCmpns.at(2) je SN
        // mPlaneStressCmpns.at(3) je SL
        // mPlaneStressCmpns.at(4) je SM
        // answer (1 az 6)

        for ( int i = 0; i < 6; i++ ) {
            answer.at(i + 1) += ( ( N [ mPlaneIndex ] [ i ] - Kronecker [ i ] / 3. ) * SD +
                                 L [ mPlaneIndex ] [ i ] * mPlaneShear_L_Stress.at(mPlaneIndex1) +
                                 M [ mPlaneIndex ] [ i ] * mPlaneShear_M_Stress.at(mPlaneIndex1) )
                                * mPlaneIntegrationWeight;
        }
    }

    SvSum *= 6.;
    //nakonec answer take *6

    // sv=min(integr(sn)/2PI,SvDash)

    if ( SvDash > SvSum / 3. ) {
        SvDash = SvSum / 3.;
        answer = zeros<6>();

        for ( int mPlaneIndex = 0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++ ) {
            int mPlaneIndex1 = mPlaneIndex + 1;

            updateVolumetricStressTo(gp, mPlaneIndex1, SvDash);

            double SD = mPlaneNormalStress.at(mPlaneIndex1) - SvDash;
            double mPlaneIntegrationWeight = this->giveMicroplaneIntegrationWeight(mPlaneIndex1);

            for ( int i = 0; i < 6; i++ ) {
                answer.at(i + 1) += ( ( N [ mPlaneIndex ] [ i ] - Kronecker [ i ] / 3. ) * SD +
                                     L [ mPlaneIndex ] [ i ] * mPlaneShear_L_Stress.at(mPlaneIndex1) +
                                     M [ mPlaneIndex ] [ i ] * mPlaneShear_M_Stress.at(mPlaneIndex1) )
                                    * mPlaneIntegrationWeight;
            }
        }
    }

    answer *= 6.0;

    //2nd constraint, addition of volumetric part
    answer.at(1) += SvDash;
    answer.at(2) += SvDash;
    answer.at(3) += SvDash;

    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
    return answer;
}
} // end namespace oofem
