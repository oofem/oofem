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

#include "microplanematerial_bazant.h"
#include "microplane.h"
#include "floatarray.h"

namespace oofem {
MicroplaneMaterial_Bazant :: MicroplaneMaterial_Bazant(int n, Domain *d) : MicroplaneMaterial(n, d)
{ }


void
MicroplaneMaterial_Bazant :: giveRealStressVector_3d(FloatArray &answer,
                                                  GaussPoint *gp,
                                                  const FloatArray &totalStrain,
                                                  TimeStep *tStep)
{
    double SvDash, SvSum = 0.;
    double SD;
    FloatArray mPlaneNormalStress(numberOfMicroplanes), mPlaneShear_L_Stress(numberOfMicroplanes),
    mPlaneShear_M_Stress(numberOfMicroplanes);
    double mPlaneIntegrationWeight;

    FloatArray mPlaneStressCmpns, mPlaneStrainCmpns;

    answer.resize(6);
    answer.zero();

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);


    for ( int mPlaneIndex = 0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++ ) {
        Microplane *mPlane = this->giveMicroplane(mPlaneIndex, gp);
        int mPlaneIndex1 = mPlaneIndex + 1;
        // compute strain projections on mPlaneIndex-th microplane
        computeStrainVectorComponents(mPlaneStrainCmpns, mPlane, totalStrain);
        // compute real stresses on this microplane
        giveRealMicroplaneStressVector(mPlaneStressCmpns, mPlane, mPlaneStrainCmpns, tStep);

        mPlaneNormalStress.at(mPlaneIndex1) = mPlaneStressCmpns.at(2);
        mPlaneShear_L_Stress.at(mPlaneIndex1) = mPlaneStressCmpns.at(3);
        mPlaneShear_M_Stress.at(mPlaneIndex1) = mPlaneStressCmpns.at(4);
        mPlaneIntegrationWeight = this->giveMicroplaneIntegrationWeight(mPlane);

        SvSum += mPlaneNormalStress.at(mPlaneIndex1) * mPlaneIntegrationWeight;
        SD = mPlaneNormalStress.at(mPlaneIndex1) - mPlaneStressCmpns.at(1);
        //SDSum +=  SD* mPlaneIntegrationWeight;


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

    SvSum = SvSum * 6.;
    //nakonec answer take *6

    SvDash = mPlaneStressCmpns.at(1);
    //volumetric stress is the same for all  mplanes
    //and does not need to be homogenized .
    //Only updating accordinging to mean normal stress must be done.
    //Use  updateVolumetricStressTo() if necessary

    // sv=min(integr(sn)/2PI,SvDash)

    if ( SvDash > SvSum / 3. ) {
        SvDash = SvSum / 3.;
        answer.zero();

        for ( int mPlaneIndex = 0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++ ) {
            Microplane *mPlane = this->giveMicroplane(mPlaneIndex, gp);
            int mPlaneIndex1 = mPlaneIndex + 1;

            updateVolumetricStressTo(mPlane, SvDash);

            SD = mPlaneNormalStress.at(mPlaneIndex1) - SvDash;
            mPlaneIntegrationWeight = this->giveMicroplaneIntegrationWeight(mPlane);

            for ( int i = 0; i < 6; i++ ) {
                answer.at(i + 1) += ( ( N [ mPlaneIndex ] [ i ] - Kronecker [ i ] / 3. ) * SD +
                                     L [ mPlaneIndex ] [ i ] * mPlaneShear_L_Stress.at(mPlaneIndex1) +
                                     M [ mPlaneIndex ] [ i ] * mPlaneShear_M_Stress.at(mPlaneIndex1) )
                                    * mPlaneIntegrationWeight;
            }
        }
    }

    answer.times(6.0);

    //2nd constraint, addition of volumetric part
    answer.at(1) += SvDash;
    answer.at(2) += SvDash;
    answer.at(3) += SvDash;

    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}
} // end namespace oofem
