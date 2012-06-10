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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "microplanematerial_bazant.h"
#include "microplane.h"
#include "flotarry.h"

namespace oofem {
MicroplaneMaterial_Bazant :: MicroplaneMaterial_Bazant(int n, Domain *d) : MicroplaneMaterial(n, d)
{ }



void
MicroplaneMaterial_Bazant :: giveRealStressVector(FloatArray &answer, MatResponseForm form,
                                                  GaussPoint *gp,
                                                  const FloatArray &totalStrain,
                                                  TimeStep *atTime)
{
    int i, mPlaneIndex, mPlaneIndex1;
    double SvDash, SvSum = 0.;
    double SD;
    FloatArray mPlaneNormalStress(numberOfMicroplanes), mPlaneShear_L_Stress(numberOfMicroplanes),
    mPlaneShear_M_Stress(numberOfMicroplanes);
    double mPlaneIntegrationWeight;

    Microplane *mPlane;
    FloatArray mPlaneStressCmpns, mPlaneStrainCmpns;
    FloatArray stressIncrement;

    answer.resize(6);
    answer.zero();

    StructuralMaterialStatus *status = ( StructuralMaterialStatus * ) this->giveStatus(gp);
    this->initTempStatus(gp);


    for ( mPlaneIndex = 0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++ ) {
        mPlane = this->giveMicroplane(mPlaneIndex, gp);
        mPlaneIndex1 = mPlaneIndex + 1;
        // compute strain projections on mPlaneIndex-th microplane
        computeStrainVectorComponents(mPlaneStrainCmpns, mPlane, totalStrain);
        // compute real stresses on this microplane
        giveRealMicroplaneStressVector(mPlaneStressCmpns, mPlane, mPlaneStrainCmpns, atTime);

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

        for ( i = 0; i < 6; i++ ) {
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

        for ( mPlaneIndex = 0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++ ) {
            mPlane = this->giveMicroplane(mPlaneIndex, gp);
            mPlaneIndex1 = mPlaneIndex + 1;

            updateVolumetricStressTo(mPlane, SvDash);

            SD = mPlaneNormalStress.at(mPlaneIndex1) - SvDash;
            mPlaneIntegrationWeight = this->giveMicroplaneIntegrationWeight(mPlane);

            for ( i = 0; i < 6; i++ ) {
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

    // uncomment this
    //status -> letStrainIncrementVectorBe (reducedStrainIncrement);
    status->letTempStrainVectorBe(totalStrain);

    // uncomment this
    // stressIncrement = answer;
    // crossSection->giveReducedCharacteristicVector(stressIncrement, gp, answer);
    // stressIncrement.subtract (status -> giveStressVector());
    // status -> letStressIncrementVectorBe (stressIncrement);
    status->letTempStressVectorBe(answer);
}
} // end namespace oofem
