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

#include "mdm.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "isolinearelasticmaterial.h"
#include "mmaclosestiptransfer.h"
#include "nonlocalmaterialext.h"
#include "contextioerr.h"

namespace oofem {
#ifndef MDM_MAPPING_DEBUG

 #ifdef MDM_USE_MMAClosestIPTransfer
MMAClosestIPTransfer MDM :: mapper;
 #endif

 #ifdef MDM_USE_MMAShapeFunctProjection
MMAShapeFunctProjection MDM :: mapper;
 #endif

 #ifdef MDM_USE_MMALeastSquareProjection
MMALeastSquareProjection MDM :: mapper;
 #endif

#else
MMAShapeFunctProjection MDM :: mapperSFT;
MMALeastSquareProjection MDM :: mapperLST;

#endif

MMAClosestIPTransfer MDM :: mapper2;

MaterialStatus *
MDM :: CreateStatus(GaussPoint *gp) const
{
    if ( gp->giveClassID() == MicroplaneClass ) {
        return NULL;
    } else {
        return new MDMStatus(1, this->nsd, this->numberOfMicroplanes, MicroplaneMaterial :: giveDomain(), gp);
    }
}


int
MDM :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _3dMat ) || ( mode == _PlaneStress ) || ( mode == _PlaneStrain ) ) {
        return 1;
    }

    return 0;
}

void
MDM :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                            const FloatArray &totalStrain, TimeStep *atTime)
{
    FloatArray reducedStrain, strainPDC, stressPDC, stress;
    FloatArray tempDamageTensorEigenVals;
    FloatMatrix tempDamageTensor, tempDamageTensorEigenVec;
    MDMStatus *status = ( MDMStatus * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain,
                                                atTime, VM_Total);

    // computeDamageTensor : locaL OR nonlocal
    this->computeDamageTensor(tempDamageTensor, totalStrain, gp, atTime);
    // compute principal direction and corresponding eigenvalues
    this->computePDC(tempDamageTensor, tempDamageTensorEigenVals, tempDamageTensorEigenVec);

    // check the sign
    for ( int ii = 1; ii <= nsd; ii++ ) {
        if ( tempDamageTensorEigenVals.at(ii) < 0.0 ) {
            char buff [ 1024 ];
            sprintf( buff, "giveRealStressVector: negative eigenvalue of damage tensor detected, element %d, ip %d",
                    gp->giveElement()->giveNumber(), gp->giveNumber() );
            _error(buff);
        }
    }

    // Transform local strain
    // into principal damage coordinates (PDC)
    transformStrainToPDC(strainPDC, reducedStrain, tempDamageTensorEigenVec, gp);

    // Evaluate effective strain in PDC
    applyDamageTranformation(strainPDC, tempDamageTensorEigenVals);

    // Compute effective stress in PDC
    computeEffectiveStress(stressPDC, strainPDC, gp, atTime);

    // Evaluate true stress in PDC
    applyDamageTranformation(stressPDC, tempDamageTensorEigenVals);

    // Transform stress into global coordinates
    transformStressFromPDC(stress, stressPDC, tempDamageTensorEigenVec, gp);


    // update status
    status->setTempDamageTensorEigenVals(tempDamageTensorEigenVals);
    status->setTempDamageTensorEigenVec(tempDamageTensorEigenVec);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stress);
    status->setTempDamageTensor(tempDamageTensor);

    /*
     * // debug test
     * for (int i=1; i<=stress.giveSize(); i++) {
     * if (!finite(stress.at(i))) {
     * char buff[1024];
     * sprintf (buff, "giveRealStressVector: INF or NAN error detected, element %d, ip %d",
     *     gp->giveElement()->giveNumber(), gp->giveNumber());
     * _error (buff);
     * }
     * }
     * // end debug test
     */

    if ( form == ReducedForm ) {
        answer = stress;
    } else {
        crossSection->giveFullCharacteristicVector(answer, gp, stress);
    }
}


void
MDM :: computeDamageTensor(FloatMatrix &damageTensor, const FloatArray &totalStrain,
                           GaussPoint *gp, TimeStep *atTime)
{
    // local or nonlocal
    if ( nonlocal ) {
        FloatMatrix nonlocalContribution, nonlocalDamageTensor(nsd, nsd);
        MDMStatus *nonlocStatus, *status = ( MDMStatus * ) this->giveStatus(gp);

        this->buildNonlocalPointTable(gp);
        this->updateDomainBeforeNonlocAverage(atTime);

        // compute nonlocal strain increment first
        dynaList< localIntegrationRecord > *list = this->giveIPIntegrationList(gp); // !
        dynaList< localIntegrationRecord > :: iterator pos;

        for ( pos = list->begin(); pos != list->end(); ++pos ) {
            nonlocStatus = ( MDMStatus * ) this->giveStatus( ( * pos ).nearGp );
            nonlocStatus->giveLocalDamageTensorForAverage(nonlocalContribution);

            if ( ndc == 3 ) {
                nonlocalDamageTensor.at(1, 1) += nonlocalContribution.at(1, 1) * ( * pos ).weight;
                nonlocalDamageTensor.at(2, 2) += nonlocalContribution.at(2, 2) * ( * pos ).weight;
                nonlocalDamageTensor.at(1, 2) += nonlocalContribution.at(1, 2) * ( * pos ).weight;
            } else {
                nonlocalDamageTensor.at(1, 1) += nonlocalContribution.at(1, 1) * ( * pos ).weight;
                nonlocalDamageTensor.at(2, 2) += nonlocalContribution.at(2, 2) * ( * pos ).weight;
                nonlocalDamageTensor.at(3, 3) += nonlocalContribution.at(3, 3) * ( * pos ).weight;
                nonlocalDamageTensor.at(1, 2) += nonlocalContribution.at(1, 2) * ( * pos ).weight;
                nonlocalDamageTensor.at(1, 3) += nonlocalContribution.at(1, 3) * ( * pos ).weight;
                nonlocalDamageTensor.at(2, 3) += nonlocalContribution.at(2, 3) * ( * pos ).weight;
            }
        }

        nonlocalDamageTensor.times( 1. / status->giveIntegrationScale() );
        nonlocalDamageTensor.symmetrized();
        this->endIPNonlocalAverage(gp);  // !

        damageTensor = nonlocalDamageTensor;
    } else {
        computeLocalDamageTensor(damageTensor, totalStrain, gp, atTime);
    }
}





void
MDM :: computeLocalDamageTensor(FloatMatrix &damageTensor, const FloatArray &totalStrain,
                                GaussPoint *gp, TimeStep *atTime)
{
    int im1;
    double PsiOld, Psi;
    FloatArray damageVector(6);
    Microplane *mPlane;
    MDMStatus *status = ( MDMStatus * ) this->giveStatus(gp);
    ;

    // Loop over microplanes.
    for ( int im = 0; im < numberOfMicroplanes; im++ ) {
        mPlane = this->giveMicroplane(im, gp);
        im1 = im + 1;

        Psi = computeDamageOnPlane(gp, mPlane, totalStrain);

        PsiOld = status->giveMicroplaneDamage(im1);
        if ( PsiOld > Psi ) {
            Psi = PsiOld;
        }

        // update damage on microplane
        status->setMicroplaneTempDamage(im1, Psi);

        if ( formulation == COMPLIANCE_DAMAGE ) {
            //for (int i=1; i<=ndc; i++) DamageVector->at(i) += MP->N[im][i] * MP->W[im] * PsiActive;
            for ( int i = 1; i <= 6; i++ ) {
                damageVector.at(i) += this->N [ im ] [ i - 1 ] * this->microplaneWeights [ im ] * Psi;
            }
        } else if ( formulation == STIFFNESS_DAMAGE ) {
            //for (int i=1; i<=ndc; i++) DamageVector->at(i) += MP->N[im][i] * MP->W[im] / PsiActive;
            for ( int i = 1; i <= 6; i++ ) {
                damageVector.at(i) += this->N [ im ] [ i - 1 ] * this->microplaneWeights [ im ] / Psi;
            }
        }
        //else MicroplaneMaterial::_error ("Unknown type of formulation");
        else {
            _error("Unknown type of formulation");
        }
    }

    if ( ndc == 3 ) {
        // 2d case
        damageTensor.resize(2, 2);
        damageVector.times(2. / M_PI);

        damageTensor.at(1, 1) = damageVector.at(1);
        damageTensor.at(2, 2) = damageVector.at(2);
        damageTensor.at(1, 2) = damageTensor.at(2, 1) = damageVector.at(6);
    } else if ( ndc == 6 ) {
        // 3d case
        damageTensor.resize(3, 3);
        damageVector.times(6.0);

        damageTensor.at(1, 1) = damageVector.at(1);
        damageTensor.at(2, 2) = damageVector.at(2);
        damageTensor.at(3, 3) = damageVector.at(3);
        damageTensor.at(2, 3) = damageTensor.at(3, 2) = damageVector.at(4);
        damageTensor.at(3, 1) = damageTensor.at(1, 3) = damageVector.at(5);
        damageTensor.at(1, 2) = damageTensor.at(2, 1) = damageVector.at(6);

        //} else MicroplaneMaterial::_error ("computeDamageTensor: unknown ndc value encountered");
    } else {
        _error("computeDamageTensor: unknown ndc value encountered");
    }
}

#define LARGE_EXPONENT   50.0
#define HUGE_RELATIVE_COMPLIANCE 1.e20

double
MDM :: computeDamageOnPlane(GaussPoint *gp, Microplane *mplane, const FloatArray &strain)
{
    int i;
    double en, em, el, Ep, Efp, ParEpp;
    double Enorm = 0.0, sv = 0.0, answer = 0.0;
    double fmicroplane;
    IntArray mask;
    FloatArray fullStrain, prevStress = ( ( StructuralMaterialStatus * ) gp->giveMaterialStatus() )->giveStressVector();
    this->giveStressStrainMask( mask, FullForm, gp->giveMaterialMode() );
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) ( gp->giveElement()->giveCrossSection() );

    crossSection->giveFullCharacteristicVector(fullStrain, gp, strain);
    en = this->computeNormalStrainComponent(mplane, fullStrain);
    em = this->computeShearMStrainComponent(mplane, fullStrain);
    el = this->computeShearLStrainComponent(mplane, fullStrain);


    // request raw parameters
    this->giveRawMDMParameters(Efp, Ep, strain, gp);

    // compute trace of stressTensor
    if ( prevStress.isNotEmpty() ) {
        for ( i = 1; i <= nsd; i++ ) {
            if ( mask.at(i) ) {
                sv += prevStress.at( mask.at(i) );
            }
        }
    }

    ParEpp = Ep / ( 1. - ParMd ); // 1d sv reduction
    fmicroplane = linearElasticMaterial->give('E', gp) * ParEpp;
    // en /= (1.-ParMd*sv);
    en /= ( 1. - ParMd * sv / ( fmicroplane ) ); // suggested by P.Grassl (ParMd is unit dependent)

    if ( type_dam == 0 ) {
        if ( en < 0. ) {
            en = 0.;    // take the positive part of normal strain
        }

        Enorm = en;
    } else if ( type_dam == 1 ) {
        if ( en < 0. ) {
            en = 0.;    // take the positive part of normal strain
        }

        Enorm = sqrt(en * en + em * em + el * el);
        //} else MicroplaneMaterial::_error ("Unknown type of damage law");
    } else {
        _error("Unknown type of damage law");
    }


    // evaluate softening law
    if ( type_soft == 0 ) {
        if ( Enorm <= ParEpp ) {
            answer = 1.;
        } else {
            double aux = ( Enorm - ParEpp ) / Efp;
            if ( aux < LARGE_EXPONENT ) {
                answer = sqrt( ( Enorm / ParEpp ) * exp(aux) );
            } else {
                answer = HUGE_RELATIVE_COMPLIANCE;
            }
        }

        //} else MicroplaneMaterial::_error ("Unknown type of softening");
    } else {
        _error("Unknown type of softening");
    }

    return answer;
}

void
MDM :: computePDC(FloatMatrix &tempDamageTensor, FloatArray &tempDamageTensorEigenVals,
                  FloatMatrix &tempDamageTensorEigenVec)
{
    FloatMatrix help = tempDamageTensor;

    // resize the results
    tempDamageTensorEigenVals.resize(nsd);
    tempDamageTensorEigenVec.resize(nsd, nsd);

#if 0
    int nrot;
    help.Jacobi(& tempDamageTensorEigenVals, & tempDamageTensorEigenVec, & nrot);
#else
    help.jaco_(tempDamageTensorEigenVals, tempDamageTensorEigenVec, 10);
#endif
}


#define N(p, q) t.at(q, p)
#define E(p) fullStrain.at(p)

void
MDM :: transformStrainToPDC(FloatArray &answer, FloatArray &strain,
                            FloatMatrix &t, GaussPoint *gp)
{
    FloatArray fullStrain;

    if ( mdmMode == mdm_3d ) {
        ( ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection() )
        ->giveFullCharacteristicVector(fullStrain, gp, strain);

        answer.resize(6);
        answer.at(1) = N(1, 1) * ( N(1, 1) * E(1) + N(1, 2) * E(6) + N(1, 3) * E(5) )
                       + N(1, 2) * ( N(1, 2) * E(2) + N(1, 3) * E(4) )
                       + N(1, 3) *  N(1, 3) * E(3);
        answer.at(2) = N(2, 1) * ( N(2, 1) * E(1) + N(2, 2) * E(6) + N(2, 3) * E(5) )
                       + N(2, 2) * ( N(2, 2) * E(2) + N(2, 3) * E(4) )
                       + N(2, 3) *  N(2, 3) * E(3);
        answer.at(3) = N(3, 1) * ( N(3, 1) * E(1) + N(3, 2) * E(6) + N(3, 3) * E(5) )
                       + N(3, 2) * ( N(3, 2) * E(2) + N(3, 3) * E(4) )
                       + N(3, 3) *  N(3, 3) * E(3);
        answer.at(4) = N(2, 1) * ( N(3, 1) * E(1)  + N(3, 2) * E(6) / 2 + N(3, 3) * E(5) / 2 )
                       + N(2, 2) * ( N(3, 1) * E(6) / 2 + N(3, 2) * E(2)  + N(3, 3) * E(4) / 2 )
                       + N(2, 3) * ( N(3, 1) * E(5) / 2 + N(3, 2) * E(4) / 2 + N(3, 3) * E(3) );
        answer.at(4) *= 2;
        answer.at(5) = N(1, 1) * ( N(3, 1) * E(1)  + N(3, 2) * E(6) / 2 + N(3, 3) * E(5) / 2 )
                       + N(1, 2) * ( N(3, 1) * E(6) / 2 + N(3, 2) * E(2)  + N(3, 3) * E(4) / 2 )
                       + N(1, 3) * ( N(3, 1) * E(5) / 2 + N(3, 2) * E(4) / 2 + N(3, 3) * E(3) );
        answer.at(5) *= 2;
        answer.at(6) = N(1, 1) * ( N(2, 1) * E(1)  + N(2, 2) * E(6) / 2 + N(2, 3) * E(5) / 2 )
                       + N(1, 2) * ( N(2, 1) * E(6) / 2 + N(2, 2) * E(2)  + N(2, 3) * E(4) / 2 )
                       + N(1, 3) * ( N(2, 1) * E(5) / 2 + N(2, 2) * E(4) / 2 + N(2, 3) * E(3) );
        answer.at(6) *= 2;
    } else if ( mdmMode == mdm_2d ) {
        fullStrain = strain;

        answer.resize(3);
        answer.at(1) = N(1, 1) * ( N(1, 1) * E(1) + N(1, 2) * E(3) )
                       + N(1, 2) *  N(1, 2) * E(2);
        answer.at(2) = N(2, 1) * ( N(2, 1) * E(1) + N(2, 2) * E(3) )
                       + N(2, 2) *  N(2, 2) * E(2);
        answer.at(3) = N(1, 1) * ( N(2, 1) * E(1)  + N(2, 2) * E(3) / 2 )
                       + N(1, 2) * ( N(2, 1) * E(3) / 2 + N(2, 2) * E(2) );
        answer.at(3) *= 2;
    }
}

void
MDM :: applyDamageTranformation(FloatArray &strainPDC, const FloatArray &tempDamageTensorEigenVals)
{
    if ( mdmMode == mdm_3d ) {
        double psi1 = tempDamageTensorEigenVals.at(1);
        double psi2 = tempDamageTensorEigenVals.at(2);
        double psi3 = tempDamageTensorEigenVals.at(3);

        if ( formulation == COMPLIANCE_DAMAGE ) {
            strainPDC.at(1) /= psi1;
            strainPDC.at(2) /= psi2;
            strainPDC.at(3) /= psi3;
            strainPDC.at(4) /= sqrt(psi2 * psi3);
            strainPDC.at(5) /= sqrt(psi1 * psi3);
            strainPDC.at(6) /= sqrt(psi1 * psi2);
        } else if ( formulation == STIFFNESS_DAMAGE ) {
            strainPDC.at(1) *= psi1;
            strainPDC.at(2) *= psi2;
            strainPDC.at(3) *= psi3;
            strainPDC.at(4) *= sqrt(psi2 * psi3);
            strainPDC.at(5) *= sqrt(psi1 * psi3);
            strainPDC.at(6) *= sqrt(psi1 * psi2);
        } else {
            //MicroplaneMaterial::_error ("Unknown type of formulation");
            _error("Unknown type of formulation");
        }
    } else if ( mdmMode == mdm_2d ) {
        double psi1 = tempDamageTensorEigenVals.at(1);
        double psi2 = tempDamageTensorEigenVals.at(2);
        if ( formulation == COMPLIANCE_DAMAGE ) {
            strainPDC.at(1) /= psi1;
            strainPDC.at(2) /= psi2;
            strainPDC.at(3) /= sqrt(psi1 * psi2);
        } else if ( formulation == STIFFNESS_DAMAGE ) {
            strainPDC.at(1) *= psi1;
            strainPDC.at(2) *= psi2;
            strainPDC.at(3) *= sqrt(psi1 * psi2);
        } else {
            //MicroplaneMaterial::_error ("Unknown type of formulation");
            _error("Unknown type of formulation");
        }

        //} else MicroplaneMaterial::_error ("Unknown type of mdm mode");
    } else {
        _error("Unknown type of mdm mode");
    }
}


void
MDM :: computeEffectiveStress(FloatArray &stressPDC, const FloatArray &strainPDC, GaussPoint *gp, TimeStep *atTime)
{
    FloatMatrix de;
    if ( mdmMode == mdm_3d ) {
        /// PDC components in 3d mode are in full 3d format, even in planeStrain situation
        this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(de, ReducedForm, TangentStiffness, gp, atTime);
    } else {
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(de, ReducedForm, TangentStiffness, gp, atTime);
    }

    stressPDC.beProductOf(de, strainPDC);
}



#define Nt(p, q) t.at(p, q)
#define S(p) stressPDC.at(p)

void
MDM :: transformStressFromPDC(FloatArray &answer, const FloatArray &stressPDC, const FloatMatrix &t, GaussPoint *gp)
{
    if ( mdmMode == mdm_3d ) {
        FloatArray fullAnswer(6);
        //answer.resize (6);

        fullAnswer.at(1) = Nt(1, 1) * ( Nt(1, 1) * S(1) + 2 * Nt(1, 2) * S(6) + 2 * Nt(1, 3) * S(5) )
                           + Nt(1, 2) * ( Nt(1, 2) * S(2) + 2 * Nt(1, 3) * S(4) )
                           + Nt(1, 3) *  Nt(1, 3) * S(3);
        fullAnswer.at(2) = Nt(2, 1) * ( Nt(2, 1) * S(1) + 2 * Nt(2, 2) * S(6) + 2 * Nt(2, 3) * S(5) )
                           + Nt(2, 2) * ( Nt(2, 2) * S(2) + 2 * Nt(2, 3) * S(4) )
                           + Nt(2, 3) *  Nt(2, 3) * S(3);
        fullAnswer.at(3) = Nt(3, 1) * ( Nt(3, 1) * S(1) + 2 * Nt(3, 2) * S(6) + 2 * Nt(3, 3) * S(5) )
                           + Nt(3, 2) * ( Nt(3, 2) * S(2) + 2 * Nt(3, 3) * S(4) )
                           + Nt(3, 3) *  Nt(3, 3) * S(3);
        fullAnswer.at(4) = Nt(2, 1) * ( Nt(3, 1) * S(1)  + Nt(3, 2) * S(6)  + Nt(3, 3) * S(5) )
                           + Nt(2, 2) * ( Nt(3, 1) * S(6)  + Nt(3, 2) * S(2)  + Nt(3, 3) * S(4) )
                           + Nt(2, 3) * ( Nt(3, 1) * S(5)  + Nt(3, 2) * S(4)  + Nt(3, 3) * S(3) );
        fullAnswer.at(5) = Nt(1, 1) * ( Nt(3, 1) * S(1)  + Nt(3, 2) * S(6)  + Nt(3, 3) * S(5) )
                           + Nt(1, 2) * ( Nt(3, 1) * S(6)  + Nt(3, 2) * S(2)  + Nt(3, 3) * S(4) )
                           + Nt(1, 3) * ( Nt(3, 1) * S(5)  + Nt(3, 2) * S(4)  + Nt(3, 3) * S(3) );
        fullAnswer.at(6) = Nt(1, 1) * ( Nt(2, 1) * S(1)  + Nt(2, 2) * S(6)  + Nt(2, 3) * S(5) )
                           + Nt(1, 2) * ( Nt(2, 1) * S(6)  + Nt(2, 2) * S(2)  + Nt(2, 3) * S(4) )
                           + Nt(1, 3) * ( Nt(2, 1) * S(5)  + Nt(2, 2) * S(4)  + Nt(2, 3) * S(3) );

        ( ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection() )
        ->giveReducedCharacteristicVector(answer, gp, fullAnswer);
    } else if ( mdmMode == mdm_2d ) {
        answer.resize(3);

        answer.at(1) = Nt(1, 1) * ( Nt(1, 1) * S(1) + Nt(1, 2) * 2. * S(3) )
                       + Nt(1, 2) *  Nt(1, 2) * S(2);
        answer.at(2) = Nt(2, 1) * ( Nt(2, 1) * S(1) + Nt(2, 2) * 2. * S(3) )
                       + Nt(2, 2) *  Nt(2, 2) * S(2);
        answer.at(3) = Nt(1, 1) * ( Nt(2, 1) * S(1) + Nt(2, 2) * S(3) )
                       + Nt(1, 2) * ( Nt(2, 1) * S(3) + Nt(2, 2) * S(2) );
    }
}



void
MDM :: giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                   MatResponseForm form,
                                   MatResponseMode mode,
                                   GaussPoint *gp,
                                   TimeStep *atTime)
{
    FloatMatrix de;
    MDMStatus *status = ( MDMStatus * ) this->giveStatus(gp);

    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(de, ReducedForm, TangentStiffness, gp, atTime);
    //answer = de;
    //return;
    // if (isVirgin()) return ;
    if ( ( mode == TangentStiffness ) || ( mode == SecantStiffness ) ) {
        // Apply damage in principal coordinates
        this->applyDamageToStiffness(de, gp);

        // Transform to global coordinates (in reduced space)
        this->transformStiffnessfromPDC( de, * ( status->giveTempDamageTensorEigenVec() ) );
    }


    if ( form == ReducedForm ) {
        answer = de;
    } else {
        IntArray mask;
        this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode() );
        answer.beSubMatrixOfSizeOf(de, mask, 6);
    }
}


#define MAX_REL_COMPL_TRESHOLD 1.e6

void
MDM :: applyDamageToStiffness(FloatMatrix &d, GaussPoint *gp)
{
    MDMStatus *status = ( MDMStatus * ) this->giveStatus(gp);

    if ( mdmMode == mdm_3d ) {
        double psi1 = 0.0, psi2 = 0.0, psi3 = 0.0;
        if ( formulation == COMPLIANCE_DAMAGE ) {
            psi1 = status->giveTempDamageTensorEigenVals()->at(1);
            psi2 = status->giveTempDamageTensorEigenVals()->at(2);
            psi3 = status->giveTempDamageTensorEigenVals()->at(3);
        } else if ( formulation == STIFFNESS_DAMAGE ) {
            psi1 = 1. / ( status->giveTempDamageTensorEigenVals()->at(1) );
            psi2 = 1. / ( status->giveTempDamageTensorEigenVals()->at(2) );
            psi3 = 1. / ( status->giveTempDamageTensorEigenVals()->at(3) );
            //} else MicroplaneMaterial::_error ("Unknown type of formulation");
        } else {
            _error("Unknown type of formulation");
        }

        //if ((psi1 > MAX_REL_COMPL_TRESHOLD) || (psi2 > MAX_REL_COMPL_TRESHOLD) || (psi3 > MAX_REL_COMPL_TRESHOLD)) printf (":");

        psi1 = min(psi1, MAX_REL_COMPL_TRESHOLD);
        psi2 = min(psi2, MAX_REL_COMPL_TRESHOLD);
        psi3 = min(psi3, MAX_REL_COMPL_TRESHOLD);

        if ( d.giveNumberOfRows() == 6 ) {
            d.at(1, 1) /= ( psi1 * psi1 );
            d.at(1, 2) /= ( psi1 * psi2 );
            d.at(1, 3) /= ( psi1 * psi3 );
            d.at(2, 1) /= ( psi2 * psi1 );
            d.at(2, 2) /= ( psi2 * psi2 );
            d.at(2, 3) /= ( psi2 * psi3 );
            d.at(3, 1) /= ( psi3 * psi1 );
            d.at(3, 2) /= ( psi3 * psi2 );
            d.at(3, 3) /= ( psi3 * psi3 );
            d.at(4, 4) /= ( psi2 * psi3 );
            d.at(5, 5) /= ( psi1 * psi3 );
            d.at(6, 6) /= ( psi1 * psi2 );
        } else if ( d.giveNumberOfRows() == 4 ) {
            d.at(1, 1) /= ( psi1 * psi1 );
            d.at(1, 2) /= ( psi1 * psi2 );
            d.at(1, 3) /= ( psi1 * psi3 );
            d.at(2, 1) /= ( psi2 * psi1 );
            d.at(2, 2) /= ( psi2 * psi2 );
            d.at(2, 3) /= ( psi2 * psi3 );
            d.at(3, 1) /= ( psi3 * psi1 );
            d.at(3, 2) /= ( psi3 * psi2 );
            d.at(3, 3) /= ( psi3 * psi3 );
            d.at(4, 4) /= ( psi1 * psi2 );
            //} else MicroplaneMaterial::_error ("Unknown type stiffness");
        } else {
            _error("Unknown type stiffness");
        }

        return;
    } else if ( mdmMode == mdm_2d ) {
        double psi1 = 0.0, psi2 = 0.0;
        if ( formulation == COMPLIANCE_DAMAGE ) {
            psi1 = status->giveTempDamageTensorEigenVals()->at(1);
            psi2 = status->giveTempDamageTensorEigenVals()->at(2);
        } else if ( formulation == STIFFNESS_DAMAGE ) {
            psi1 = 1. / ( status->giveTempDamageTensorEigenVals()->at(1) );
            psi2 = 1. / ( status->giveTempDamageTensorEigenVals()->at(2) );
            //} else MicroplaneMaterial::_error ("Unknown type of formulation");
        } else {
            _error("Unknown type of formulation");
        }


        //if ((psi1 > MAX_REL_COMPL_TRESHOLD) || (psi2 > MAX_REL_COMPL_TRESHOLD)) printf (":");
        psi1 = min(psi1, MAX_REL_COMPL_TRESHOLD);
        psi2 = min(psi2, MAX_REL_COMPL_TRESHOLD);

        d.at(1, 1) /= psi1 * psi1;
        d.at(1, 2) /= psi1 * psi2;
        d.at(2, 1) /= psi2 * psi1;
        d.at(2, 2) /= psi2 * psi2;
        d.at(3, 3) /= psi1 * psi2;

        return;
    }
}


void
MDM :: transformStiffnessfromPDC(FloatMatrix &de, const FloatMatrix &t)
{
    this->rotateTensor4(de, t);
}


void
MDM :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                     MatResponseForm form,
                                     MatResponseMode mode,
                                     GaussPoint *gp,
                                     TimeStep *atTime)
{
    this->giveMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
}

void
MDM :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode,
                                GaussPoint *gp, TimeStep *atTime)
{
    this->giveMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
}


void
MDM :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode,
                                GaussPoint *gp, TimeStep *atTime)
{
    this->giveMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
}



int
MDM :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    MDMStatus *status = ( MDMStatus * ) this->giveStatus(aGaussPoint);

    if ( type == IST_DamageTensor ) {
        // DamageTensor = I-Phi*Phi \approx I - Phi^{-1}*Phi*^{-1}
        FloatMatrix d, c;
        if ( formulation == COMPLIANCE_DAMAGE ) {
            status->giveDamageTensor(d); // d corresponds damage compliance
            c.beInverseOf(d);    // c corresponds to damage
        } else {
            status->giveDamageTensor(c);
        }

        d.beProductOf(c, c);
        if ( ndc == 3 ) {
            answer.resize(3);
            answer.at(1) = 1. - d.at(1, 1);
            answer.at(2) = 1. - d.at(2, 2);
            answer.at(3) = d.at(1, 2);
        } else {
            answer.resize(6);
            answer.at(1) = 1. - d.at(1, 1);
            answer.at(2) = 1. - d.at(2, 2);
            answer.at(3) = 1. - d.at(3, 3);
            answer.at(4) = d.at(2, 3);
            answer.at(5) = d.at(1, 3);
            answer.at(6) = d.at(1, 2);
        }

        return 1;
    } else if ( type == IST_DamageTensorTemp ) {
        FloatMatrix d, c;
        if ( formulation == COMPLIANCE_DAMAGE ) {
            status->giveTempDamageTensor(d);
            c.beInverseOf(d);
        } else {
            status->giveTempDamageTensor(c);
        }

        d.beProductOf(c, c);
        if ( ndc == 3 ) {
            answer.resize(3);
            answer.at(1) = 1. - d.at(1, 1);
            answer.at(2) = 1. - d.at(2, 2);
            answer.at(3) = d.at(1, 2);
        } else {
            answer.resize(6);
            answer.at(1) = 1. - d.at(1, 1);
            answer.at(2) = 1. - d.at(2, 2);
            answer.at(3) = 1. - d.at(3, 3);
            answer.at(4) = d.at(2, 3);
            answer.at(5) = d.at(1, 3);
            answer.at(6) = d.at(1, 2);
        }

        return 1;
    } else if ( type == IST_PrincipalDamageTensor ) {
        int i;
        FloatArray const *d;

        answer.resize(nsd);
        d = status->giveDamageTensorEigenVals();
        if ( formulation == COMPLIANCE_DAMAGE ) {
            for ( i = 1; i <= nsd; i++ ) {
                answer.at(i) = 1. - 1. / d->at(i) / d->at(i);
            }
        } else {
            for ( i = 1; i <= nsd; i++ ) {
                answer.at(i) = 1. - d->at(i) * d->at(i);
            }
        }

        return 1;
    } else if ( type == IST_PrincipalDamageTempTensor ) {
        int i;
        FloatArray const *d;

        answer.resize(nsd);
        d = status->giveTempDamageTensorEigenVals();
        if ( formulation == COMPLIANCE_DAMAGE ) {
            for ( i = 1; i <= nsd; i++ ) {
                answer.at(i) = 1. - 1. / d->at(i) / d->at(i);
            }
        } else {
            for ( i = 1; i <= nsd; i++ ) {
                answer.at(i) = 1. - d->at(i) * d->at(i);
            }
        }

        return 1;
    } else if ( type == IST_MicroplaneDamageValues ) {
        status->giveMicroplaneDamageValues(answer);
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}




InternalStateValueType
MDM :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_DamageTensor ) || ( type == IST_DamageTensorTemp ) ||
        ( type == IST_DamageInvTensor ) || ( type == IST_DamageInvTensorTemp ) ||
        ( type == IST_PrincipalDamageTensor ) || ( type == IST_PrincipalDamageTempTensor ) ) {
        return ISVT_TENSOR_S3;
    } else if ( type == IST_MicroplaneDamageValues ) {
        return ISVT_VECTOR;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}


int
MDM :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( ( type == IST_DamageTensor ) || ( type == IST_DamageTensorTemp ) ||
        ( type == IST_DamageTensor ) || ( type == IST_DamageTensorTemp ) ) {
        answer.resize(6);
        if ( this->ndc == 3 ) {
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(6) = 3;
        } else { // ndc == 6
            for ( int i = 1; i <= 6; i++ ) {
                answer.at(i) = i;
            }
        }

        return 1;
    } else if ( ( type == IST_PrincipalDamageTensor ) || ( type == IST_PrincipalDamageTempTensor ) ) {
        answer.resize(6);
        answer.zero();
        for ( int i = 1; i <= nsd; i++ ) {
            answer.at(i) = i;
        }

        return 1;
    } else if ( type == IST_MicroplaneDamageValues ) {
        answer.resize(numberOfMicroplanes);
        answer.zero();
        for ( int i = 1; i <= numberOfMicroplanes; i++ ) {
            answer.at(i) = i;
        }

        return 1;
    } else {
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
MDM :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( ( type == IST_DamageTensor ) || ( type == IST_DamageTensor ) ||
        ( type == IST_DamageTensor ) || ( type == IST_DamageTensor ) ) {
        return ndc;
    } else if ( ( type == IST_PrincipalDamageTensor ) || ( type == IST_PrincipalDamageTempTensor ) ) {
        return nsd;
    } else if ( type == IST_MicroplaneDamageValues ) {
        return numberOfMicroplanes;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}


void
MDM :: giveThermalDilatationVector(FloatArray &answer,
                                   GaussPoint *gp,  TimeStep *tStep)
{
    answer.resize(6);
    answer.zero();
    answer.at(1) = this->tempDillatCoeff;
    answer.at(2) = this->tempDillatCoeff;
    answer.at(3) = this->tempDillatCoeff;
}


IRResultType
MDM :: initializeFrom(InputRecord *ir)
//
// initializes according to string
//
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->tempDillatCoeff, IFT_MDM_talpha, "talpha"); // Macro
    IR_GIVE_FIELD(ir, this->ParMd, IFT_MDM_parmd, "parmd"); // Macro
    IR_GIVE_FIELD(ir, this->nonlocal, IFT_MDM_nonloc, "nonloc"); // Macro

    if ( this->nonlocal ) {
        IR_GIVE_FIELD(ir, R, IFT_MDM_r, "r"); // Macro
        if ( R < 0.0 ) {
            R = 0.0;
        }

        if ( ( ir->hasField(IFT_MDM_efp, "efp") ) && ( ir->hasField(IFT_MDM_ep, "ep") ) ) {
            // read raw_params if available
            IR_GIVE_FIELD(ir, this->mdm_Efp, IFT_MDM_efp, "efp"); // Macro
            IR_GIVE_FIELD(ir, this->mdm_Ep, IFT_MDM_ep, "ep"); // Macro
        } else if ( ( ir->hasField(IFT_MDM_gf, "gf") ) && ( ir->hasField(IFT_MDM_ft, "ft") ) ) {
            IR_GIVE_FIELD(ir, this->Gf, IFT_MDM_gf, "gf"); // Macro
            IR_GIVE_FIELD(ir, this->Ft, IFT_MDM_ft, "ft"); // Macro
        } else {
            _error("instanciateFrom: unknown set of parameters");
        }
    } else { // local case
        if ( ( ir->hasField(IFT_MDM_efp, "efp") ) && ( ir->hasField(IFT_MDM_ep, "ep") ) ) {
            // read raw_params if available
            IR_GIVE_FIELD(ir, this->mdm_Efp, IFT_MDM_efp, "efp"); // Macro
            IR_GIVE_FIELD(ir, this->mdm_Ep, IFT_MDM_ep, "ep"); // Macro
        } else if ( ( ir->hasField(IFT_MDM_gf, "gf") ) && ( ir->hasField(IFT_MDM_ep, "ep") ) ) {
            IR_GIVE_FIELD(ir, this->Gf, IFT_MDM_gf, "gf"); // Macro
            IR_GIVE_FIELD(ir, this->mdm_Ep, IFT_MDM_ep, "ep"); // Macro
        } else {
            _error("instanciateFrom: unknown set of parameters");
        }
    }

    // read formulation
    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_MDM_formulation, "formulation"); // Macro
    this->formulation = ( MDMFormulatrionType ) _val;

    IR_GIVE_FIELD(ir, _val, IFT_MDM_mode, "mode"); // Macro
    this->mdmMode     = ( MDMModeType ) _val;

    if ( this->mdmMode == mdm_3d ) {
        this->ndc = 6;
        this->nsd = 3;
    } else if ( this->mdmMode == mdm_2d ) {
        this->ndc = 3;
        this->nsd = 2;
    }


#ifdef MDM_MAPPING_DEBUG
    _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_MDM_mapper, "mapper"); // Macro
    this->mapperType = ( MDMMapperType ) _val;
    OOFEM_LOG_INFO("MDM: using optional mapper %d\n", mapperType);
#endif



    MicroplaneMaterial :: initializeFrom(ir);
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);

    linearElasticMaterial = new IsotropicLinearElasticMaterial( 1, MicroplaneMaterial :: giveDomain() );
    linearElasticMaterial->initializeFrom(ir);

#ifdef MDM_MAPPING_DEBUG
    mapperSFT.initializeFrom(ir);
    mapperLST.initializeFrom(ir);
#else
    mapper.initializeFrom(ir);
#endif
    mapper2.initializeFrom(ir);

    return IRRT_OK;
}


int
MDM :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    MicroplaneMaterial :: giveInputRecordString(str, keyword);
    this->giveLinearElasticMaterial()->giveInputRecordString(str, false);
    sprintf(buff, " talpha %e", this->tempDillatCoeff);
    str += buff;
    sprintf(buff, " parmd %e", this->ParMd);
    str += buff;
    sprintf(buff, " nonloc %d", this->nonlocal);
    str += buff;
    if ( this->nonlocal ) {
        sprintf(buff, " r %e", this->R);
        str += buff;

        if ( this->mdm_Ep >= 0.0 && this->mdm_Efp >= 0.0 ) {
            sprintf(buff, " efp %e", this->mdm_Efp);
            str += buff;
            sprintf(buff, " ep %e", this->mdm_Ep);
            str += buff;
        } else {
            sprintf(buff, " gf %e", this->Gf);
            str += buff;
            sprintf(buff, " ft %e", this->Ft);
            str += buff;
        }
    } else {             // local case
        if ( this->mdm_Ep >= 0.0 && this->mdm_Efp >= 0.0 ) {
            sprintf(buff, " efp %e", this->mdm_Efp);
            str += buff;
            sprintf(buff, " ep %e", this->mdm_Ep);
            str += buff;
        } else {
            sprintf(buff, " gf %e", this->Gf);
            str += buff;
            if ( this->mdm_Ep >= 0.0 ) {
                sprintf(buff, " ep %e", this->mdm_Ep);
                str += buff;
            }
        }
    }

    if ( this->formulation ) {
        sprintf( buff, " formulation %d", ( int ) ( this->formulation ) );
        str += buff;
    }

    sprintf( buff, " mode %d", ( int ) ( this->mdmMode ) );
    str += buff;

#ifdef MDM_MAPPING_DEBUG
    sprintf( buff, " mapper %d", ( int ) ( this->mapperType ) );
    str += buff;
#endif

    StructuralNonlocalMaterialExtensionInterface :: giveInputRecordString(str, false);

#ifdef MDM_MAPPING_DEBUG
    this->mapperSFT.giveInputRecordString(str, false);
    this->mapperLST.giveInputRecordString(str, false);
#else
    this->mapper.giveInputRecordString(str, false);
#endif
    this->mapper2.giveInputRecordString(str, false);

    return 1;
}


void
MDM :: rotateTensor4(FloatMatrix &Dlocal, const FloatMatrix &t)
//
// Interpreting "t" as the rotation matrix containing in its
// columns the local base vectors, transform a fourth-order tensor
// represented by a matrix from local to global coordinates.
{
    FloatMatrix tt, td;
    this->formTransformationMatrix( tt, t, Dlocal.giveNumberOfRows() );
    td.beProductOf(tt, Dlocal);
    Dlocal.beProductTOf(td, tt);
}

#define FORMT33(i, j, ij) answer.at(ij, 1) = t.at(i, 1) * t.at(j, 1); answer.at(ij, 2) = t.at(i, 2) * t.at(j, 2); answer.at(ij, 3) = t.at(i, 1) * t.at(j, 2) + t.at(i, 2) * t.at(j, 1)

#define FORMT44(i, j, ij) answer.at(ij, 1) = t.at(i, 1) * t.at(j, 1); answer.at(ij, 2) = t.at(i, 2) * t.at(j, 2); answer.at(ij, 3) = t.at(i, 3) * t.at(j, 3); answer.at(ij, 4) = t.at(i, 1) * t.at(j, 2) + t.at(i, 2) * t.at(j, 1)

#define FORMT66(i, j, ij) answer.at(ij, 1) = t.at(i, 1) * t.at(j, 1); answer.at(ij, 2) = t.at(i, 2) * t.at(j, 2); answer.at(ij, 3) = t.at(i, 3) * t.at(j, 3); answer.at(ij, 4) = t.at(i, 2) * t.at(j, 3) + t.at(i, 3) * t.at(j, 2); answer.at(ij, 5) = t.at(i, 1) * t.at(j, 3) + t.at(i, 3) * t.at(j, 1); answer.at(ij, 6) = t.at(i, 1) * t.at(j, 2) + t.at(i, 2) * t.at(j, 1)


void
MDM :: formTransformationMatrix(FloatMatrix &answer, const FloatMatrix &t, int n)
//
// Interpreting "t" as the rotation matrix containing in its
// columns the local base vectors, form the transformation matrix
// for the transformation of a fourth-order tensor represented by
// a matrix from local to global coordinates.
{
    switch ( n ) {
    case 1:
        answer.resize(1, 1);
        answer.at(1, 1) = 1.0;
        return;

    case 3:
        answer.resize(3, 3);
        answer.zero();

        FORMT33(1, 1, 1);
        FORMT33(2, 2, 2);
        FORMT33(1, 2, 3);
        return;

    case 4:
        answer.resize(4, 4);
        answer.zero();

        FORMT44(1, 1, 1);
        FORMT44(2, 2, 2);
        FORMT44(3, 3, 3);
        FORMT44(1, 2, 4);
        return;

    case 6:
        answer.resize(6, 6);
        answer.zero();

        FORMT66(1, 1, 1);
        FORMT66(2, 2, 2);
        FORMT66(3, 3, 3);
        FORMT66(2, 3, 4);
        FORMT66(1, 3, 5);
        FORMT66(1, 2, 6);
        return;

    default:
        //MicroplaneMaterial::_error ("Stress transformation matrix format not implemented");
        _error("Stress transformation matrix format not implemented");
    }
}

void
MDM :: giveRawMDMParameters(double &Efp, double &Ep, const FloatArray &reducedStrain, GaussPoint *gp)
{
    // test if raw parameters are given
    if ( this->mdm_Efp > 0.0 ) {
        Efp = this->mdm_Efp;
        Ep  = this->mdm_Ep;
        return;
    }

    // determine params from macroscopic ones
    if ( nonlocal ) {
        // formulas derived for 3d case
        double EModulus = linearElasticMaterial->give('E', gp);
        double gammaf = ( EModulus * this->Gf ) / ( this->R * this->Ft * this->Ft );
        double gamma  = gammaf / ( 1.47 - 0.0014 * gammaf );
        double f = this->Ft / ( 1.56 + 0.006 * gamma ); // microplane tensile strength
        Efp = this->mdm_Efp = ( gamma * f ) / EModulus;
        //double mdtilda = this->ParMd/f;
        Ep = this->mdm_Ep   = f / EModulus;
        return;
    } else { // local model
        StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();
        FloatArray strainVector, principalStrain;
        FloatMatrix dirs;
        double h;

        crossSection->giveFullCharacteristicVector(strainVector, gp, reducedStrain);
        this->computePrincipalValDir(principalStrain, dirs,
                                     strainVector,
                                     principal_strain);

        // find EigenVector With Largest EigenValue
        int indx = 1;
        double max = principalStrain.at(1);
        if ( principalStrain.at(2) > max ) {
            indx = 2;
        }

        if ( principalStrain.at(3) > max ) {
            indx = 3;
        }

        FloatArray dir(3);
        dir.at(1) = dirs.at(1, indx);
        dir.at(2) = dirs.at(2, indx);
        dir.at(3) = dirs.at(3, indx);

        h  = gp->giveElement()->giveCharacteristicLenght(gp, dir);

        double E  = this->giveLinearElasticMaterial()->give(Ex, gp);
        Ep = this->mdm_Ep;
        if ( nsd == 2 ) {
            Efp = ( Gf / ( h * E * Ep ) + 1.2 * Ep ) / 1.75 - Ep;
        } else {
            Efp = ( Gf / ( h * E * Ep ) + ( 2.13 + ParMd ) * Ep ) / ( 2.73 - ParMd ) - Ep;
        }

        if ( Efp <= 0. ) {
            //MicroplaneMaterial::_error("Warning: negative Efp encountered");
            _error("Warning: negative Efp encountered");
        }

        return;
    }
}

/*
 * double
 * MDM::giveParameterEfp(const FloatArray& reducedStrain, GaussPoint* gp)
 * {
 * if (Efp > 0.)
 *  return Efp;
 *
 *
 * StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();
 * FloatArray strainVector, principalStrain;
 * FloatMatrix dirs;
 * double h;
 *
 * if (nonlocal) h = 1.0;
 * else {
 * crossSection->giveFullCharacteristicVector(strainVector, gp, reducedStrain);
 * this->computePrincipalValDir (principalStrain, dirs,
 *               strainVector,
 *               principal_strain);
 *
 * // find EigenVector With Largest EigenValue
 * int indx = 1;
 * double max = principalStrain.at(1);
 * if (principalStrain.at(2) > max) indx = 2;
 * if (principalStrain.at(3) > max) indx = 3;
 *
 * FloatArray dir (3);
 * dir.at(1) = dirs.at(1,indx);
 * dir.at(2) = dirs.at(2,indx);
 * dir.at(3) = dirs.at(3,indx);
 *
 * h  = gp->giveElement()->giveCharacteristicLenght (gp, dir);
 * }
 *
 * double E  = this -> giveLinearElasticMaterial()->give(Ex);
 * if (nsd==2)
 * Efp = (Gf/(h*E*Ep)+1.2*Ep)/1.75 - Ep;
 * else
 * Efp = (Gf/(h*E*Ep)+(2.13+ParMd)*Ep)/(2.73-ParMd) - Ep;
 * if (Efp<=0.)
 * MicroplaneMaterial::_error("Warning: negative Efp encountered");
 *
 * return Efp;
 * }
 */

void
MDM :: initializeData(int numberOfMicroplanes)
{
    if ( this->mdmMode == mdm_3d ) {
        MicroplaneMaterial :: initializeData(numberOfMicroplanes);
    } else if ( this->mdmMode == mdm_2d ) {
        if ( numberOfMicroplanes > MAX_NUMBER_OF_MICROPLANES ) {
            //MicroplaneMaterial::_error ("initializeData: required number of microplanes too big");
            _error("initializeData: required number of microplanes too big");
        }

        int iplane;
        int i, ii, jj;
        double alpha = M_PI / numberOfMicroplanes;
        FloatArray n(3), m(3), l(3);

        int ij [ 6 ] [ 2 ] = { { 1, 1 }, { 2, 2 }, { 3, 3 }, { 2, 3 }, { 3, 1 }, { 1, 2 } };

        for ( iplane = 0; iplane < numberOfMicroplanes; iplane++ ) {
            microplaneWeights [ iplane ] = alpha;

            n.at(1) = microplaneNormals [ iplane ] [ 0 ] = cos(iplane * alpha);
            n.at(2) = microplaneNormals [ iplane ] [ 1 ] = sin(iplane * alpha);
            n.at(3) = microplaneNormals [ iplane ] [ 2 ] = 0.0;

            // compute projection tensors for each microplane
            m.at(1) = n.at(2);
            m.at(2) = -n.at(1);
            m.at(3) = 0.0;

            l.at(1) = 0.0;
            l.at(2) = 0.0;
            l.at(3) = 1.0;

            for ( i = 0; i < 6; i++ ) {
                ii = ij [ i ] [ 0 ];
                jj = ij [ i ] [ 1 ];

                N [ iplane ] [ i ] = n.at(ii) * n.at(jj);
                M [ iplane ] [ i ] = 0.5 * ( m.at(ii) * n.at(jj) + m.at(jj) * n.at(ii) );
                L [ iplane ] [ i ] = 0.5 * ( l.at(ii) * n.at(jj) + l.at(jj) * n.at(ii) );
            }
        } // end loop over mplanes

        //} else MicroplaneMaterial::_error ("initializeData: Unknown MDMModeType ecountered");
    } else {
        _error("initializeData: Unknown MDMModeType ecountered");
    }
}


void
MDM :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime)
{
    /*  Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process.
     * The current implementation computes the localDamageTensor and stores it in
     * corresponding status variable.
     * This service is declared at StructuralNonlocalMaterial level.
     */
    MDMStatus *status = ( MDMStatus * ) this->giveStatus(gp);
    FloatMatrix tempDamageTensor;

    this->computeLocalDamageTensor(tempDamageTensor, strainVector, gp, atTime);
    status->setLocalDamageTensorForAverage(tempDamageTensor);
}


double
MDM :: computeWeightFunction(const FloatArray &src, const FloatArray &coord)
{
    // Bell shaped function decaying with the distance.

    double dist = src.distance(coord);

    if ( ( dist >= 0. ) && ( dist <= this->R ) ) {
        double help = ( 1. - dist * dist / ( R * R ) );
        return help * help;
    }

    return 0.0;
}


int
MDM :: MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep)
{
    int result = 0;
    FloatArray intVal, strainIncr(3);
    IntArray toMap(1);
    MDMStatus *status = ( MDMStatus * ) this->giveStatus(gp);

    toMap.at(1) = ( int ) IST_MicroplaneDamageValues;

#ifndef MDM_MAPPING_DEBUG
    this->mapper.init(oldd, toMap, gp, tStep);
    result = mapper.mapVariable(intVal, gp, IST_MicroplaneDamageValues, tStep);
#else
    if ( mapperType == mdm_cpt ) {
        this->mapper2.init(oldd, toMap, gp, tStep);
        result = mapper2.mapVariable(intVal, gp, IST_MicroplaneDamageValues, tStep);
    } else if ( mapperType == mdm_sft ) {
        this->mapperSFT.init(oldd, toMap, gp, tStep);
        result = mapperSFT.mapVariable(intVal, gp, IST_MicroplaneDamageValues, tStep);
    } else if ( mapperType == mdm_lst ) {
        this->mapperLST.init(oldd, toMap, gp, tStep);
        result = mapperLST.mapVariable(intVal, gp, IST_MicroplaneDamageValues, tStep);
    } else {
        _error("MMI_map: unsupported Mapper id");
    }

#endif

    if ( formulation == COMPLIANCE_DAMAGE ) {
        for ( int i = 1; i <= intVal.giveSize(); i++ ) {
            if ( intVal.at(i) < 1.0 ) {
                intVal.at(i) = 1.0;
            }
        }
    } else {
        for ( int i = 1; i <= intVal.giveSize(); i++ ) {
            if ( intVal.at(i) < 0.0 ) {
                intVal.at(i) = 0.0;
            }

            if ( intVal.at(i) > 1.0 ) {
                intVal.at(i) = 1.0;
            }
        }
    }

    if ( result ) {
        status->setMicroplaneTempDamageValues(intVal);
    }

    // map stress, since it is necessary for keeping the
    // trace of stress (sv)

    toMap.resize(2);
    toMap.at(1) = ( int ) IST_StrainTensor;
    toMap.at(2) = ( int ) IST_StressTensor;
    this->mapper2.init(oldd, toMap, gp, tStep);

    result = mapper2.mapVariable(intVal, gp, IST_StressTensor, tStep);
    if ( result ) {
        status->letTempStressVectorBe(intVal);
    }

    result = mapper2.mapVariable(intVal, gp, IST_StrainTensor, tStep);
    if ( result ) {
        status->letTempStrainVectorBe(intVal);
    }

    status->updateYourself(tStep);

    return result;
}




int
MDM :: MMI_update(GaussPoint *gp,  TimeStep *tStep, FloatArray *estrain)
{
    int result = 1;
    FloatArray intVal, strain;
    MDMStatus *status = ( MDMStatus * ) this->giveStatus(gp);

    // now update all internal vars accordingly
    strain = status->giveStrainVector();
    this->giveRealStressVector(intVal, ReducedForm, gp, strain, tStep);
    this->updateYourself(gp, tStep);
    return result;
}

int
MDM :: MMI_finish(TimeStep *tStep)
{
#ifndef MDM_MAPPING_DEBUG
    this->mapper.finish(tStep);
#else
    this->mapperSFT.finish(tStep);
    this->mapperLST.finish(tStep);
#endif
    this->mapper2.finish(tStep);
    return 1;
}




Interface *
MDM :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return ( StructuralNonlocalMaterialExtensionInterface * ) this;
    } else if ( type == MaterialModelMapperInterfaceType ) {
        return ( MaterialModelMapperInterface * ) this;
    } else {
        return NULL;
    }
}


#ifdef __PARALLEL_MODE
int
MDM :: packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
{
    MDMStatus *status = ( MDMStatus * ) this->giveStatus(ip);

    this->buildNonlocalPointTable(ip);
    this->updateDomainBeforeNonlocAverage(stepN);

    return status->giveLocalDamageTensorForAveragePtr()->packToCommBuffer(buff);
}

int
MDM :: unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
{
    int result;
    MDMStatus *status = ( MDMStatus * ) this->giveStatus(ip);
    FloatMatrix _LocalDamageTensorForAverage;

    result = _LocalDamageTensorForAverage.unpackFromCommBuffer(buff);
    status->setLocalDamageTensorForAverage(_LocalDamageTensorForAverage);
    return result;
}

int
MDM :: estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip)
{
    //
    // Note: status localStrainVectorForAverage memeber must be properly sized!
    //
    if ( nonlocal ) {
        FloatMatrix _help(nsd, nsd);
        return _help.givePackSize(buff);
    } else {
        return 0;
    }
}

double
MDM :: predictRelativeComputationalCost(GaussPoint *gp)
{
    //
    // The values returned come from mesurement
    // do not change them unless you know what are you doing
    //
    double cost = 1.5;

    if ( nsd == 2 ) {
        cost = 1.5;
    } else if ( nsd == 3 )  {
        cost = 1.8;
    }

    if ( nonlocal ) {
        MDMStatus *status = ( MDMStatus * ) this->giveStatus(gp);
        int size = status->giveIntegrationDomainList()->size();
        // just a guess (size/10) found optimal
        // cost *= (1.0 + (size/10)*0.5);
        cost *= ( 1.0 + size / 15.0 );
    }

    return cost;
}

#endif



















MDMStatus :: MDMStatus(int n, int nsd, int nmplanes, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g), StructuralNonlocalMaterialStatusExtensionInterface(), Psi(nmplanes), PsiTemp(nmplanes), DamageTensor(nsd, nsd), DamageTensorTemp(nsd, nsd), tempDamageTensorEigenValues(nsd), damageTensorEigenValues(nsd), tempDamageTensorEigenVectors(nsd, nsd), damageTensorEigenVectors(nsd, nsd)
{
    int i;
    for ( i = 1; i <= nsd; i++ ) {
        damageTensorEigenValues.at(i) = tempDamageTensorEigenValues.at(i) = 1.0;
        tempDamageTensorEigenVectors.at(i, i) = damageTensorEigenVectors.at(i, i) = 1.0;
    }
}


MDMStatus :: ~MDMStatus() { }

contextIOResultType
MDMStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    //if (stream == NULL) StructuralMaterialStatus::_error ("saveContex : can't write into NULL stream");
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = Psi.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = DamageTensor.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = damageTensorEigenValues.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = damageTensorEigenVectors.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
MDMStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = Psi.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = DamageTensor.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = damageTensorEigenValues.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = damageTensorEigenVectors.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

void
MDMStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);

    Psi = PsiTemp;
    DamageTensor = DamageTensorTemp;
    damageTensorEigenValues = tempDamageTensorEigenValues;
    damageTensorEigenVectors = tempDamageTensorEigenVectors;
}

void
MDMStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    PsiTemp = Psi;
    DamageTensorTemp = DamageTensor;
    tempDamageTensorEigenValues = damageTensorEigenValues;
    tempDamageTensorEigenVectors = damageTensorEigenVectors;
}

void
MDMStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int i, j, n;

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");

    n = Psi.giveSize();
    fprintf(file, " compliance on microplanes: ");
    for ( i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", Psi.at(i) );
    }

    n = damageTensorEigenVectors.giveNumberOfRows();
    fprintf(file, ", complianceTensorEigenVectors ");
    for ( i = 1; i <= n; i++ ) {
        fprintf(file, "{");
        for ( j = 1; j <= n; j++ ) {
            fprintf( file, " % .4e", damageTensorEigenVectors.at(j, i) );
        }

        fprintf(file, "}");
    }

    fprintf(file, "}");

    fprintf(file, ", complianceTensorEigenValues ");
    for ( i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", damageTensorEigenValues.at(i) );
    }

    fprintf(file, "}\n");
}

Interface *
MDMStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return this;
    } else {
        return NULL;
    }
}
} // end namespace oofem
//