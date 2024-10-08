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

#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "intmatbilinearcz.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(IntMatBilinearCZ);


///@todo - need to rearrange traction and stiffness matrix so the first component is normal
IntMatBilinearCZStatus :: IntMatBilinearCZStatus(GaussPoint *g) : StructuralInterfaceMaterialStatus(g)
{}


void IntMatBilinearCZStatus :: initTempStatus()
{ }

void IntMatBilinearCZStatus :: updateYourself(TimeStep *tStep)
{
    mTractionOld = mTractionNew;
    mJumpOld     = mJumpNew;
    mDamageOld   = mDamageNew;

    jump = mJumpNew;

    mPlastMultIncOld = mPlastMultIncNew;

    StructuralInterfaceMaterialStatus ::updateYourself(tStep);
}


void IntMatBilinearCZStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    StructuralInterfaceMaterialStatus :: copyStateVariables(iStatus);

    const IntMatBilinearCZStatus &structStatus = static_cast< const IntMatBilinearCZStatus & >(iStatus);

    mDamageNew   = structStatus.mDamageNew;
    mDamageOld   = structStatus.mDamageOld;
    mTractionOld = structStatus.mTractionOld;
    mTractionNew = structStatus.mTractionNew;
    mJumpOld     = structStatus.mJumpOld;
    mJumpNew     = structStatus.mJumpNew;

    mPlastMultIncNew = structStatus.mPlastMultIncNew;
    mPlastMultIncOld = structStatus.mPlastMultIncOld;
}

void IntMatBilinearCZStatus :: addStateVariables(const MaterialStatus &iStatus)
{
    OOFEM_ERROR("not implemented.");
}


IntMatBilinearCZ :: IntMatBilinearCZ(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{ }


int IntMatBilinearCZ :: checkConsistency()
{
    return 1;
}

FloatArrayF<3> IntMatBilinearCZ :: giveFirstPKTraction_3d(const FloatArrayF<3> &jump, const FloatMatrixF<3,3> &F, GaussPoint *gp, TimeStep *tStep) const
{
    double maxDamage = 0.99999999;

    IntMatBilinearCZStatus *status = static_cast< IntMatBilinearCZStatus * >( this->giveStatus(gp) );

    status->mJumpNew = jump;

    auto jumpInc = status->mJumpNew - status->mJumpOld;

    auto tractionTrial = status->mTractionOld + mPenaltyStiffness * jumpInc;

    double TTrNormal = tractionTrial.at(3);
    double TTrTang = sqrt( pow(tractionTrial.at(1), 2.0) + pow(tractionTrial.at(2), 2.0) );
    double phiTr = computeYieldFunction(TTrNormal, TTrTang);

    if ( status->mDamageOld > maxDamage ) {
        status->mDamageNew = maxDamage;
        status->mDamageOld = maxDamage;
        status->mPlastMultIncNew = 0.0;
        FloatArrayF<3> answer;
        status->mTractionNew = answer;

        status->letTempJumpBe(jump);
        status->letTempFirstPKTractionBe(answer);
        status->letTempTractionBe(answer);

        return answer;
    }


    if ( phiTr < 0.0 ) {
        status->mDamageNew = status->mDamageOld;
        status->mPlastMultIncNew = 0.0;
        status->mTractionNew = tractionTrial;


        status->letTempJumpBe(jump);
        status->letTempFirstPKTractionBe(tractionTrial);
        status->letTempTractionBe(tractionTrial);

        auto answer = tractionTrial;
//        answer.beScaled( ( 1.0 - status->mDamageNew ), answer );
        answer.at(1) *= ( 1.0 - status->mDamageNew );
        answer.at(2) *= ( 1.0 - status->mDamageNew );

        if ( answer.at(3) > 0.0 ) {
            answer.at(3) *= ( 1.0 - status->mDamageNew );
        }

        return answer;
    } else {
        // Iterate to find plastic strain increment.
        int maxIter = 50;
        int minIter = 1;
        double absTol = 1.0e-9; // Absolute error tolerance
        double relTol = 1.0e-9; // Relative error tolerance
        double eps = 1.0e-12; // Small value for perturbation when computing numerical Jacobian
        double plastMultInc = 0.0;
        double initialRes = 0.0;

        for ( int iter = 0; iter < maxIter; iter++ ) {
            // Evaluate residual (i.e. yield function)
            auto answer = computeTraction(tractionTrial, plastMultInc);

            double TNormal = answer.at(3);
            double TTang = sqrt( pow(answer.at(1), 2.0) + pow(answer.at(2), 2.0) );
            double phi = computeYieldFunction(TNormal, TTang);

            //          if(iter > 20) {
            //              printf("iter: %d res: %e\n", iter, fabs(phi) );
            //          }

            if ( iter == 0 ) {
                initialRes = fabs(phi);
                initialRes = max(initialRes, 1.0e-12);
            }

            if ( (iter >= minIter && fabs(phi) < absTol) || ( iter >= minIter && ( fabs(phi) / initialRes ) < relTol ) ) {
                // Add damage evolution
                double S = mGIc / mSigmaF;
                status->mPlastMultIncNew = plastMultInc;

                double damageInc = status->mPlastMultIncNew / S;
                status->mDamageNew = status->mDamageOld + damageInc;

                if ( status->mDamageNew > maxDamage ) {
                    status->mDamageNew = maxDamage;
                }

                status->mTractionNew = answer;

                // Jim
                status->letTempJumpBe(jump);
                status->letTempFirstPKTractionBe(answer);
                status->letTempTractionBe(answer);

                if ( mSemiExplicit ) {
                    answer = computeTraction(tractionTrial, status->mPlastMultIncOld);

//                    answer.beScaled( ( 1.0 - status->mDamageOld ), answer );

                    answer.at(1) *= ( 1.0 - status->mDamageOld );
                    answer.at(2) *= ( 1.0 - status->mDamageOld );

                    if ( answer.at(3) > 0.0 ) {
                        answer.at(3) *= ( 1.0 - status->mDamageOld );
                    }
               } else {
//                    answer.beScaled( ( 1.0 - status->mDamageNew ), answer );

                    answer.at(1) *= ( 1.0 - status->mDamageNew );
                    answer.at(2) *= ( 1.0 - status->mDamageNew );

                    if ( answer.at(3) > 0.0 ) {
                        answer.at(3) *= ( 1.0 - status->mDamageNew );
                    }
                }

                return answer;
            }

            // Numerical Jacobian
            auto tractionPert = computeTraction(tractionTrial, plastMultInc + eps);
            double TNormalPert = tractionPert.at(3);
            double TTangPert = sqrt( pow(tractionPert.at(1), 2.0) + pow(tractionPert.at(2), 2.0) );
            double phiPert = computeYieldFunction(TNormalPert, TTangPert);

            double Jac = ( phiPert - phi ) / eps;
            plastMultInc -= ( 1.0 / Jac ) * phi;
        }
    }

    OOFEM_ERROR("No convergence in.");
}

FloatMatrixF<3,3> IntMatBilinearCZ :: give3dStiffnessMatrix_dTdj(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_WARNING("not implemented. Use numerical Jacobian instead.");
    return this->give3dStiffnessMatrix_dTdj_Num(gp, tStep);
}

double IntMatBilinearCZ :: computeYieldFunction(double iTractionNormal, double iTractionTang) const
{
    return mSigmaF * pow(fabs(iTractionTang) / ( mGamma * mSigmaF ), 2.0)
            + ( mSigmaF / mGamma ) * ( mGamma - 2.0 * mMu ) * pow( ( max(iTractionNormal, 0.0) ) / mSigmaF, 2.0 )
            - ( 1.0 / mGamma ) * ( mGamma * mSigmaF - 2.0 * mMu * iTractionNormal );
}

FloatArrayF<3> IntMatBilinearCZ :: computeTraction(const FloatArrayF<3> &iTTrial, double iPlastMultInc) const
{
    // Note that the traction vector is assumed to be decomposed into its tangential and normal part,
    // with tangential directions (1,2) and normal direction (3).

    double gammaF = mGIIc / mGIc;

    FloatArrayF<3> oT;
    // Tangential part
    oT.at(1) = iTTrial.at(1) / ( 1.0 + ( 2.0 * mPenaltyStiffness * gammaF ) * iPlastMultInc / ( mGamma * mGamma * mSigmaF ) );
    oT.at(2) = iTTrial.at(2) / ( 1.0 + ( 2.0 * mPenaltyStiffness * gammaF ) * iPlastMultInc / ( mGamma * mGamma * mSigmaF ) );

    // Normal part
    if ( iTTrial.at(3) <= 0.0 ) {    // TODO: Check if-statement
        oT.at(3) = iTTrial.at(3);
    } else {
        oT.at(3) = iTTrial.at(3) / ( 1.0 + 2.0 * mPenaltyStiffness * iPlastMultInc / mSigmaF );
    }
    return oT;
}

void IntMatBilinearCZ :: initializeFrom(InputRecord &ir)
{
    StructuralInterfaceMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, mPenaltyStiffness, _IFT_IntMatBilinearCZ_PenaltyStiffness);

    IR_GIVE_FIELD(ir, mGIc, _IFT_IntMatBilinearCZ_g1c);

    mGIIc = mGIc;                                               //Defaults to GIc
    IR_GIVE_OPTIONAL_FIELD(ir, mGIIc, _IFT_IntMatBilinearCZ_g2c);

    IR_GIVE_FIELD(ir, mSigmaF, _IFT_IntMatBilinearCZ_sigf);

    IR_GIVE_FIELD(ir, mMu, _IFT_IntMatBilinearCZ_mu);

    IR_GIVE_FIELD(ir, mGamma, _IFT_IntMatBilinearCZ_gamma);

    if ( ir.hasField(_IFT_IntMatBilinearCZ_semiexplicit) ) {
        mSemiExplicit = true;
        printf("In IntMatBilinearCZ::initializeFrom: Semi-explicit time integration activated.\n");
    }
}

void IntMatBilinearCZ :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(mPenaltyStiffness, _IFT_IntMatBilinearCZ_PenaltyStiffness);

    input.setField(mGIc, _IFT_IntMatBilinearCZ_g1c);
    input.setField(mGIIc, _IFT_IntMatBilinearCZ_g2c);

    input.setField(mSigmaF, _IFT_IntMatBilinearCZ_sigf);
    input.setField(mMu, _IFT_IntMatBilinearCZ_mu);
    input.setField(mGamma, _IFT_IntMatBilinearCZ_gamma);
}

void IntMatBilinearCZ :: printYourself()
{
    printf("\nInitializing IntMatBilinearCZ:\n");
    printf("mPenaltyStiffness: %e\n", mPenaltyStiffness);
    printf("mGIc: %e\n", mGIc);
    printf("mGIIc: %e\n", mGIIc);
    printf("mSigmaF: %e\n", mSigmaF);
    printf("mMu: %e\n", mMu);
    printf("mGamma: %e\n\n", mGamma);
}
} /* namespace oofem */
