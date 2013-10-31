/*
 * intmatbilinearcz.C
 *
 *  Created on: Oct 20, 2013
 *      Author: svennine
 */

#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "intmatbilinearcz.h"

namespace oofem {

REGISTER_Material( IntMatBilinearCZ );

//IntMatBilinearCZFagerstromStatus :: IntMatBilinearCZFagerstromStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
IntMatBilinearCZStatus :: IntMatBilinearCZStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g),
mDamageNew(0.0),
mDamageOld(0.0)
{
	mTractionOld.resize(3);
	mTractionOld.zero();

	mTractionNew.resize(3);
	mTractionNew.zero();

	mJumpOld.resize(3);
	mJumpOld.zero();

	mJumpNew.resize(3);
	mJumpNew.zero();

	mDamageOld = 0.0;
	mDamageNew = 0.0;
}


IntMatBilinearCZStatus :: ~IntMatBilinearCZStatus()
{

}

void IntMatBilinearCZStatus :: initTempStatus()
{

}

void IntMatBilinearCZStatus :: updateYourself(TimeStep *atTime)
{
	mTractionOld 	= mTractionNew;
	mJumpOld 		= mJumpNew;
	mDamageOld		= mDamageNew;
}


IntMatBilinearCZ::IntMatBilinearCZ(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{

}

IntMatBilinearCZ::~IntMatBilinearCZ()
{

}

int IntMatBilinearCZ::checkConsistency()
{
	return 1;
}

/*
void IntMatBilinearCZ::give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                                 GaussPoint *gp, TimeStep *atTime)
{

}
*/

void IntMatBilinearCZ::giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                     const FloatMatrix &F, TimeStep *tStep)
{
    IntMatBilinearCZStatus *status = static_cast< IntMatBilinearCZStatus * >( this->giveStatus(gp) );

    status->mJumpNew = jump;

    FloatArray jumpInc;
    jumpInc.beDifferenceOf( status->mJumpNew, status->mJumpOld );

    FloatArray tractionTrial = status->mTractionOld;
    tractionTrial.add(mPenaltyStiffness, jumpInc);

    double TTrNormal 	= tractionTrial.at(3);
    double TTrTang 		= sqrt( pow(tractionTrial.at(1),2.0) + pow(tractionTrial.at(2),2.0) );
    double phiTr = computeYieldFunction(TTrNormal, TTrTang);

    double damageTol = 1.0e-6;
    if( status->mDamageOld > (1.0-damageTol) ) {
    	status->mDamageNew = 1.0;
    	answer.zero();
    	status->mTractionNew = answer;
    	return;
    }


	answer = tractionTrial;

    if( phiTr < 0.0 ) {
    	status->mDamageNew = status->mDamageOld;
    	answer.beScaled((1.0-status->mDamageNew), answer);

    	status->mTractionNew = answer;
    	return;
    }
    else {

    	// Iterate to find plastic strain increment.
    	int maxIter = 50;
    	double absTol = 1.0e-9;
    	double eps = 1.0e-9; // Small value for perturbation when computing numerical Jacobian
		double plastMultInc = 0.0;

		for(int iter = 0; iter < maxIter; iter++) {

			// Evaluate residual (i.e. yield function)
			computeTraction( answer, tractionTrial, plastMultInc );

    	    double TNormal 	= answer.at(3);
    	    double TTang 		= sqrt( pow(answer.at(1),2.0) + pow(answer.at(2),2.0) );
    	    double phi = computeYieldFunction(TNormal, TTang);

//    	    printf("iter: %d res: %e\n", iter, fabs(phi) );

    	    if( fabs(phi) < absTol ) {

    	    	// Add damage evolution
    	    	double S = mGIc/mSigmaF;
    	    	double damageInc = plastMultInc/S;
    	    	status->mDamageNew = status->mDamageOld + damageInc;

    	    	if(status->mDamageNew > 1.0) {
    	    		status->mDamageNew = 1.0;
    	    	}

    	    	answer.beScaled((1.0-status->mDamageNew), answer);
    	    	status->mTractionNew = answer;

    	    	return;
    	    }

    	    // Numerical Jacobian
    	    FloatArray tractionPert(3);
			computeTraction( tractionPert, tractionTrial, plastMultInc+eps );
    	    double TNormalPert 		= tractionPert.at(3);
    	    double TTangPert 		= sqrt( pow(tractionPert.at(1),2.0) + pow(tractionPert.at(2),2.0) );
    	    double phiPert = computeYieldFunction(TNormalPert, TTangPert);

    	    double Jac = (phiPert - phi)/eps;
    	    plastMultInc -= (1.0/Jac)*phi;
    	}

    }

    OOFEM_ERROR("Warning: No convergence in IntMatBilinearCZ::giveFirstPKTraction_3d().\n");
}

void IntMatBilinearCZ::give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
	OOFEM_ERROR("IntMatBilinearCZ::give3dStiffnessMatrix_dTdj is not implemented. Use numerical Jacobian instead.");
}

double IntMatBilinearCZ::computeYieldFunction(const double &iTractionNormal, const double &iTractionTang) const
{
	return (
			mSigmaF*pow(fabs(iTractionTang)/(mGamma*mSigmaF) , 2.0)
		 + (mSigmaF/mGamma)*(mGamma-2.0*mMu)*pow( (max(iTractionNormal,0.0))/mSigmaF, 2.0)
		 - (1.0/mGamma)*( mGamma*mSigmaF - 2.0*mMu*iTractionNormal )  );
}

void IntMatBilinearCZ::computeTraction( FloatArray &oT, const FloatArray &iTTrial, const double &iPlastMultInc ) const
{
	// Note that the traction vector is assumed to be decomposed into its tangential and normal part,
	// with tangential directions (1,2) and normal direction (3).

//	double gammaF = mGamma*mGIIc/mGIc;
	double gammaF =  mGIIc/mGIc;

	// Tangential part
	oT.at(1) = iTTrial.at(1)/( 1.0 + (2.0*mPenaltyStiffness*gammaF)*iPlastMultInc/(mGamma*mGamma*mSigmaF) );
	oT.at(2) = iTTrial.at(2)/( 1.0 + (2.0*mPenaltyStiffness*gammaF)*iPlastMultInc/(mGamma*mGamma*mSigmaF) );

	// Normal part
	if( iTTrial.at(3) <= 0.0 ) { // TODO: Check if-statement
		oT.at(3) = iTTrial.at(3);
	}
	else {
		oT.at(3) = iTTrial.at(3)/(1.0 + 2.0*mPenaltyStiffness*iPlastMultInc/mSigmaF);
	}

}

int IntMatBilinearCZ::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    IntMatBilinearCZStatus *status = static_cast< IntMatBilinearCZStatus * >( this->giveStatus(gp) );
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = status->mDamageNew;
        return 1;
    } else {
        return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
    }

}

InternalStateValueType IntMatBilinearCZ::giveIPValueType(InternalStateType type)
{
    return StructuralInterfaceMaterial :: giveIPValueType(type);
}

IRResultType IntMatBilinearCZ::initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";  // Required by IR_GIVE_FIELD macro
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, mPenaltyStiffness, _IFT_IntMatBilinearCZ_PenaltyStiffness);

    IR_GIVE_FIELD(ir, mGIc, _IFT_IntMatBilinearCZ_g1c);

    mGIIc = mGIc;						//Defaults to GIc
    IR_GIVE_OPTIONAL_FIELD(ir, mGIIc, _IFT_IntMatBilinearCZ_g2c);

    IR_GIVE_FIELD(ir, mSigmaF, _IFT_IntMatBilinearCZ_sigf);

    IR_GIVE_FIELD(ir, mMu, _IFT_IntMatBilinearCZ_mu);

    IR_GIVE_FIELD(ir, mGamma, _IFT_IntMatBilinearCZ_gamma);

    this->checkConsistency();
    this->printYourself();
    return IRRT_OK;
}

void IntMatBilinearCZ::giveInputRecord(DynamicInputRecord &input)
{
	StructuralInterfaceMaterial::giveInputRecord(input);

	input.setField(mPenaltyStiffness, _IFT_IntMatBilinearCZ_PenaltyStiffness);

	input.setField(mGIc, _IFT_IntMatBilinearCZ_g1c);
	input.setField(mGIIc, _IFT_IntMatBilinearCZ_g2c);

	input.setField(mSigmaF, _IFT_IntMatBilinearCZ_sigf);
	input.setField(mMu, _IFT_IntMatBilinearCZ_mu);
	input.setField(mGamma, _IFT_IntMatBilinearCZ_gamma);

}

void IntMatBilinearCZ::printYourself()
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
