/*
 * plmaterialforce.C
 *
 *  Created on: Nov 14, 2014
 *      Author: svennine
 */

#include "plmaterialforce.h"

#include "xfem/propagationlaw.h"
#include "xfem/tipinfo.h"
#include "classfactory.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "spatiallocalizer.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "Materials/structuralms.h"
#include "xfem/enrichmentitem.h"
#include "feinterpol.h"
#include "xfem/xfemmanager.h"
#include "xfem/matforceevaluator.h"

#include "engngm.h"

namespace oofem {
REGISTER_PropagationLaw(PLMaterialForce)

PLMaterialForce :: PLMaterialForce():
    mRadius(0.0),
    mIncrementLength(0.0),
    mCrackPropThreshold(0.0),
    mpMaterialForceEvaluator( new MaterialForceEvaluator() )
{}

PLMaterialForce :: ~PLMaterialForce()
{}

IRResultType PLMaterialForce :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, mRadius,                          _IFT_PLMaterialForce_Radius);
    IR_GIVE_FIELD(ir, mIncrementLength,                 _IFT_PLMaterialForce_IncLength);
    IR_GIVE_OPTIONAL_FIELD(ir, mCrackPropThreshold,     _IFT_PLMaterialForce_CrackPropThreshold);

    return IRRT_OK;
}

void PLMaterialForce :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);

    input.setField(mRadius,                     _IFT_PLMaterialForce_Radius);
    input.setField(mIncrementLength,            _IFT_PLMaterialForce_IncLength);
    input.setField(mCrackPropThreshold,         _IFT_PLMaterialForce_CrackPropThreshold);
}

bool PLMaterialForce :: propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp)
{
//    printf("Entering PLMaterialForce :: propagateInterface().\n");

    if ( !iEnrFront.propagationIsAllowed() ) {
        return false;
    }

    // Fetch crack tip data
    const TipInfo &tipInfo = iEnrFront.giveTipInfo();

    // Check if the tip is located in the domain
    SpatialLocalizer *localizer = iDomain.giveSpatialLocalizer();
    FloatArray lCoords, closest;
    localizer->giveElementClosestToPoint(lCoords, closest, tipInfo.mGlobalCoord);

    if(closest.distance(tipInfo.mGlobalCoord) > 1.0e-9) {
//        printf("Tip is outside all elements.\n");
        return false;
    }


    FloatArray matForce;
    TimeStep *tStep = iDomain.giveEngngModel()->giveCurrentStep();
    mpMaterialForceEvaluator->computeMaterialForce(matForce, iDomain, tipInfo, tStep, mRadius);

//    printf("matForce: "); matForce.printYourself();

    if(matForce.giveSize() == 0) {
        return false;
    }

    double forceNorm = matForce.computeNorm();
//    printf("forceNorm: %e\n", forceNorm);

    if(forceNorm < mCrackPropThreshold || forceNorm < 1.0e-20) {
        return false;
    }

//    printf("Propagating crack.\n");
//    printf("Tip coord: "); tipInfo.mGlobalCoord.printYourself();

    FloatArray dir(matForce);
    dir.times(1.0/forceNorm);
//    printf("dir: "); dir.printYourself();

    const double cosAngTol = 1.0/sqrt(2.0);
    if(tipInfo.mTangDir.dotProduct(dir) < cosAngTol) {
        // Do not allow sharper turns than 45 degrees

        if( tipInfo.mNormalDir.dotProduct(dir) > 0.0 ) {
            dir = tipInfo.mTangDir;
            dir.add(tipInfo.mNormalDir);
            dir.normalize();
        }
        else {
//            dir = tipInfo.mNormalDir;
//            dir.times(-1.0);
            dir = tipInfo.mTangDir;
            dir.add(-1.0,tipInfo.mNormalDir);
            dir.normalize();
        }

        printf("//////////////////////////////////////////// Resticting crack propagation direction.\n");
//        printf("tipInfo.mTangDir: "); tipInfo.mTangDir.printYourself();
//        printf("dir: "); dir.printYourself();
    }

    // Fill up struct
    oTipProp.mTipIndex = tipInfo.mTipIndex;
    oTipProp.mPropagationDir = dir;
    oTipProp.mPropagationLength = mIncrementLength;


    return true;
}

} /* namespace oofem */
