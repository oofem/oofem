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



#include "propagationlaw.h"
#include "enrichmentdomain.h"
#include "tipinfo.h"
#include "classfactory.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "spatiallocalizer.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "enrichmentitem.h"
#include "feinterpol.h"
#include "xfemmanager.h"

#include "XFEMDebugTools.h"

namespace oofem {
REGISTER_PropagationLaw(PLDoNothing)
REGISTER_PropagationLaw(PLCrackPrescribedDir)

PropagationLaw :: PropagationLaw() { }

PropagationLaw :: ~PropagationLaw() { }

void PLDoNothing :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);
}

IRResultType PLCrackPrescribedDir :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, mAngle, _IFT_PLCrackPrescribedDir_Dir);
    IR_GIVE_FIELD(ir, mIncrementLength, _IFT_PLCrackPrescribedDir_IncLength);

    //    printf("In PLCrackPrescribedDir :: initializeFrom: mAngle: %e mIncrementLength: %e\n", mAngle, mIncrementLength);

    return IRRT_OK;
}

void PLCrackPrescribedDir :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);

    input.setField(mAngle, _IFT_PLCrackPrescribedDir_Dir);
    input.setField(mIncrementLength, _IFT_PLCrackPrescribedDir_IncLength);
}

void PLCrackPrescribedDir :: propagateInterfaces(Domain &iDomain, EnrichmentDomain &ioEnrDom)
{
    // Fetch crack tip data
    std :: vector< TipInfo >tipInfo;
    ioEnrDom.giveTipInfos(tipInfo);

    int tipIndex = 1;
    double angleRad = mAngle * M_PI / 180.0;
    FloatArray dir = {cos(angleRad), sin(angleRad)};
    dir.normalize(); ///@todo Why? /Mikael


    std :: vector< TipPropagation >tipPropagations;
    TipPropagation tipProp;
    tipProp.mTipIndex = tipIndex;
    tipProp.mPropagationDir = dir;
    tipProp.mPropagationLength = mIncrementLength;
    tipPropagations.push_back(tipProp);

    ioEnrDom.propagateTips(tipPropagations);
}

} // end namespace oofem
