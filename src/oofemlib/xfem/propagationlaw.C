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



#include "propagationlaw.h"
#include "enrichmentdomain.h"
#include "tipinfo.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {

REGISTER_PropagationLaw(PLDoNothing)
REGISTER_PropagationLaw(PLCrackPrescribedDir)

PropagationLaw::PropagationLaw() {

}

PropagationLaw::~PropagationLaw() {

}

IRResultType PLCrackPrescribedDir :: initializeFrom(InputRecord *ir) {

    const char *__proc = "initializeFrom";
    IRResultType result;

    IR_GIVE_FIELD(ir, mAngle, _IFT_PLCrackPrescribedDir_Dir);
    IR_GIVE_FIELD(ir, mIncrementLength, _IFT_PLCrackPrescribedDir_IncLength);

    printf("In PLCrackPrescribedDir :: initializeFrom: mAngle: %e mIncrementLength: %e\n", mAngle, mIncrementLength );

	return IRRT_OK;
}

void PLCrackPrescribedDir::propagateInterfaces(EnrichmentDomain &ioEnrDom) {
	printf("Entering PLCrackPrescribedDir::propagateInterfaces().\n");

	// Fetch crack tip data
	std::vector<TipInfo> tipInfo;
	ioEnrDom.giveTipInfos(tipInfo);
//	printf("tipInfo.size(): %lu\n", tipInfo.size());


	int tipIndex = 1;
	FloatArray dir;
	double angleRad = mAngle*M_PI/180.0;
	dir.setValues(2, cos(angleRad), sin(angleRad));
	dir.normalize();


	std::vector<TipPropagation> tipPropagations;
	TipPropagation tipProp;
	tipProp.mTipIndex = tipIndex;
	tipProp.mPropagationDir = dir;
	tipProp.mPropagationLength = mIncrementLength;
	tipPropagations.push_back(tipProp);

	ioEnrDom.propagateTips(tipPropagations);
}

} // end namespace oofem

