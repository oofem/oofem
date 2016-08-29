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

#ifndef SRC_SM_XFEM_NUCLEATIONCRITERIA_NCPRINCIPALSTRESS_H_
#define SRC_SM_XFEM_NUCLEATIONCRITERIA_NCPRINCIPALSTRESS_H_

#define _IFT_NCPrincipalStress_Name "ncprincipalstress"
#define _IFT_NCPrincipalStress_StressThreshold "stressthreshold"
#define _IFT_NCPrincipalStress_InitialCrackLength "initialcracklength"

#include "xfem/nucleationcriterion.h"

#include <memory>

namespace oofem {

class NCPrincipalStress : public NucleationCriterion {
public:
	NCPrincipalStress(Domain *ipDomain);
	virtual ~NCPrincipalStress();

	virtual std::vector<std::unique_ptr<EnrichmentItem>> nucleateEnrichmentItems();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void appendInputRecords(DynamicDataReader &oDR);

    /// @return Class name of the receiver.
    virtual const char *giveClassName() const {return "NCPrincipalStress";}
    /// @return Input record name of the receiver.
    virtual const char *giveInputRecordName() const {return _IFT_NCPrincipalStress_Name;};

protected:
    double mStressThreshold;
    double mInitialCrackLength;

    /// If the initiated crack should cut exactly one element.
    bool mCutOneEl;

    /// Index of the cross section that the nucleation criterion applies to.
    int mCrossSectionInd;
};

} /* namespace oofem */

#endif /* SRC_SM_XFEM_NUCLEATIONCRITERIA_NCPRINCIPALSTRESS_H_ */
