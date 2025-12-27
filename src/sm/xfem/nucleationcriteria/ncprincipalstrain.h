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

#ifndef SRC_SM_XFEM_NUCLEATIONCRITERIA_NCPRINCIPALSTRAIN_H_
#define SRC_SM_XFEM_NUCLEATIONCRITERIA_NCPRINCIPALSTRAIN_H_

#define _IFT_NCPrincipalStrain_Name "ncprincipalstrain"
#define _IFT_NCPrincipalStrain_StrainThreshold "strainthreshold"
#define _IFT_NCPrincipalStrain_IncrementLength "incrementlength"
#define _IFT_NCPrincipalStrain_PropStrainThreshold "propagationstrainthreshold"
#define _IFT_NCPrincipalStrain_InitialCrackLength "initialcracklength"
#define _IFT_NCPrincipalStrain_CrossSectionIndex "csindex"

#include "xfem/nucleationcriterion.h"
#include <memory>

namespace oofem {

class NCPrincipalStrain : public NucleationCriterion
{
public:
    NCPrincipalStrain(Domain *ipDomain);
    virtual ~NCPrincipalStrain();

    std::vector<std::unique_ptr<EnrichmentItem>> nucleateEnrichmentItems() override;

    void initializeFrom(InputRecord &ir) override;

    void appendInputRecords(DynamicDataReader &oDR) override;

    const char *giveClassName() const override { return "NCPrincipalStrain"; }
    const char *giveInputRecordName() const override { return _IFT_NCPrincipalStrain_Name; }

protected:
    double mStrainThreshold;
    double mInitialCrackLength;
    double mIncrementLength;
    double mPropStrainThreshold;

    /// If the initiated crack should cut exactly one element.
    bool mCutOneEl;

    /// Index of the cross section that the nucleation criterion applies to.
    int mCrossSectionInd;
};

} /* namespace oofem */

#endif /* SRC_SM_XFEM_NUCLEATIONCRITERIA_NCPRINCIPALSTRAIN_H_ */
