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

#ifndef XFEMSTRUCTUREMANAGER_H_
#define XFEMSTRUCTUREMANAGER_H_

#include "xfem/xfemmanager.h"

#include <memory>

///@name Input fields for XfemManager
//@{
#define _IFT_XfemStructureManager_Name "xfemstructuremanager"
#define _IFT_XfemStructureManager_splitCracks "splitcracks"
#define _IFT_XfemStructureManager_nonstandardCZ "nonstandardcz"
#define _IFT_XfemStructureManager_minCrackLength "mincracklength"
#define _IFT_XfemStructureManager_crackMergeTol "crackmergetol"
//@}

namespace oofem {

class MaterialForceEvaluator;

/**
 * XfemStructureManager: XFEM manager with extra functionality
 * specific for the sm module.
 *
 * @author Erik Svenning
 * @date Apr 28, 2014
 */
class OOFEM_EXPORT XfemStructureManager : public XfemManager
{
public:
    XfemStructureManager(Domain *domain);
    virtual ~XfemStructureManager();

    /// Initializes receiver according to object description stored in input record.
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    int instanciateYourself(DataReader &dr) override;
    const char *giveClassName() const override { return "XfemStructureManager"; }
    const char *giveInputRecordName() const override { return _IFT_XfemStructureManager_Name; }
    virtual void postInitialize() override;

    void propagateFronts(bool &oAnyFronHasPropagated) override;

    /**
     * Update enrichment items (level sets).
     */
    void updateYourself(TimeStep *tStep) override;

    void splitCracks();

    void removeShortCracks();

    bool tipsHaveOppositeDirection(EnrichmentFront *iEf1, EnrichmentFront *iEf2);
    void mergeCloseCracks();

    bool giveUseNonStdCz() const { return mNonstandardCz; }

    /// Compute the total length of all cracks in the domain.
    double computeTotalCrackLength();

protected:

    /**
     * If cracks should be splitted at intersections as a pre-processing step.
     */
    bool mSplitCracks;

    /**
     * If a non-standard cohesive zone formulation should be used.
     */
    bool mNonstandardCz;

    /**
     * Cracks shorter than this length are automatically removed.
     */
    double mMinCrackLength;

    /**
     * Cracks with tips closer than this distance are automatically merged.
     */
    double mCrackMergeTol;

    /**
     * Evaluator for material forces.
     */
    std :: unique_ptr< MaterialForceEvaluator > mpMatForceEvaluator;
};
} /* namespace oofem */

#endif /* XFEMSTRUCTUREMANAGER_H_ */
