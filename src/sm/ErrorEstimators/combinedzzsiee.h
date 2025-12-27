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

#ifndef combinedzzsiee_h
#define combinedzzsiee_h

#include "sm/ErrorEstimators/zzerrorestimator.h"
#include "sm/ErrorEstimators/scalarerrorindicator.h"
#include "sm/ErrorEstimators/directerrorindicatorrc.h"

#define _IFT_CombinedZZSIErrorEstimator_Name "zzscalar"

namespace oofem {
/**
 * The implementation of combined criteria: Zienkiewicz Zhu Error Estimator for elastic regime and
 * scalar error indicator in non-linear regime.
 * The basic task is to evaluate the stress error on associated domain.
 * The algorithm is written in general way, so it is possible to to evaluate
 * different errors (for example temperature error). Then corresponding
 * member attribute identifying the type of quantity used should be declared and initialized
 * (for example in instanciateYourself() service). Then the modification is required
 * only when requesting element contributions.
 *
 */
class CombinedZZSIErrorEstimator : public ErrorEstimator
{
protected:
    ZZErrorEstimator zzee;
    ScalarErrorIndicator siee;

public:
    /// Constructor
    CombinedZZSIErrorEstimator(int n, Domain * d) : ErrorEstimator(n, d), zzee(n, d), siee(n, d) {
        eeType = EET_CZZSI;
    }
    /// Destructor
    virtual ~CombinedZZSIErrorEstimator() { }

    double giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep) override;
    double giveValue(EE_ValueType type, TimeStep *tStep) override;

    int estimateError(EE_ErrorMode mode, TimeStep *tStep) override;
    RemeshingCriteria *giveRemeshingCrit() override;
    void initializeFrom(InputRecord &ir) override;
    const char *giveInputRecordName() const override { return _IFT_CombinedZZSIErrorEstimator_Name; }
    const char *giveClassName() const override { return "CombinedZZSIErrorEstimator"; }
    void setDomain(Domain *d) override;
};

/**
 * The class represent the corresponding remeshing criteria to CombinedZZSIErrorEstimator.
 * In regions, where the error indicator is larger than given threshold, the
 * mesh is refined according this indicator (currently linear interpolation).
 * Otherwise, the mesh size is determined from Zienkiewicz-Zhu remeshing criteria.
 * (Assumes that error is equally distributed between elements, then the requirement for max. permissible error
 * can be translated into placing a limit on the error on each element.)
 * The basic task is to evaluate the required mesh density (at nodes) on given domain,
 * based on informations provided by the compatible error estimator.
 *
 * The remeshing criteria is maintained by the corresponding error estimator. This is mainly due to fact, that is
 * necessary for given EE to create compatible RC. In our concept, the EE is responsible.
 */
class CombinedZZSIRemeshingCriteria : public RemeshingCriteria
{
protected:
    ZZRemeshingCriteria zzrc;
    DirectErrorIndicatorRC dirc;

public:
    /// Constructor
    CombinedZZSIRemeshingCriteria(int n, ErrorEstimator * e);
    /// Destructor
    virtual ~CombinedZZSIRemeshingCriteria() { }

    double giveRequiredDofManDensity(int num, TimeStep *tStep, int relative = 0) override;
    double giveDofManDensity(int num) override;
    RemeshingStrategy giveRemeshingStrategy(TimeStep *tStep) override;
    int estimateMeshDensities(TimeStep *tStep) override;
    void initializeFrom(InputRecord &ir) override;
    const char *giveInputRecordName() const override { return nullptr; }
    const char *giveClassName() const override { return "CombinedZZSIRemeshingCriteria"; }
    void setDomain(Domain *d) override;
};
} // end namespace oofem
#endif // combinedzzsiee_h
