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

#ifndef adaptlinearstatic_h
#define adaptlinearstatic_h

#include "sm/EngineeringModels/linearstatic.h"
#include "meshpackagetype.h"
#include "remeshingcrit.h"
#include "mesherinterface.h"

#include <memory>

///@name Input fields for AdaptiveLinearStatic
//@{
#define _IFT_AdaptiveLinearStatic_Name "adaptlinearstatic"
#define _IFT_AdaptiveLinearStatic_meshpackage "meshpackage"
#define _IFT_AdaptiveLinearStatic_adapt "adapt"
#define _IFT_AdaptiveLinearStatic_maxadaptsteps "maxadaptsteps"
//@}

namespace oofem {
/**
 * This class implements an adaptive linear static engineering problem.
 * Only single loading case is supported.
 *
 * The problem is analyzed on the initial mesh, the error is estimated and the required
 * mesh density is evaluated by the remeshing criteria and passed to the mesher
 * (typically written into mesher specific density file). This is the default behavior.
 *
 * When adaptivity is enabled (adapt 1) and the error exceeds the required threshold, the
 * mesher is asked to create the new mesh. If the mesher can provide the new mesh
 * (directly, as Subdivision, or by running external mesher and converter, see MesherInterface),
 * the new solution step is created on the new mesh, and the analysis followed by error
 * estimation is repeated, until the error is acceptable or the maximum number of adaptive
 * steps is reached. Due to linearity of a problem, the complete reanalysis is done on
 * each new mesh, no mapping of the solution is needed.
 * Solution steps represent a series of adaptive analyses (all steps have the same target time).
 */
class AdaptiveLinearStatic : public LinearStatic
{
protected:
    /// Meshing package used for refinements.
    MeshPackageType meshPackage;
    /// Mesher interface.
    std :: unique_ptr< MesherInterface >mesher;
    /// Flag indicating whether the adaptive remeshing loop is enabled.
    bool adaptFlag;
    /// Maximum number of adaptive remeshings.
    int maxAdaptSteps;
    /// Remeshing strategy determined by the remeshing criteria in the current step.
    RemeshingStrategy remeshingStrategy;

public:
    AdaptiveLinearStatic(int i, EngngModel *master = nullptr);
    virtual ~AdaptiveLinearStatic();

    void solveYourself() override;
    void updateYourself(TimeStep *tStep) override;
    TimeStep *giveNextStep() override;

    /**
     * Initializes the newly generated discretization state according to previous solution.
     * This process should typically include restoring old solution, instanciating newly
     * generated domain(s) and by mapping procedure.
     */
    int initializeAdaptive(int tStepNumber) override;
    void printOutputAt(FILE *file, TimeStep *tStep) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

    void updateDomainLinks() override;

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;

    // identification
    const char *giveClassName() const override { return "AdaptiveLinearStatic"; }
    const char *giveInputRecordName() const override { return _IFT_AdaptiveLinearStatic_Name; }

protected:
    /// Replaces the current domain by the new one (the old domain is deleted).
    void replaceDomain(Domain *dNew);
};
} // end namespace oofem
#endif // adaptlinearstatic_h
