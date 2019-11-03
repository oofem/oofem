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

#ifndef adaptlinearstatic_h
#define adaptlinearstatic_h

#include "sm/EngineeringModels/linearstatic.h"
#include "meshpackagetype.h"

///@name Input fields for AdaptiveLinearStatic
//@{
#define _IFT_AdaptiveLinearStatic_Name "adaptlinearstatic"
#define _IFT_AdaptiveLinearStatic_meshpackage "meshpackage"
//@}

namespace oofem {
/**
 * This class implements an adaptive linear static engineering problem.
 * Multiple loading cases are not supported.
 * Due to linearity of a problem, the complete reanalysis from the beginning
 * is done after adaptive remeshing.
 * Solution steps represent a series of adaptive analyses.
 */
class AdaptiveLinearStatic : public LinearStatic
{
protected:
    /// Meshing package used for refinements.
    MeshPackageType meshPackage;

public:
    AdaptiveLinearStatic(int i, EngngModel *master = nullptr) : LinearStatic(i, master) { }
    virtual ~AdaptiveLinearStatic() { }

    void updateYourself(TimeStep *tStep) override;

    /**
     * Initializes the newly generated discretization state according to previous solution.
     * This process should typically include restoring old solution, instanciating newly
     * generated domain(s) and by mapping procedure.
     */
    int initializeAdaptive(int tStepNumber) override;
    void printOutputAt(FILE *file, TimeStep *tStep) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

    void updateDomainLinks() override;

    void initializeFrom(InputRecord &ir) override;

    // identification
    const char *giveClassName() const override { return "AdaptiveLinearStatic"; }
    const char *giveInputRecordName() const override { return _IFT_AdaptiveLinearStatic_Name; }
};
} // end namespace oofem
#endif // adaptlinearstatic_h
