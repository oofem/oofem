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

#ifndef fluidmaterialevaluator_h
#define fluidmaterialevaluator_h

#include "engngm.h"

#include <fstream>

///@name Input fields for material evaluator
//@{
#define _IFT_FluidMaterialEvaluator_Name "fluidmaterialevaluator"
#define _IFT_FluidMaterialEvaluator_deltat "deltat"
#define _IFT_FluidMaterialEvaluator_numberOfTimeSteps "nsteps"
#define _IFT_FluidMaterialEvaluator_nDimensions "ndim" ///< Number of dimensions (2 or 3)
#define _IFT_FluidMaterialEvaluator_componentFunctions "componentfunctions" ///< Integer list of time functions for each component
#define _IFT_FluidMaterialEvaluator_stressControl "stresscontrol" ///< Integer list of the stress components which are controlled
#define _IFT_FluidMaterialEvaluator_volFunction "volfunction" ///< Integer of time function for volumetric part
#define _IFT_FluidMaterialEvaluator_pressureControl "pressurecontrol" ///< Bool(Integer) determining if pressure or volumetric strain-rate is controlled.
#define _IFT_FluidMaterialEvaluator_outputVariables "vars" ///< Variables (from integration point) to be written.
//@}

namespace oofem {
class GaussPoint;

/**
 * For testing material behavior, particularly useful for multiscale modeling where one can test a single RVE.
 * The deviatoric and volumetric parts are split. No nodes or elements are used.
 *
 * @note The user must take care to ensure that a deviatoric stress is possible, i.e. is all components are controlled for the stress,
 * the time functions must evaluate to something deviatoric (e.g. all zeros).
 *
 * This model will output data in its own way since it does not contain any actual FE-results.
 * It will not call the output manager at all.
 * @author Mikael Öhman
 */
class FluidMaterialEvaluator : public EngngModel
{
protected:
    double deltaT; ///< Time increment.

    int ndim; /// Number of spatial dimensions.
    IntArray cmpntFunctions; /// Time functions controlling each component of the deviatoric part of the stress.
    int volFunction; /// Time function controlling the volumetric/pressure part
    IntArray sControl, eControl;
    bool pressureControl;

    IntArray vars;

    std::vector< std :: unique_ptr< GaussPoint > >gps;

    std :: ofstream outfile;

public:
    FluidMaterialEvaluator(int i, EngngModel * _master = NULL);
    virtual ~FluidMaterialEvaluator();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void solveYourself();

    virtual int checkConsistency();
    virtual void doStepOutput(TimeStep *tStep);
    virtual TimeStep *giveNextStep();

    virtual const char *giveClassName() const { return "FluidMaterialEvaluator"; }
    virtual const char *giveInputRecordName() const { return _IFT_FluidMaterialEvaluator_Name; }
};
} // end namespace oofem

#endif // fluidmaterialevaluator_h
