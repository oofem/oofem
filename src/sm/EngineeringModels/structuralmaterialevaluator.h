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

#ifndef structuralmaterialevaluator_h
#define structuralmaterialevaluator_h

#include "engngm.h"

#include <fstream>

///@name Input fields for structural material evaluator
//@{
#define _IFT_StructuralMaterialEvaluator_Name "structuralmaterialevaluator"
#define _IFT_StructuralMaterialEvaluator_deltat "deltat"
#define _IFT_StructuralMaterialEvaluator_numberOfTimeSteps "nsteps"
#define _IFT_StructuralMaterialEvaluator_componentFunctions "componentfunctions" ///< Integer list of time functions for each component
#define _IFT_StructuralMaterialEvaluator_stressControl "stresscontrol" ///< Integer list of the stress components which are controlled
#define _IFT_StructuralMaterialEvaluator_outputVariables "vars" ///< Variables (from integration point) to be written.
#define _IFT_StructuralMaterialEvaluator_tolerance "tolerance" ///< Tolerance for stress control
//@}

namespace oofem {
class GaussPoint;

/**
 * For testing material behavior, particularly useful for multiscale modeling where one can test a single RVE.
 * The deviatoric and volumetric parts are split. No nodes or elements are used.
 *
 * This model will output data in its own way since it does not contain any actual FE-results so no export modules are called.
 * @author Mikael Öhman
 */
class StructuralMaterialEvaluator : public EngngModel
{
protected:
    double deltaT; ///< Time increment.

    IntArray cmpntFunctions; /// Time functions controlling each component of the deviatoric part of the stress.
    IntArray sControl, eControl;
    IntArray vars;

    std::vector< std :: unique_ptr< GaussPoint > >gps;

    std :: ofstream outfile;

    double tolerance;

public:
    StructuralMaterialEvaluator(int i, EngngModel * _master = NULL);
    virtual ~StructuralMaterialEvaluator();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void solveYourself();

    virtual int checkConsistency();
    virtual void doStepOutput(TimeStep *tStep);
    virtual TimeStep *giveNextStep();

    virtual const char *giveClassName() const { return "StructuralMaterialEvaluator"; }
    virtual const char *giveInputRecordName() const { return _IFT_StructuralMaterialEvaluator_Name; }
};
} // end namespace oofem

#endif // structuralmaterialevaluator_h
