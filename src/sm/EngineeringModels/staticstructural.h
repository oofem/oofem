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

#ifndef staticstructural_h
#define staticstructural_h

#include "sm/EngineeringModels/structengngmodel.h"
#include "sm/EngineeringModels/xfemsolverinterface.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrxtype.h"
#include "staggeredsolver.h"

#define _IFT_StaticStructural_Name "staticstructural"
#define _IFT_StaticStructural_deltat "deltat"
#define _IFT_StaticStructural_prescribedTimes "prescribedtimes"
#define _IFT_StaticStructural_solvertype "solvertype"
#define _IFT_StaticStructural_stiffmode "stiffmode"
#define _IFT_StaticStructural_nonlocalExtension "nonlocalext"

#define _IFT_StaticStructural_recomputeaftercrackpropagation "recomputeaftercrackprop"
namespace oofem {
class SparseMtrx;
class DofDistributedPrimaryField;

/**
 * Solves a static structural problem.
 * @author Mikael Öhman
 */
class StaticStructural : public StructuralEngngModel, public XfemSolverInterface
{
protected:
    FloatArray solution;
    FloatArray internalForces;
    FloatArray eNorm;
    FloatArray referenceForces, externalForces;

public:
    std :: unique_ptr< SparseMtrx >stiffnessMatrix;
protected:

    std :: unique_ptr< DofDistributedPrimaryField >field;

    SparseMtrxType sparseMtrxType;

    std :: unique_ptr< SparseNonLinearSystemNM >nMethod;
    std :: string solverType;
    MatResponseMode stiffMode;

    double loadLevel;
    double deltaT;
    FloatArray prescribedTimes;

    InitialGuess initialGuessType;

    bool mRecomputeStepAfterPropagation;

public:
    StaticStructural(int i, EngngModel *master=nullptr);
    virtual ~StaticStructural();
    void initializeFrom(InputRecord &ir) override;
    void updateAttributes(MetaStep *mStep) override;

    void solveYourself() override;
    void solveYourselfAt(TimeStep *tStep) override;

    void terminate(TimeStep *tStep) override;

    void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d) override;
    void updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d) override;
    void updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm) override;
    void updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d) override;

    double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof) override;
    bool newDofHandling() override { return true; }

    void updateDomainLinks() override;

    int forceEquationNumbering() override;

    double giveLoadLevel() override { return loadLevel; }
    TimeStep *giveNextStep() override;
    double giveEndOfTimeOfInterest() override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;

    //fMode giveFormulation() override { return TL; }

    bool requiresEquationRenumbering(TimeStep *tStep) override;

    int requiresUnknownsDictionaryUpdate() override { return true; }
    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override;
    void computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di) override;
    // identification
    const char *giveInputRecordName() const { return _IFT_StaticStructural_Name; }
    const char *giveClassName() const override { return "StaticStructural"; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    int estimateMaxPackSize(IntArray &commMap, DataStream &buff, int packUnpackType) override;
};
} // end namespace oofem
#endif // staticstructural_h
