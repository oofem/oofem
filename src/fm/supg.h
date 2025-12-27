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

#ifndef supg_h
#define supg_h

#include "fluidmodel.h"
#include "sparsemtrxtype.h"
#include "sparselinsystemnm.h"
#include "primaryfield.h"
#include "materialinterface.h"
#include "assemblercallback.h"

#include <memory>

///@name Input fields for SUPG
//@{
#define _IFT_SUPG_Name "supg"
#define _IFT_SUPG_deltat "deltat"
#define _IFT_SUPG_deltatFunction "deltatltf"
#define _IFT_SUPG_cmflag "cmflag"
#define _IFT_SUPG_alpha "alpha"
#define _IFT_SUPG_scaleflag "scaleflag"
#define _IFT_SUPG_lscale "lscale"
#define _IFT_SUPG_uscale "uscale"
#define _IFT_SUPG_dscale "dscale"
#define _IFT_SUPG_miflag "miflag"
#define _IFT_SUPG_rtolv "rtolv"
#define _IFT_SUPG_atolv "atolv"
#define _IFT_SUPG_maxiter "maxiter"
#define _IFT_SUPG_stopmaxiter "stopmaxiter"
#define _IFT_SUPG_fsflag "fsflag"
//@}

namespace oofem {
class SparseMtrx;
class SparseNonLinearSystemNM;

/**
 * Callback class for assembling SUPG internal forces
 * @author Mikael Öhman
 */
class SUPGInternalForceAssembler : public VectorAssembler
{
protected:
    double lscale, dscale, uscale;

public:
    SUPGInternalForceAssembler(double l, double d, double u);
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
};

/**
 * Callback class for assembling SUPG tangent matrices
 * @author Mikael Öhman
 */
class SUPGTangentAssembler : public MatrixAssembler
{
protected:
    MatResponseMode rmode;
    double lscale, dscale, uscale;
    double alpha;

public:
    SUPGTangentAssembler(MatResponseMode m, double l, double d, double u, double a);
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};

/**
 * This class represents transient incompressible flow problem. Solution is based on
 * algorithm with SUPG/PSPG stabilization.
 */
class SUPG : public FluidModel
{
protected:
    /// Numerical method used to solve the problem
    std :: unique_ptr< SparseLinearSystemNM >nMethod;

    LinSystSolverType solverType=LinSystSolverType::ST_Direct;
    SparseMtrxType sparseMtrxType=SparseMtrxType::SMT_Skyline;

    std :: unique_ptr< SparseMtrx > lhs;
    std :: unique_ptr< PrimaryField > VelocityPressureField;
    //PrimaryField VelocityField;
    FloatArray accelerationVector; //, previousAccelerationVector;
    FloatArray incrementalSolutionVector;

    FloatArray internalForces;
    FloatArray eNorm;

    ///@todo Use ScalarFunction here!
    double deltaT = 1.;
    int deltaTF = 0;
    /// Convergence tolerance.
    double atolv = 1., rtolv = 1.;
    /// Max number of iterations.
    int maxiter = 0;
    /**
     * Flag if set to true (default), then when max number of iteration reached, computation stops
     * otherwise computation continues with next step.
     */
    bool stopmaxiter = true;
    /// Integration constant.
    double alpha = 0.;

    int initFlag = 1;
    int consistentMassFlag = 0;

    bool equationScalingFlag = false;
    // length, velocity, and density scales
    double lscale = 1., uscale = 1., dscale = 1.;
    /// Reynolds number.
    double Re = 0.;

    // material interface representation for multicomponent flows
    std :: unique_ptr< MaterialInterface > materialInterface;
    // map of active dofmans for problems with free surface and only one fluid considered
    // IntArray __DofManActivityMask;
    // free surface flag -> we solve free surface problem by single reference fluid
    // int fsflag;

public:
    SUPG(int i, EngngModel * _master = nullptr);

    void solveYourselfAt(TimeStep *tStep) override;
    void updateYourself(TimeStep *tStep) override;

    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override;
    void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d) override;
    void updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d) override;
    void updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm) override;
    void updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d) override;
    double giveReynoldsNumber() override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    void updateDomainLinks() override;

    TimeStep *giveNextStep() override;
    TimeStep *giveSolutionStepWhenIcApply(bool force = false) override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;

    void initializeFrom(InputRecord &ir) override;

    int checkConsistency() override;
    // identification
    const char *giveClassName() const override { return "SUPG"; }
    const char *giveInputRecordName() const { return _IFT_SUPG_Name; }

    fMode giveFormulation() override { return TL; }

    void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep) override;

    int requiresUnknownsDictionaryUpdate() override { return renumberFlag; }

    bool giveEquationScalingFlag() override { return equationScalingFlag; }
    double giveVariableScale(VarScaleType varId) override;

    void updateDofUnknownsDictionary(DofManager *dman, TimeStep *tStep) override;
    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override;

    MaterialInterface *giveMaterialInterface(int n) override { return materialInterface.get(); }

protected:
    void updateInternalState(TimeStep *tStep);
    void applyIC(TimeStep *tStep);
    void evaluateElementStabilizationCoeffs(TimeStep *tStep);
    void updateElementsForNewInterfacePosition(TimeStep *tStep);

    void updateDofUnknownsDictionary_predictor(TimeStep *tStep);
    void updateDofUnknownsDictionary_corrector(TimeStep *tStep);

    void updateSolutionVectors(FloatArray &solutionVector, FloatArray &accelerationVector, FloatArray &incrementalSolutionVector, TimeStep *tStep);
    void updateSolutionVectors_predictor(FloatArray &solutionVector, FloatArray &accelerationVector, TimeStep *tStep);

    //void initDofManActivityMap ();
    //void updateDofManActivityMap (TimeStep* tStep);
    void updateDofManVals(TimeStep *tStep);
    //void imposeAmbientPressureInOuterNodes(SparseMtrx* lhs, FloatArray* rhs, TimeStep* tStep);
    //void __debug(TimeStep* tStep);
};
} // end namespace oofem
#endif // supg_h
