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
    virtual void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const;
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
    virtual void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const;
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

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;

    std :: unique_ptr< SparseMtrx > lhs;
    std :: unique_ptr< PrimaryField > VelocityPressureField;
    //PrimaryField VelocityField;
    FloatArray accelerationVector; //, previousAccelerationVector;
    FloatArray incrementalSolutionVector;

    FloatArray internalForces;
    FloatArray eNorm;

    ///@todo Use ScalarFunction here!
    double deltaT;
    int deltaTF;
    /// Convergence tolerance.
    double atolv, rtolv;
    /// Max number of iterations.
    int maxiter;
    /**
     * Flag if set to true (default), then when max number of iteration reached, computation stops
     * otherwise computation continues with next step.
     */
    bool stopmaxiter;
    /// Integration constant.
    double alpha;

    int initFlag;
    int consistentMassFlag;

    bool equationScalingFlag;
    // length, velocity, and density scales
    double lscale, uscale, dscale;
    /// Reynolds number.
    double Re;

    // material interface representation for multicomponent flows
    std :: unique_ptr< MaterialInterface > materialInterface;
    // map of active dofmans for problems with free surface and only one fluid considered
    // IntArray __DofManActivityMask;
    // free surface flag -> we solve free surface problem by single reference fluid
    // int fsflag;

public:
    SUPG(int i, EngngModel * _master = NULL);
    virtual ~SUPG();

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);

    virtual double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof);
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    virtual double giveReynoldsNumber();

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void updateDomainLinks();

    virtual TimeStep *giveNextStep();
    virtual TimeStep *giveSolutionStepWhenIcApply(bool force = false);
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int checkConsistency();
    // identification
    virtual const char *giveClassName() const { return "SUPG"; }
    virtual const char *giveInputRecordName() const { return _IFT_SUPG_Name; }

    virtual fMode giveFormulation() { return TL; }

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep);

    virtual int requiresUnknownsDictionaryUpdate() { return renumberFlag; }

    virtual bool giveEquationScalingFlag() { return equationScalingFlag; }
    virtual double giveVariableScale(VarScaleType varId);

    virtual void updateDofUnknownsDictionary(DofManager *dman, TimeStep *tStep);
    virtual int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep);

    virtual MaterialInterface *giveMaterialInterface(int n) { return materialInterface.get(); }


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
