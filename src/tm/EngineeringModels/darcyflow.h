/*
 * DarcyFlow.h
 *
 *  Created on: Feb 25, 2010
 *      Author: carl
 */

#ifndef darcyflow_h
#define darcyflow_h

#include "engngm.h"
#include "inputrecord.h"
#include "sparsemtrxtype.h"
#include "linsystsolvertype.h"
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "primaryfield.h"

#define _IFT_DarcyFlow_Name "darcyflow"

namespace oofem {
class CommunicatorBuff;
class ProblemCommunicator;

/**
 * Class describing an extended Darcy flow. A Darcy flow is a linear relation between the seepage velocity and the pressure gradient. By 'extended',
 * we imply that there is a possibility for a nonlinear relation.
 * @author Carl Sandström
 */
class DarcyFlow : public EngngModel
{
private:
    LinSystSolverType solverType;

protected:
    std :: unique_ptr< PrimaryField > PressureField;
    SparseMtrxType sparseMtrxType;
    std :: unique_ptr< SparseNonLinearSystemNM > nMethod;

    std :: unique_ptr< SparseMtrx > stiffnessMatrix;
    FloatArray internalForces;
    FloatArray externalForces;
    FloatArray incrementOfSolution;

    int currentIterations;
    FloatArray ebeNorm;
    FloatArray incrementOfDisplacement;
    FloatArray incrementalBCLoadVector;
    FloatArray incrementalLoadVector;
    FloatArray initialLoad;
    SparseNonLinearSystemNM :: referenceLoadInputModeType refLoadInputMode;
    bool hasAdvanced;

public:
    DarcyFlow(int i, EngngModel * _master);
    virtual ~DarcyFlow();

    void solveYourselfAt(TimeStep *tStep) override;
    void updateYourself(TimeStep *tStep) override;

    contextIOResultType saveContext(DataStream &stream, ContextMode mode) override { return CIO_IOERR; }
    contextIOResultType restoreContext(DataStream &stream, ContextMode mode) override { return CIO_IOERR; }
    int checkConsistency() override { return 1; }
    fMode giveFormulation() override { return TL; }

    void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d) override;

    double giveUnknownComponent(ValueModeType, TimeStep *tSteo, Domain *d, Dof *dof) override;

    IRResultType initializeFrom(InputRecord *ir) override;
    void DumpMatricesToFile(FloatMatrix *LHS, FloatArray *RHS, FloatArray *SolutionVector);

    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;
    TimeStep *giveNextStep() override;

    int forceEquationNumbering(int id) override;

    const char *giveInputRecordName() const { return _IFT_DarcyFlow_Name; }
    const char *giveClassName() const override { return "DarcyFlow"; }
};
}

#endif // darcyflow_h
