/*
 * DarcyFlow.h
 *
 *  Created on: Feb 25, 2010
 *      Author: carl
 */



#ifndef DARCYFLOW_H_
#define DARCYFLOW_H_

#include "engngm.h"
#include "inputrecord.h"
#include "sparsemtrxtype.h"
#include "linsystsolvertype.h"
#include "sparselinsystemnm.h"
// #include "LinearFlowBase.h"
#include "sparsenonlinsystemnm.h"

#ifdef __PARALLEL_MODE
#include "problemcomm.h"
#include "processcomm.h"
#endif

namespace oofem {

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
    PrimaryField *PressureField;
    SparseMtrxType sparseMtrxType;
    SparseNonLinearSystemNM *nMethod;


    SparseMtrx *stiffnessMatrix;
    FloatArray internalForces;
    FloatArray externalForces;
    FloatArray incrementOfSolution;

    int currentIterations;
    double ebenorm, loadLevel;
    FloatArray incrementOfDisplacement;
    FloatArray incrementalBCLoadVector;
    FloatArray incrementalLoadVector;
    FloatArray initialLoad;
    SparseNonLinearSystemNM :: referenceLoadInputModeType refLoadInputMode;
    bool hasAdvanced;

public:
    DarcyFlow(int i, EngngModel *_master);
    virtual ~DarcyFlow();

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL) {return CIO_OK; };
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL) {return CIO_OK; };
    virtual void updateDomainLinks() { return; };
    virtual int checkConsistency() { return 1; };
    virtual fMode giveFormulation() { return TL; }
    virtual void updateInternalState(TimeStep *stepN) {return; };
    virtual void giveInternalForces(FloatArray &answer, const FloatArray &DeltaR, Domain *domain, TimeStep *stepN) {};

    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);

    virtual double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);
    void DumpMatricesToFile(FloatMatrix *LHS, FloatArray *RHS, FloatArray *SolutionVector);

    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);
    virtual TimeStep * giveNextStep();

    virtual int forceEquationNumbering(int id);

    virtual const char *giveClassName() const { return "DarcyFlowEngngModel"; }
    virtual classType giveClassID() const { return DarcyFlowClass; }

#ifdef __PARALLEL_MODE
    CommunicatorBuff *commBuff; //new CommunicatorBuff(this->giveNumberOfProcesses(), CBT_static);
    ProblemCommunicator *communicator;
#endif

#ifdef __PETSC_MODULE
    virtual void initPetscContexts();
#endif

};

}

#endif /* DARCYFLOW_H_ */
