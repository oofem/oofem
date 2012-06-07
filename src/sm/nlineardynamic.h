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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef nlineardynamic_h
#define nlineardynamic_h

#include "structengngmodel.h"
#include "sparsemtrxtype.h"
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "processcomm.h"
#endif

namespace oofem {
/**
 * This class implements nonlinear dynamic engineering problem.
 * Solution of this problem is performed  as a series of increments (loading or displacement).
 * At start of Each increment we assemble new tangent stiffness, and iteratively trying
 * to fulfill balance of external and real internal forces
 * at end of load step (see numerical method ).
 * The loading applied can bo of two types:
 * - proportional incremental  loading
 * - non-proportional fixed loading, reflecting the previous history,
 *   but could not be scaled (like dead weight).
 *
 * Tasks:
 * - Creating Numerical method for solving nonlinear problem.
 * - Assembling tangent stiffness matrix.
 * - Interfacing Numerical method to Elements.
 * - Managing time steps.
 *
 * Solution proceedure described in:
 * A SURVEY OF DIRECT TIME-INTEGRATION METHODS IN COMPUTATIONAL STRUCTURAL DYNAMICS - II. IMPLICIT METHODS
 * K. Subbaraj and M. A. Dokainish
 * Computers & Structures Vol. 32. No. 6. pp. 1387-1401, 1989
 *
 * @author Andreas Feymark
 * @author Alper Cesur
 */

enum NonLinearDynamic_ddtScheme {
    newmark = 0, ///< Newmark-beta method
    euler2  = 1, ///< Two-point Backward method
    euler3  = 2, ///< Three-point Backward Euler method
};

class NonLinearDynamic : public StructuralEngngModel
{
protected:
    SparseMtrx *stiffnessMatrix;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;

    int initFlag;
    NonLinearDynamic_ddtScheme ddtScheme;
    double dumpingCoef, alpha, beta;
    double a0, a1, a2, a3, a4, a5, a6, a7;

    FloatArray velocityVector, accelerationVector, previousLoadVector;
    FloatArray help, rhs, rhs2, previousInternalForces;
    FloatArray previousIncrementOfDisplacement;
    FloatArray previousTotalDisplacement, totalDisplacement,  incrementOfDisplacement;
    FloatArray internalForces;

    /// A load vector already applied, which does not scales.
    FloatArray initialLoadVector;
    FloatArray incrementalLoadVector;

    /// A load vector which does not scale for prescribed DOFs.
    FloatArray initialLoadVectorOfPrescribed;

    /// incremental Load Vector for prescribed DOFs.
    FloatArray incrementalLoadVectorOfPrescribed;

    int currentIterations;
    int commInitFlag;
    int nonlocalStiffnessFlag;
    NM_Status numMetStatus;
    /// Numerical method used to solve the problem.
    SparseNonLinearSystemNM *nMethod;
    /// Intrinsic time increment.
    double deltaT;

    SparseNonLinearSystemNM :: referenceLoadInputModeType refLoadInputMode;

    virtual void giveElementCharacteristicMatrix(FloatMatrix &answer, int num,
                                                 CharType type, TimeStep *tStep, Domain *domain);

public:
    NonLinearDynamic(int i, EngngModel *_master = NULL);
    virtual ~NonLinearDynamic();

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void terminate(TimeStep *tStep);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    virtual void updateYourself(TimeStep *tStep);
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn, Domain *d);
    virtual void updateAttributes(MetaStep *mStep);

    virtual double giveUnknownComponent(EquationID eid, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void updateDomainLinks();
    virtual int checkConsistency();

    // Identification
    virtual const char *giveClassName() const { return "NonLinearDynamic"; }
    virtual classType giveClassID() const { return NonLinearDynamicClass; }
    virtual int isIncremental() { return 1; }
    virtual fMode giveFormulation() { return nonLinFormulation; }
    virtual int useNonlocalStiffnessOption() { return this->nonlocalStiffnessFlag; }
    /// For load balancing purposes we store all values with same EquationID; so hash is computed from mode value only
    virtual int giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN)
    { return ( int ) mode; }
    void timesMassMtrx(FloatArray &answer, FloatArray &vec, Domain *domain, TimeStep *tStep);

#ifdef __OOFEG
    void showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime);
#endif

#ifdef __PARALLEL_MODE
    int estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType);
    /**
     * Initializes communication maps of the receiver.
     */
    void initializeCommMaps(bool forceInit = false);

    virtual LoadBalancer *giveLoadBalancer();
    virtual LoadBalancerMonitor *giveLoadBalancerMonitor();

#endif

#ifdef __PETSC_MODULE
    virtual void initPetscContexts();
#endif

protected:
    void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID ut, CharType type,
                  const UnknownNumberingScheme &, Domain *domain);

    void proceedStep(int di, TimeStep *tStep);
    void computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di);
    void assembleIncrementalReferenceLoadVectors(FloatArray &_incrementalLoadVector,
                                                 FloatArray &_incrementalLoadVectorOfPrescribed,
                                                 SparseNonLinearSystemNM :: referenceLoadInputModeType _refMode,
                                                 Domain *sourceDomain, EquationID ut, TimeStep *tStep);
#ifdef __PARALLEL_MODE
    /** Packs receiver data when rebalancing load. When rebalancing happens, the local numbering will be lost on majority of processors.
     *  Instead of identifying values of solution vectors that have to be send/received and then performing renumbering, all solution vectors
     *  are assumed to be stored in dof dictionaries before data migration. Then dofs will take care themselves for packing and unpacking. After
     *  data migration and local renumbering, the solution vectors will be restored from dof dictionary data back.
     */
    virtual void packMigratingData(TimeStep *tStep);
    /** Unpacks receiver data when rebalancing load. When rebalancing happens, the local numbering will be lost on majority of processors.
     *  Instead of identifying values of solution vectors that have to be send/received and then performing renumbering, all solution vectors
     *  are assumed to be stored in dof dictionaries before data migration. Then dofs will take care themselves for packing and unpacking. After
     *  data migration and local renumbering, the solution vectors will be restored from dof dictionary data back.
     */
    virtual void unpackMigratingData(TimeStep *tStep);
#endif
};
} // end namespace oofem
#endif // nlineardynamic_h
