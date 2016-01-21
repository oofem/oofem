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

#ifndef nlineardynamic_h
#define nlineardynamic_h

#include "../sm/EngineeringModels/structengngmodel.h"
#include "sparsemtrxtype.h"
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "timediscretizationtype.h"

///@name Input fields for NonLinearDynamic
//@{
#define _IFT_NonLinearDynamic_Name "nonlineardynamic"
#define _IFT_NonLinearDynamic_deltat "deltat"
#define _IFT_NonLinearDynamic_refloadmode "refloadmode"
#define _IFT_NonLinearDynamic_nonlocstiff "nonlocstiff"
#define _IFT_NonLinearDynamic_nonlocalext "nonlocalext"
#define _IFT_NonLinearDynamic_ddtScheme "ddtscheme"
#define _IFT_NonLinearDynamic_gamma "gamma"
#define _IFT_NonLinearDynamic_beta "beta"
#define _IFT_NonLinearDynamic_eta "eta"
#define _IFT_NonLinearDynamic_delta "delta"
//@}

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

class NonLinearDynamic : public StructuralEngngModel
{
protected:
    std :: unique_ptr< SparseMtrx > effectiveStiffnessMatrix, massMatrix;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;

    int initFlag;
    TimeDiscretizationType initialTimeDiscretization;
    double gamma, beta;
    double eta, delta;
    double a0, a1, a2, a3, a4, a5, a6, a7;

    FloatArray velocityVector, accelerationVector, previousLoadVector;
    FloatArray previousVelocityVector, previousAccelerationVector;
    FloatArray help, rhs, rhs2, previousInternalForces;
    FloatArray previousIncrementOfDisplacement;
    FloatArray previousTotalDisplacement, totalDisplacement,  incrementOfDisplacement;
    FloatArray internalForces, forcesVector;

    int currentIterations, totIterations, MANRMSteps;
    int commInitFlag;
    int nonlocalStiffnessFlag;
    /// Numerical method used to solve the problem.
    std :: unique_ptr< SparseNonLinearSystemNM > nMethod;
    /// Intrinsic time increment.
    double deltaT;

public:
    NonLinearDynamic(int i, EngngModel * _master = NULL);
    virtual ~NonLinearDynamic();

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void terminate(TimeStep *tStep);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep);

    virtual void updateYourself(TimeStep *tStep);
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn, Domain *d);
    virtual void updateAttributes(MetaStep *mStep);
	virtual void initializeYourself(TimeStep *tStep);

    virtual double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void updateDomainLinks();

    // Identification
    virtual const char *giveInputRecordName() const { return _IFT_NonLinearDynamic_Name; }
    virtual const char *giveClassName() const { return "NonLinearDynamic"; }
    virtual fMode giveFormulation() { return nonLinFormulation; }
    virtual int useNonlocalStiffnessOption() { return this->nonlocalStiffnessFlag; }
    virtual int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) { return ( int ) mode; }
    void timesMtrx(FloatArray &answer, FloatArray &vec, CharType type, Domain *domain, TimeStep *tStep);

	TimeDiscretizationType giveInitialTimeDiscretization() { return initialTimeDiscretization; }

#ifdef __OOFEG
    void showSparseMtrxStructure(int type, oofegGraphicContext &gc, TimeStep *tStep);
#endif

    virtual int estimateMaxPackSize(IntArray &commMap, DataStream &buff, int packUnpackType);
#ifdef __PARALLEL_MODE
    virtual LoadBalancer *giveLoadBalancer();
    virtual LoadBalancerMonitor *giveLoadBalancerMonitor();
#endif

protected:
    void assemble(SparseMtrx &answer, TimeStep *tStep, const MatrixAssembler &ma,
                  const UnknownNumberingScheme &, Domain *domain);

    void proceedStep(int di, TimeStep *tStep);
    void determineConstants(TimeStep *tStep);

    virtual void packMigratingData(TimeStep *tStep);
    virtual void unpackMigratingData(TimeStep *tStep);
};
} // end namespace oofem
#endif // nlineardynamic_h
