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

#ifndef mpmproblem_h
#define mpmproblem_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "function.h"
#include "dofdistributedprimaryfield.h"

///@name Input fields for mpmproblem
//@{
#define _IFT_MPMProblem_Name "mpmproblem"
#define _IFT_MPMProblem_initt "initt"
#define _IFT_MPMProblem_deltat "deltat"
#define _IFT_MPMProblem_deltatfunction "deltatfunction"
#define _IFT_MPMProblem_prescribedtimes "prescribedtimes"
#define _IFT_MPMProblem_alpha "alpha"
#define _IFT_MPMProblem_keepTangent "keeptangent" ///< Fixes the tangent to be reused on each step.
#define _IFT_MPMProblem_exportFields "exportfields" ///< Fields to export for staggered problems.
#define _IFT_MPMProblem_problemType "ptype" 
//@}

namespace oofem {

/**
 * Callback class for assembling mid point effective tangents. 
 * @todo Need to parametrize individual contributing terms, ther locations and multilication factors.
 */
class UPLhsAssembler : public MatrixAssembler
{
protected:
    double alpha;
    double deltaT;

public:
    UPLhsAssembler(double alpha, double deltaT);
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};

/**
 * Callback class for assembling residuals
 */
class UPResidualAssembler : public VectorAssembler
{
    protected:
    double alpha;
    double deltaT;
public:
    UPResidualAssembler(double alpha, double deltaT) : VectorAssembler(), alpha(alpha), deltaT(deltaT) {}
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
};

/**
 * Callback class for assembling mid point effective tangents. 
 * @todo Need to parametrize individual contributing terms, ther locations and multilication factors.
 */
class TMLhsAssembler : public MatrixAssembler
{
protected:
    double alpha;
    double deltaT;

public:
    TMLhsAssembler(double alpha, double deltaT);
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};

/**
 * Callback class for assembling residuals
 */
class TMResidualAssembler : public VectorAssembler
{
    protected:
    double alpha;
    double deltaT;
public:
    TMResidualAssembler(double alpha, double deltaT) : VectorAssembler(), alpha(alpha), deltaT(deltaT) {}
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
};



/**
 * This class represents generic nonlinear multi-physics problem. The problem can be growing/decreasing, signalized by flag "changingproblemsize"
 * in the problem description. The solution is stored in UnknownsField, which can obtain/ project solution from/to DOFs (nodes). If the problem
 * keeps the same equation numbers, solution is taken from UnknownsField without any projection, which is more efficient. See the matlibmanual
 * for solution strategy of balance equations and the solution algorithm.
 *
 * @todo Documentation errors (there is no "UnknownsField" used here).
 */
class MPMProblem : public EngngModel
{
protected:
    SparseMtrxType sparseMtrxType = SMT_Skyline;
    std :: unique_ptr< DofDistributedPrimaryField > field;

    std :: unique_ptr< SparseMtrx > effectiveMatrix;

    FloatArray solution;
    FloatArray internalForces;
    FloatArray eNorm;

    /// Numerical method used to solve the problem
    std :: unique_ptr< SparseNonLinearSystemNM > nMethod;

    /// Initial time from which the computation runs. Default is zero.
    double initT = 0.;
    /// Length of time step.
    double deltaT = 0.;
    double alpha = 0.;
    /// Associated time function for time step increment.
    int dtFunction = 0;
    /// Specified times where the problem is solved
    FloatArray prescribedTimes;
    bool keepTangent = false, hasTangent = false;
    IntArray exportFields;
    /// identifies what problem to solve (UP, UPV, etc) 
    std::string problemType; 

public:
    MPMProblem(int i, EngngModel * _master);


    void solveYourselfAt(TimeStep *tStep) override;
    void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d) override;
    bool newDofHandling() override { return true; }
    void updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d) override;
    void updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm) override;
    void updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d) override;
    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    virtual void applyIC();

    int requiresUnknownsDictionaryUpdate() override;
    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override;
    void updateDomainLinks() override;

    Function *giveDtFunction();
    double giveDeltaT(int n);
    double giveDiscreteTime(int iStep);

    TimeStep *giveNextStep() override;
    TimeStep *giveSolutionStepWhenIcApply(bool force = false) override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;

    void initializeFrom(InputRecord &ir) override;

    bool requiresEquationRenumbering(TimeStep *tStep) override;
    int forceEquationNumbering() override;

    void printOutputAt(FILE *file, TimeStep *tStep) override;
    
    void updateYourself(TimeStep *tStep) override;
    
    int checkConsistency() override;
    FieldPtr giveField (FieldType key, TimeStep *tStep) override;
    // identification
    const char *giveInputRecordName() const { return _IFT_MPMProblem_Name; }
    const char *giveClassName() const override { return "MPMProblem"; }
    fMode giveFormulation() override { return TL; }

  /** nlinear statics number starts simulation at time = 0
   */
  double giveFinalTime() //override
  {
    if(prescribedTimes.giveSize()) {
      return prescribedTimes.at(prescribedTimes.giveSize());
    } else {
      return deltaT * numberOfSteps;
    }
  }

};

  /*
    TimeStep* giveNextStep() override;
    void solveYourselfAt(TimeStep *tStep) override;
    void updateYourself(TimeStep *tStep) override;
    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override;

    void initializeFrom(InputRecord &ir) override;

    // identification
    const char *giveInputRecordName() const { return _IFT_MPMProblem_Name; }
    const char *giveClassName() const override { return "MPMProblem"; }
    fMode giveFormulation() override { return nonLinFormulation; }
    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override;
    void updateDofUnknownsDictionary(DofManager *dman, TimeStep *tStep) override;

    int giveCurrentNumberOfIterations() override { return currentIterations; }
    TimeStep *giveSolutionStepWhenIcApply(bool force = false) override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;


protected:
    void updateInternalState(TimeStep *tStep) ;
    void applyIC(TimeStep *tStep) ;
    void createPreviousSolutionInDofUnknownsDictionary(TimeStep *tStep);
    double giveDeltaT(int n);
    virtual void copyUnknownsInDictionary(ValueModeType mode, TimeStep *fromTime, TimeStep *toTime);
    Function * giveDtFunction();
    double giveDiscreteTime(int n);
};
*/
} // end namespace oofem
#endif // mpmproblem_h
