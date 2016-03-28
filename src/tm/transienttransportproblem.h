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

#ifndef transienttransportproblem_h
#define transienttransportproblem_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"

#include <memory>

///@name Input fields for TransientTransportProblem
//@{
#define _IFT_TransientTransportProblem_Name "transienttransport"
#define _IFT_TransientTransportProblem_alpha "alpha" ///< Defines the time discretization of the value: @f$ u = u_n (1-\alpha) + u_{n+1} \alpha @f$.
#define _IFT_TransientTransportProblem_deltaT "deltat" ///< Fixed timestep.
#define _IFT_TransientTransportProblem_dtFunction "dtfunction" ///< Function that determines size of time step.
#define _IFT_TransientTransportProblem_prescribedTimes "prescribedtimes" ///< Discrete times for each time step.
#define _IFT_TransientTransportProblem_keepTangent "keeptangent" ///< Fixes the tangent to be reused on each step.
#define _IFT_TransientTransportProblem_lumped "lumped" ///< Use of lumped "mass" matrix
#define _IFT_TransientTransportProblem_exportFields "exportfields" ///< Fields to export for staggered problems.
//@}

namespace oofem {
class SparseNonLinearSystemNM;
class PrimaryField;
class Function;

/**
 * Solves general nonlinear transient transport problems.
 * @author Mikael Öhman
 */
class TransientTransportProblem : public EngngModel
{
protected:
    SparseMtrxType sparseMtrxType;
    std :: unique_ptr< PrimaryField > field;

    std :: unique_ptr< SparseMtrx > capacityMatrix;
    FloatArray capacityDiag; /// In case of a lumped matrix, the diagonal entries are stored here.
    std :: unique_ptr< SparseMtrx > effectiveMatrix;

    FloatArray solution;
    FloatArray internalForces;
    FloatArray eNorm;

    /// Numerical method used to solve the problem
    std :: unique_ptr< SparseNonLinearSystemNM > nMethod;

    double alpha;
    int dtFunction;
    FloatArray prescribedTimes;
    double deltaT;
    bool keepTangent;
    bool lumped;

    IntArray exportFields;

public:
    /// Constructor.
    TransientTransportProblem(int i, EngngModel * _master);
    /// Destructor.
    virtual ~TransientTransportProblem();

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    virtual double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    
    virtual void applyIC();

    virtual int requiresUnknownsDictionaryUpdate();
    virtual int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep);
    virtual void updateDomainLinks();

    Function *giveDtFunction();
    double giveDeltaT(int n);
    double giveDiscreteTime(int iStep);

    virtual TimeStep *giveNextStep();
    virtual TimeStep *giveSolutionStepWhenIcApply(bool force = false);
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual bool requiresEquationRenumbering(TimeStep *tStep);
    virtual int forceEquationNumbering();

    virtual int checkConsistency();

    // identification
    virtual const char *giveInputRecordName() const { return _IFT_TransientTransportProblem_Name; }
    virtual const char *giveClassName() const { return "TransientTransportProblem"; }
    virtual fMode giveFormulation() { return TL; }
};
} // end namespace oofem
#endif // transienttransportproblem_h
