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

#ifndef stationarytransportproblem_h
#define stationarytransportproblem_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"

///@name Input fields for StationaryTransportProblem
//@{
#define _IFT_StationaryTransportProblem_Name "stationaryproblem"
#define _IFT_StationaryTransportProblem_exportfields "exportfields"
#define _IFT_StationaryTransportProblem_keepTangent "keeptangent"
//@}

namespace oofem {
class SparseNonLinearSystemNM;

/**
 * This class represents stationary transport problem.
 */
class StationaryTransportProblem : public EngngModel
{
protected:
    SparseMtrxType sparseMtrxType;
    /// This field stores solution vector. For fixed size of problem, the PrimaryField is used, for growing/decreasing size, DofDistributedPrimaryField applies.
    PrimaryField *UnknownsField;

    SparseMtrx *conductivityMatrix;
    FloatArray internalForces;
    FloatArray eNorm;

    /// Numerical method used to solve the problem
    SparseNonLinearSystemNM *nMethod;

    bool keepTangent;

public:
    /// Constructor.
    StationaryTransportProblem(int i, EngngModel *_master);
    /// Destructor.
    virtual ~StationaryTransportProblem();

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    virtual double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    void assembleDirichletBcRhsVector(FloatArray &answer, TimeStep *tStep, EquationID ut,
                                                           ValueModeType mode, CharType lhsType,
                                                           const UnknownNumberingScheme &ns, Domain *d);
    virtual void updateDomainLinks();

    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int checkConsistency();

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    // identification
    virtual const char *giveClassName() const { return "StationaryTransportProblem"; }
    virtual classType giveClassID() const { return StationaryTransportProblemClass; }
    virtual fMode giveFormulation() { return TL; }

#ifdef __PETSC_MODULE
    virtual void initPetscContexts();
#endif

protected:
    /**
     * Updates IP values on elements
     * @param stepN Solution step.
     */
    virtual void updateInternalState(TimeStep *stepN);
};
} // end namespace oofem
#endif // stationarytransportproblem_h
