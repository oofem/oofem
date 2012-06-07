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

#ifndef primaryfield_h
#define primaryfield_h

#include "field.h"
#include "interface.h"

#include "flotarry.h"
#include "equationid.h"
#include "valuemodetype.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#ifndef __MAKEDEPEND
 #include <vector>
#endif

namespace oofem {
class PrimaryField;
class Dof;

/**
 * Element interface class. Declares the functionality required to support PrimaryField element interpolation.
 */
class EIPrimaryFieldInterface : public Interface
{
public:

    EIPrimaryFieldInterface() : Interface() { }
    /**
     * @name The element interface required by EIPrimaryFieldInterface
     */
    //@{
    /**
     * Evaluates the value of field at given point of interest (should be located inside receiver's volume) using
     * element interpolation.
     * @param answer Field evaluated at coordinate.
     * @param pf Field to use for evaluation.
     * @param coords Coordinate.
     * @param dofId IDs of DOFs to evaluate.
     * @param mode Mode of field.
     * @param atTime Time step to evaluate at.
     * @return Zero if ok, nonzero when error encountered.
     */
    virtual int EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
        FloatArray &coords, IntArray &dofId, ValueModeType mode, TimeStep *atTime) = 0;
    //@}
};


/**
 * Abstract class representing field of primary variables (those, which are unknown and are typically
 * associated to nodes).
 * In the current design the primary field is understood as simple database, that allows to keep track
 * of the history of a solution vector representing primary field. The history is kept as sequence of
 * solution vectors. The history depth kept can be selected. PrimaryField class basically provides access to
 * time-dependent vectors of the field of unknowns. It adds the possibility to
 * further interpolate the field values using element interpolation functions.
 * The prescribed values of the field are not maintained, since they can be obtained directly from
 * corresponding DOFs of associated domain.
 *
 * As the PrimaryField stores the state directly in solution vectors that are usually directly
 * updated by EngngModel, it may contain a mix of different fields (this is especially true for
 * strongly coupled problems). Then masked primary field can be used to select only certain DOFs
 * (based on DofID) from its master PrimaryField.
 */
class PrimaryField : public Field
{
protected:
    int actualStepNumber;
    int actualStepIndx;
    int nHistVectors;
    AList< FloatArray >solutionVectors;
    AList< TimeStep >solStepList;
    EngngModel *emodel;
    int domainIndx;
    EquationID ut;

public:
    /**
     * Constructor. Creates a field of given type associated to given domain.
     * Not using pointer to domain, because this will prevent the use of PrimaryField as an
     * EngngModel attribute. This is because the domain does not exists when
     * PrimaryField is created (this is when EngngModel is created).
     * @param a Engineering model which field belongs to.
     * @param idomain Index of domain for field.
     * @param ft Type of stored field.
     * @param ut Equation ID for unknowns in field.
     * @param nHist Number of old time steps to store.
     */
    PrimaryField(EngngModel *a, int idomain, FieldType ft, EquationID ut, int nHist);
    virtual ~PrimaryField();

    /**
     * Copy unknowns from previous solution or DOF's dictionary to the solution vector
     * @param mode what the unknown describes (increment, total value etc.).
     * @param atTime Time of interest.
     * @param answer Resulting vector.
     */
    virtual void initialize(ValueModeType mode, TimeStep *atTime, FloatArray &answer);

    /**
     * @param dof Pointer to DOF.
     * @param mode What the unknown describes (increment, total value etc.).
     * @param atTime Time step of interest.
     * @return Value of interest at given DOF.
     */
    virtual double giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *atTime);
    /**
     * Evaluates the field at given point
     * @param answer Evaluated field at point.
     * @param coords Coordinates of the point of interest.
     * @param mode Mode of evaluated unknowns.
     * @param atTime Time step of interest.
     * @return Error code (0-ok, 1-point not found in domain).
     */
    virtual int evaluateAt(FloatArray &answer, FloatArray &coords,
                           ValueModeType mode, TimeStep *atTime);
    /**
     * Evaluates the field at given DofManager
     * @param answer Evaluated field at dman.
     * @param dman DOF manager to evaluate at.
     * @param mode Mode of evaluated unknowns.
     * @param atTime Time step of interest.
     * @return Error code (0-ok, 1-point not found in domain).
     */
    virtual int evaluateAt(FloatArray &answer, DofManager* dman,
                           ValueModeType mode, TimeStep *atTime);
    /**
     * Evaluates the field at given DOF manager, allows to select specific
     * dofs using mask
     * @param answer Evaluated field at dman.
     * @param dman DOF manager of interest.
     * @param mode Mode of evaluated unknowns.
     * @param atTime Time step of interest.
     * @param dofId Dof mask, id set to NULL, all Dofs evaluated.
     * @return Error code (0=ok, 1=point not found in domain).
     */
    virtual int __evaluateAt(FloatArray &answer, DofManager* dman,
                    ValueModeType mode, TimeStep *atTime, IntArray *dofId);
    /**
     * Evaluates the field at given point, allows to select specific
     * dofs using mask.
     * @param answer Evaluated field at coords.
     * @param coords Coordinates of the point of interest.
     * @param mode Mode of evaluated unknowns.
     * @param atTime Time step of interest.
     * @param dofId Dof mask, id set to NULL, all Dofs evaluated.
     * @return Error code (0=ok, 1=point not found in domain)
     */
    virtual int __evaluateAt(FloatArray &answer, FloatArray& coords,
                    ValueModeType mode, TimeStep *atTime, IntArray *dofId);
    /**
     * @param atTime Time step to take solution for.
     * @return Solution vector for requested time step.
     */
    virtual FloatArray *giveSolutionVector(TimeStep *atTime);

    /**
     * Project vectorToStore back to DOF's dictionary
     * @param vectorToStore Vector with the size of number of equations.
     * @param mode Mode of the unknown (increment, total value etc.)
     * @param atTime Time step unknowns belong to.
     */
    virtual void update(ValueModeType mode, TimeStep *atTime, FloatArray &vectorToStore) { };

    /**
     * Brings up a new solution vector for given time step.
     * @param atTime Time step for new solution vector.
     */
    virtual void advanceSolution(TimeStep *atTime);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode);

    /// @return Equation ID coupled to the field.
    EquationID giveEquationID() { return this->ut; }

protected:
    int resolveIndx(TimeStep *atTime, int shift);
    virtual FloatArray *giveSolutionVector(int);
};

} // end namespace oofem
#endif // primaryfield_h
