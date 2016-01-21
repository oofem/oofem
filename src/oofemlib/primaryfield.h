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

#ifndef primaryfield_h
#define primaryfield_h

#include "field.h"
#include "interface.h"
#include "floatarray.h"
#include "valuemodetype.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "timestep.h"

#include <vector>

namespace oofem {
class PrimaryField;
class Dof;
class BoundaryCondition;
class InitialCondition;
class UnknownNumberingScheme;

/**
 * Element interface class. Declares the functionality required to support PrimaryField element interpolation.
 */
class OOFEM_EXPORT EIPrimaryFieldInterface : public Interface
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
     * @todo This should use local coordinates instead of having all elements search for it manually.
     * @todo Shouldn't this just be part of Element? It's very much the core of functionality for elements.
     * 
     * @param answer Field evaluated at coordinate.
     * @param pf Field to use for evaluation.
     * @param coords Coordinate.
     * @param dofId IDs of DOFs to evaluate.
     * @param mode Mode of field.
     * @param tStep Time step to evaluate at.
     * @return Zero if ok, nonzero when error encountered.
     */
    virtual int EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                      FloatArray &coords, IntArray &dofId, ValueModeType mode, TimeStep *tStep) = 0;
    //@}
};


/**
 * Abstract class representing field of primary variables (those, which are unknown and are typically associated to nodes).
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
 * 
 * @note This primary field will always assume default numbering schemes.
 */
class OOFEM_EXPORT PrimaryField : public Field
{
protected:
    int actualStepNumber;
    int actualStepIndx;
    int nHistVectors;
    std :: vector< FloatArray >solutionVectors;
    std :: vector< FloatArray >prescribedVectors;
    std :: vector< TimeStep >solStepList;
    EngngModel *emodel;
    int domainIndx;

public:
    /**
     * Constructor. Creates a field of given type associated to given domain.
     * Not using pointer to domain, because this will prevent the use of PrimaryField as an
     * EngngModel attribute. This is because the domain does not exists when
     * PrimaryField is created (this is when EngngModel is created).
     * @param a Engineering model which field belongs to.
     * @param idomain Index of domain for field.
     * @param ft Type of stored field.
     * @param nHist Number of old time steps to store.
     */
    PrimaryField(EngngModel * a, int idomain, FieldType ft, int nHist);
    virtual ~PrimaryField();

    /**
     * Copy unknowns from previous solution or DOF's dictionary to the solution vector
     * @param mode what the unknown describes (increment, total value etc.).
     * @param tStep Time of interest.
     * @param answer Resulting vector.
     */
    virtual void initialize(ValueModeType mode, TimeStep *tStep, FloatArray &answer, const UnknownNumberingScheme &s);

    // These functions are hardcoded to assume the default numbering scheme (as is the rest of primaryfield)
    void storeDofManager(TimeStep *tStep, DofManager &dman);
    void storeInDofDictionaries(TimeStep *tStep);
    void readDofManager(TimeStep *tStep, DofManager &dman);
    void readFromDofDictionaries(TimeStep *tStep);

    /**
     * Applies the default initial values values for all DOFs (0) in given domain.
     * @param domain Domain number
     */
    virtual void applyDefaultInitialCondition();
    /**
     * Applies initial condition to all DOFs.
     * @param ic Initial condition for DOFs
     */
    void applyInitialCondition(InitialCondition &ic);
    /**
     * Applies all boundary conditions to all prescribed DOFs.
     * @param tStep Current time step.
     */
    virtual void applyBoundaryCondition(TimeStep *tStep);
    /**
     * Applies the boundary condition to all prescribed DOFs in given domain.
     * @param bc Boundary condition.
     * @param domain tStep Time step for when bc applies.
     */
    void applyBoundaryCondition(BoundaryCondition &bc, TimeStep *tStep);

    /**
     * @param dof Pointer to DOF.
     * @param mode What the unknown describes (increment, total value etc.).
     * @param tStep Time step of interest.
     * @return Value of interest at given DOF.
     */
    virtual double giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *tStep);
    /**
     * Evaluates the field at given point
     * @param answer Evaluated field at point.
     * @param coords Coordinates of the point of interest.
     * @param mode Mode of evaluated unknowns.
     * @param tStep Time step of interest.
     * @return Error code (0-ok, 1-point not found in domain).
     */
    virtual int evaluateAt(FloatArray &answer, FloatArray &coords,
                           ValueModeType mode, TimeStep *tStep);
    /**
     * Evaluates the field at given DofManager
     * @param answer Evaluated field at dman.
     * @param dman DOF manager to evaluate at.
     * @param mode Mode of evaluated unknowns.
     * @param tStep Time step of interest.
     * @return Error code (0-ok, 1-point not found in domain).
     */
    virtual int evaluateAt(FloatArray &answer, DofManager *dman,
                           ValueModeType mode, TimeStep *tStep);
    /**
     * Evaluates the field at given DOF manager, allows to select specific
     * dofs using mask
     * @param answer Evaluated field at dman.
     * @param dman DOF manager of interest.
     * @param mode Mode of evaluated unknowns.
     * @param tStep Time step of interest.
     * @param dofId Dof mask, id set to NULL, all Dofs evaluated.
     * @return Error code (0=ok, 1=point not found in domain).
     */
    virtual int __evaluateAt(FloatArray &answer, DofManager *dman,
                             ValueModeType mode, TimeStep *tStep, IntArray *dofId);
    /**
     * Evaluates the field at given point, allows to select specific
     * dofs using mask.
     * @param answer Evaluated field at coords.
     * @param coords Coordinates of the point of interest.
     * @param mode Mode of evaluated unknowns.
     * @param tStep Time step of interest.
     * @param dofId Dof mask, id set to NULL, all Dofs evaluated.
     * @return Error code (0=ok, 1=point not found in domain)
     */
    virtual int __evaluateAt(FloatArray &answer, FloatArray &coords,
                             ValueModeType mode, TimeStep *tStep, IntArray *dofId);
    /**
     * @param tStep Time step to take solution for.
     * @return Solution vector for requested time step.
     */
    virtual FloatArray *giveSolutionVector(TimeStep *tStep);

    /**
     * Project vectorToStore back to DOF's dictionary
     * @param vectorToStore Vector with the size of number of equations.
     * @param mode Mode of the unknown (increment, total value etc.)
     * @param tStep Time step unknowns belong to.
     */
    virtual void update(ValueModeType mode, TimeStep *tStep, const FloatArray &vectorToStore, const UnknownNumberingScheme &s);

    /**
     * Brings up a new solution vector for given time step.
     * @param tStep Time step for new solution vector.
     */
    virtual void advanceSolution(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode);

    virtual const char *giveClassName() const { return "PrimaryField"; }


    int giveActualStepNumber() { return actualStepNumber; }
protected:
    int resolveIndx(TimeStep *tStep, int shift);
    FloatArray *giveSolutionVector(int);
    FloatArray *givePrescribedVector(int);
};
} // end namespace oofem
#endif // primaryfield_h
