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

#ifndef activebc_h
#define activebc_h

#include "generalboundarycondition.h"
#include "intarray.h"
#include "inputrecord.h"
#include "chartype.h"
#include "valuemodetype.h"
#include "error.h"

#include <vector>

///@name Input fields for active boundary condition
//@{
#define _IFT_ActiveBoundaryCondition_elementSides "elementsides"
//@}

namespace oofem {
class SparseMtrx;
class UnknownNumberingScheme;
class ActiveDof;
class Dof;
class PrimaryField;

/**
 * Abstract base class for all active boundary conditions.
 * Design of active boundary conditions are subject to change.
 */
class OOFEM_EXPORT ActiveBoundaryCondition : public GeneralBoundaryCondition
{
public:
    /**
     * Constructor. Creates boundary an active condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    ActiveBoundaryCondition(int n, Domain * d) : GeneralBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~ActiveBoundaryCondition() { }

    IRResultType initializeFrom(InputRecord *ir)
    {
        GeneralBoundaryCondition :: initializeFrom(ir);

        IRResultType result;
        IntArray tempA, tempB, tempC;
        IR_GIVE_OPTIONAL_FIELD(ir, tempB, _IFT_ActiveBoundaryCondition_elementSides);
        for ( int i = 0; i < tempB.giveSize() / 2; ++i ) {
            this->addElementSide( tempB(i * 2), tempB(i * 2 + 1) );
        }

        return IRRT_OK;
    }

    ///  @name Methods supporting classical input files.
    //@{
    /**
     * Adds element for active boundary condition.
     * @param elem Element number.
     * @param side Side number.
     */
    virtual void addElementSide(int elem, int side) { OOFEM_ERROR("Not supported"); }
    //@}

    /**
     * Assembles B.C. contributions to specified matrix.
     * @param[in,out] answer Matrix to assemble to.
     * @param tStep Active time step.
     * @param tStep Active time step.
     * @param type Type of matrix to assemble.
     * @param r_s Row numbering scheme.
     * @param c_s Column numbering scheme.
     */
    virtual void assemble(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) { }

    /**
     * Assembles B.C. contributions to specified vector.
     * @param[in,out] answer Vector to assemble to.
     * @param tStep Active time step.
     * @param type Type of matrix to assemble.
     * @param mode Mode of value.
     * @param s Numbering scheme.
     * @param eNorms Norms for each dofid (optional).
     * @return Equivalent of the sum of the element norm (squared) of assembled vector.
     */
    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorms = NULL) { }

    /**
     * Gives a list of location arrays that will be assembled.
     * This should only be used to construct zero structure in sparse matrices.
     * The rows and columns location arrays returned in tuples (stored in vector),
     * allowing to efficiently assemble and allocate off-diagonal blocks.
     * The nonzero entries are assembled and allocated for entries at (rows[i], cols[i]) positions.
     * @param rows List of location arrays for r_s.
     * @param cols List of location arrays for c_s.
     * @param type Type of matrix to assemble.
     * @param r_s Row numbering scheme.
     * @param c_s Column numbering scheme.
     */
    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) { }


    /**
     * @name Functions related to boundary conditions which have to deal with special slave DOFs.
     * Boundary conditions that do not require slave dofs do not need to overload these functions.
     * @see ActiveDof The active DOF passes these functions here.
     */
    //@{
    /**
     * Checks to see if active boundary condition requires special DOFs.
     * @return True if ActiveDof should be created.
     */
    virtual bool requiresActiveDofs() { return false; }
    /**
     * Checks to see if the dof is a primary DOF.
     * @return True if ActiveDof is a primary DOF.
     */
    virtual bool isPrimaryDof(ActiveDof *dof) { return true; }
    /**
     * Returns the prescribed value of a dof (if any).
     */
    virtual double giveBcValue(Dof *dof, ValueModeType mode, TimeStep *tStep) { return 0.0; }
    /**
     * Returns the prescribed value of a dof (if any).
     * @return True if dof is prescribed.
     */
    virtual bool hasBc(Dof *dof, TimeStep *tStep) { return false; }
    /**
     * Allows for active boundary conditions to handle their own special DOF.
     * @param dof Active dof belonging to receiver.
     * @return Number of primary master DOFs.
     */
    virtual int giveNumberOfMasterDofs(ActiveDof *dof)
    {
        OOFEM_ERROR("Not supported by bc.");
        return 0;
    }
    /**
     * Give the pointer to master dof belonging to active DOF.
     * @param dof Active dof belonging to receiver.
     * @param mdof Local master dof number.
     * @return Master dof.
     */
    virtual Dof *giveMasterDof(ActiveDof *dof, int mdof)
    {
        OOFEM_ERROR("Not supported by bc.");
        return NULL;
    }
    virtual void computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs)
    {
        OOFEM_ERROR("Not supported by bc.");
    }
    /**
     * Computes the value of the dof.
     * @param field Field to take value from.
     * @param mode Mode of unknown value.
     * @param tStep Time step.
     * @param dof Active dof for which to obtain the value.
     * @return Value of dof.
     */
    virtual double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
    {
        OOFEM_ERROR("Not supported by bc.");
        return 0.0;
    }
    /**
     * Computes the value of the dof.
     * @param mode Mode of unknown value.
     * @param tStep Time step.
     * @param dof Active dof for which to obtain the value.
     * @return Value of dof.
     */
    virtual double giveUnknown(ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
    {
        OOFEM_ERROR("Not supported by bc.");
        return 0.0;
    }
    //@}
};
} // end namespace oofem
#endif // activebc_h
