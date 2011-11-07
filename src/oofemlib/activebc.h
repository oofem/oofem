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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef activebc_h
#define activebc_h

#include "generalbc.h"
#include "alist.h"
#include "intarray.h"
#include "equationid.h"
#include "chartype.h"
#include "valuemodetype.h"

namespace oofem {
class SparseMtrx;
class UnknownNumberingScheme;

/**
 * Abstract base class for all active boundary conditions.
 * Design of active boundary conditions are subject to change.
 */
class ActiveBoundaryCondition : public GeneralBoundaryCondition
{
public:
    /**
     * Constructor. Creates boundary an active condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    ActiveBoundaryCondition(int n, Domain *d) : GeneralBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~ActiveBoundaryCondition() { }

    IRResultType initializeFrom(InputRecord *ir)
    {
        GeneralBoundaryCondition :: initializeFrom(ir);

        const char *__proc = "initializeFrom";
        IRResultType result;
        IntArray tempA, tempB, tempC;
        IR_GIVE_OPTIONAL_FIELD(ir, tempA, IFT_ActiveBoundaryCondition_elements, "elements");
        for (int i = 0; i < tempA.giveSize(); ++i) {
            this->addElement(tempA(i));
        }
        IR_GIVE_OPTIONAL_FIELD(ir, tempB, IFT_ActiveBoundaryCondition_elementSides, "elementsides");
        for (int i = 0; i < tempB.giveSize()/2; ++i) {
            this->addElementSide(tempB(i*2+1),tempB(i*2));
        }
        IR_GIVE_OPTIONAL_FIELD(ir, tempC, IFT_ActiveBoundaryCondition_dofManagers, "dofmans");
        for (int i = 0; i < tempB.giveSize()/2; ++i) {
            this->addElementSide(tempB(i*2+1),tempB(i*2));
        }

        return IRRT_OK;
    }

    /**@name Methods supporting classical input files*/
    //@{
    /**
     * Adds element for active boundary condition.
     * @param elem Element number.
     */
    virtual void addElement(int elem) { OOFEM_ERROR2("%s :: addElement - Not supported", giveClassName()); }

    /**
     * Adds element for active boundary condition.
     * @param elem Element number.
     * @param side Side number.
     */
    virtual void addElementSide(int elem, int side) { OOFEM_ERROR2("%s :: addElement - Not supported", giveClassName()); }

    /**
     * Adds dof manager for active boundary condition.
     * @param dman Dof manager number.
     */
    virtual void addDofman(int dman) { OOFEM_ERROR2("%s :: addElement - Not supported", giveClassName()); }
    //@}

    /**
     * Assembles B.C. contributions to specified matrix.
     * @param[in,out] answer Matrix to assemble to.
     * @param tStep Active time step.
     * @param tStep Active time step.
     * @param eid Equation ID.
     * @param type Type of matrix to assemble.
     * @param r_s Row numbering scheme.
     * @param c_s Column numbering scheme.
     * @param domain Domain to assemble from.
     */
    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain) {};

    /**
     * Assembles B.C. contributions to specified vector.
     * @param[in,out] answer Vector to assemble to.
     * @param tStep Active time step.
     * @param eid Equation ID.
     * @param type Type of matrix to assemble.
     * @param mode Mode of value.
     * @param s Numbering scheme.
     * @param domain Domain to assemble from.
     */
    virtual void assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain) {};

    /**
     * Gives a list of location arrays that will be assembled.
     * This should only be used to construct zero structure in sparse matrices.
     * @param rows List of location arrays for r_s.
     * @param cols List of location arrays for c_s.
     * @param eid Equation ID.
     * @param type Type of matrix to assemble.
     * @param r_s Row numbering scheme.
     * @param c_s Column numbering scheme.
     * @param domain Domain to assemble from.
     */
    virtual void giveLocationArrays(AList<IntArray> &rows, AList<IntArray> &cols, EquationID eid, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain) {};

    classType giveClassID() const { return ActiveBoundaryConditionClass; }
    const char *giveClassName() const { return "ActiveBoundaryCondition"; }
};
} // end namespace oofem
#endif // activebc_h

