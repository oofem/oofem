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

#ifndef unknownnumberingscheme_h
#define unknownnumberingscheme_h

#include "dof.h"
#include "intarray.h"

namespace oofem {
/**
 * Abstract base class allowing to control the way, how equations are assigned to individual DOFs.
 * The instances are typically used in EngngModel to assemble characteristic contributions and they
 * allow to control the numbering of unknowns.
 */
class OOFEM_EXPORT UnknownNumberingScheme
{
public:
    UnknownNumberingScheme(void) { }
    virtual ~UnknownNumberingScheme() { }

    /**
     * Initializes the receiver, if necessary.
     */
    virtual void init() { }
    /**
     * Returns true, if receiver is the default engngModel equation numbering scheme;
     * This is useful for some components (typically elements), that cache their code numbers
     * for default numbering to avoid repeated evaluation.
     */
    virtual bool isDefault() const { return false; }
    /**
     * Returns the equation number for corresponding DOF. The numbering should return nonzero value if
     * the equation is assigned to the given DOF, zero otherwise.
     */
    virtual int giveDofEquationNumber(Dof *dof) const = 0;

    /**
     * Returns required number of domain equation. Number is always less or equal to the sum of all DOFs gathered from all nodes.
     */
    virtual int giveRequiredNumberOfDomainEquation() const { return 0; }
};

/**
 * The representation of EngngModel default unknown numbering. The equation numbers are assigned
 * by the engineering model itself to individual DOFs. Therefore, this call is a simple shell around
 * DofEquationNumbering interface, forwarding all the requests to individual DOFs.
 */
class OOFEM_EXPORT EModelDefaultEquationNumbering : public UnknownNumberingScheme
{
public:
    EModelDefaultEquationNumbering(void) : UnknownNumberingScheme() { }

    virtual bool isDefault() const { return true; }
    virtual int giveDofEquationNumber(Dof *dof) const {
        return dof->__giveEquationNumber();
    }
};

/**
 * The representation of EngngModel default prescribed unknown numbering.
 * The equation numbers are assigned by the engineering model itself to individual DOFs.
 * Therefore, this call is a simple shell around
 * DofEquationNumbering interface, forwarding all the requests to individual DOFs.
 */
class OOFEM_EXPORT EModelDefaultPrescribedEquationNumbering : public UnknownNumberingScheme
{
public:
    EModelDefaultPrescribedEquationNumbering(void) : UnknownNumberingScheme() { }

    virtual int giveDofEquationNumber(Dof *dof) const {
        return dof->__givePrescribedEquationNumber();
    }
};


/**
 * Specialized numbering scheme for assembling only specified DofIDs.
 * This can be useful for computing the reaction forces for a subset of dofids.
 */
class DofIDEquationNumbering : public UnknownNumberingScheme
{
protected:
    IntArray dofids;
    bool prescribed;

public:
    DofIDEquationNumbering(bool prescribed, IntArray dofids) : 
        UnknownNumberingScheme(), dofids(std :: move(dofids)), prescribed(prescribed) { }

    virtual bool isDefault() const { return !prescribed; }
    virtual int giveDofEquationNumber(Dof *dof) const {
        DofIDItem id = dof->giveDofID();
        if ( dofids.contains(id) ) {
            return prescribed ? dof->__givePrescribedEquationNumber() : dof->__giveEquationNumber();
        }

        return 0;
    }
};
} // end namespace oofem
#endif // unknownnumberingscheme_h
