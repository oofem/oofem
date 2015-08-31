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

#ifndef pfemnumberingschemes_h
#define pfemnumberingschemes_h

#include "unknownnumberingscheme.h"
#include "dof.h"
#include "domain.h"
#include "dofmanager.h"

namespace oofem {
/**
 * Numbering scheme that takes into account only pressure DOFs in PFEM problems.
 * Pressure is scalar unknown in each node.
 *
 * @author David Krybus
 */
class PressureNumberingScheme : public UnknownNumberingScheme
{
protected:
    Domain *domain;
    /// Container storing particular equation numbers for each node
    IntArray nodalPressureEquationNumbers;
    /// Last given number of equation
    int neq;
    /// Last given number of prescribed equation
    int pres_neq;
    /// Flag controlling wether the numbering has been initialized or not
    bool isInitialized;

public:
    /// Constructor
    PressureNumberingScheme();
    /// Destructor
    virtual ~PressureNumberingScheme();

    /**
     * Initializes the receiver
     */
    virtual void init(Domain *domain, TimeStep *tStep);
    virtual bool isDefault() const { return false; }
    virtual int giveDofEquationNumber(Dof *dof) const;
    virtual int giveRequiredNumberOfDomainEquation() const;

    /// Returns total number of equations
    virtual int giveTotalNumberOfEquations() const;

    /// Returns total number of prescribed equations
    virtual int giveTotalNumberOfPrescribedEquations() const;

    /// Resets the numbering in order to start numbering again from 1
    virtual void reset();
};

/**
 * Velocity numbering scheme for PFEM purposes
 *
 * @author David Krybus
 */
class VelocityNumberingScheme : public UnknownNumberingScheme
{
protected:
    /// Last given equation number
    int numEqs;
    /// prescribed equations or not
    bool prescribed;

public:
    /// Constructor
    VelocityNumberingScheme(bool prescribed);
    /// Destructor
    virtual ~VelocityNumberingScheme();

    /// Resets the numbering in order to start numbering again from 1
    void reset() { numEqs = 0; }
    virtual bool isDefault() const { return !prescribed; }
    virtual int giveDofEquationNumber(Dof *dof) const;
    virtual int giveRequiredNumberOfDomainEquation() const { return numEqs; }

    /// Asks new equation number
    int askNewEquationNumber() { return ++numEqs; }
    /// Returns Dof with passed equation number
    Dof *giveDofToEquationNumber(Domain *d, int equationNumber);
};

/**
 * Numbering scheme for auxiliary velocity in PFEM problems.
 * Nodes are not subjected to any boundary conditions considering auxiliary velocities,
 * so every single node has full set of equations.
 *
 * @author David Krybus
 */
class AuxVelocityNumberingScheme : public UnknownNumberingScheme
{
protected:
    Domain *domain;
    /// Nodal equation numbers are stored in an IntArray
    IntArray nodalAuxVelocityEquationNumbers;
    /// Number of equations
    int neq;

public:
    /// Constructor
    AuxVelocityNumberingScheme();
    /// Destructor
    virtual ~AuxVelocityNumberingScheme();

    /// Resets the numbering in order to start numbering again from 1
    virtual void reset() { neq = 0; }

    /// Initializes the numbering scheme
    virtual void init(Domain *domain);

    virtual int giveDofEquationNumber(Dof *dof) const;
    virtual int giveRequiredNumberOfDomainEquation() const;
};
} // end namespace oofem
#endif // pfemnumberingschemes_h
