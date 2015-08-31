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

#include "pfemnumberingschemes.h"

namespace oofem {
PressureNumberingScheme :: PressureNumberingScheme() :
    UnknownNumberingScheme()
    , neq(0)
    , pres_neq(0)
    , isInitialized(false)
{ }

PressureNumberingScheme :: ~PressureNumberingScheme()
{ }

void
PressureNumberingScheme :: init(Domain *domain, TimeStep *tStep)
{
    isInitialized = true;
    int inode;
    int nnode = domain->giveNumberOfDofManagers();
    DofManager *idofman;

    this->nodalPressureEquationNumbers.resize(nnode);
    for ( inode = 1; inode <= nnode; inode++ ) {
        idofman = domain->giveDofManager(inode);
        for ( Dof *jDof: *idofman ) {
            if ( jDof->giveDofID() == P_f ) {
                if ( jDof->hasBc(tStep) ) {
                    this->nodalPressureEquationNumbers.at(inode) = --pres_neq;
                } else {
                    this->nodalPressureEquationNumbers.at(inode) = ++neq;
                }
            }
        }
    }
}

void
PressureNumberingScheme :: reset()
{
    neq = 0;
    pres_neq = 0;
}

int
PressureNumberingScheme :: giveDofEquationNumber(Dof *dof) const
{
    int dofEqNum = 0;
    DofIDItem id = dof->giveDofID();
    if ( id == P_f ) {
        dofEqNum = this->nodalPressureEquationNumbers.at( dof->giveDofManNumber() );
	if ( dofEqNum < 0 )
	    dofEqNum = 0;
    }
    return dofEqNum;
}


int
PressureNumberingScheme :: giveTotalNumberOfEquations() const
{
    return neq;
}


int
PressureNumberingScheme :: giveRequiredNumberOfDomainEquation() const
{
    return neq;
}

int
PressureNumberingScheme :: giveTotalNumberOfPrescribedEquations() const
{
    return -1 * pres_neq;
}


VelocityNumberingScheme :: VelocityNumberingScheme(bool prescribed) :
    UnknownNumberingScheme()
    , numEqs(0)
    , prescribed(prescribed)
{ }

VelocityNumberingScheme :: ~VelocityNumberingScheme()
{ }

int
VelocityNumberingScheme :: giveDofEquationNumber(Dof *dof) const
{
    DofIDItem id = dof->giveDofID();
    if ( id == V_u || id == V_v || id == V_w ) {
        return prescribed ? dof->__givePrescribedEquationNumber() : dof->__giveEquationNumber();
    }

    return 0;
}

Dof *
VelocityNumberingScheme :: giveDofToEquationNumber(Domain *d, int equationNumber)
{
    bool found = false;
    Dof *foundDof = NULL;
    for ( int i = 1; i <= d->giveNumberOfDofManagers(); i++ ) {
        DofManager *dman = d->giveDofManager(i);
        for ( Dof *dof : *dman ) {//= d->giveDofManager(i)->giveDof(j);
            DofIDItem id = dof->giveDofID();
            if ( id == V_u || id == V_v || id == V_w ) {
                if ( dof->__giveEquationNumber() == equationNumber ) {
                    foundDof = dof;
                    found = true;
                    break;
                }
            }
        }

        if ( found ) {
            break;
        }
    }

    return foundDof;
}

AuxVelocityNumberingScheme :: AuxVelocityNumberingScheme() :
    UnknownNumberingScheme()
    , neq(0)
{ }

AuxVelocityNumberingScheme :: ~AuxVelocityNumberingScheme()
{ }

void
AuxVelocityNumberingScheme :: init(Domain *domain)
{
    neq = domain->giveNumberOfSpatialDimensions() * domain->giveNumberOfDofManagers();
}

int
AuxVelocityNumberingScheme :: giveDofEquationNumber(Dof *dof) const
{
    DofIDItem type = dof->giveDofID();
    int c;
    int n = dof->giveDofManager()->giveGlobalNumber();

    if ( type == V_u ) {
        c = 1;
    } else if ( type == V_v ) {
        c = 2;
    } else if ( type == V_w ) {
        c = 3;
    } else {
        return 0;
    }

    return dof->giveDofManager()->giveDomain()->giveNumberOfSpatialDimensions() * ( n - 1 ) + c;
}

int
AuxVelocityNumberingScheme :: giveRequiredNumberOfDomainEquation() const
{
    return neq;
}
} // end namespace oofem
