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

#include "dofdistributedprimaryfield.h"
#include "dofmanager.h"
#include "dof.h"
#include "timestep.h"
#include "flotarry.h"

namespace oofem {
DofDistributedPrimaryField :: DofDistributedPrimaryField(EngngModel *a, int idomain,
                                                         FieldType ft, EquationID ut, int nHist) :
    PrimaryField(a, idomain, ft, ut, nHist)
{ }

DofDistributedPrimaryField :: ~DofDistributedPrimaryField()
{ }

double
DofDistributedPrimaryField :: giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *atTime)
{
    return dof->giveUnknown(mode, atTime);
}

FloatArray *
DofDistributedPrimaryField :: giveSolutionVector(TimeStep *atTime)
{
    return PrimaryField :: giveSolutionVector(atTime);
}

void
DofDistributedPrimaryField :: initialize(ValueModeType mode, TimeStep *atTime, FloatArray &answer, const UnknownNumberingScheme &s)
{
    Domain *domain = emodel->giveDomain(domainIndx);
    int neq =  emodel->giveNumberOfDomainEquations(domainIndx, s);
    int nnodes = domain->giveNumberOfDofManagers();

    answer.resize(neq);
    answer.zero();

    for ( int j = 1; j <= nnodes; j++ ) {
        DofManager *inode = domain->giveDofManager(j);
        int ndofs = inode->giveNumberOfDofs();
        for ( int i = 1; i <= ndofs; i++ ) {
            Dof *iDof = inode->giveDof(i);
            int eqNum = iDof->__giveEquationNumber();
            double val;
            if ( eqNum ) {
                iDof->giveUnknownsDictionaryValue(atTime, mode, val);
                answer.at(eqNum) = val;
                //answer.at(eqNum) = iDof->giveUnknown( mode, atTime);
            }
        }
    }
}

// project solutionVector to DoF unknowns dictionary
void
DofDistributedPrimaryField :: update(ValueModeType mode, TimeStep *atTime, FloatArray &vectorToStore)
{
    Domain *domain = emodel->giveDomain(domainIndx);
    int nnodes = domain->giveNumberOfDofManagers();

    for ( int j = 1; j <= nnodes; j++ ) {
        DofManager *inode = domain->giveDofManager(j);
        int ndofs = inode->giveNumberOfDofs();
        for ( int i = 1; i <= ndofs; i++ ) {
            Dof *iDof = inode->giveDof(i);
            int eqNum = iDof->__giveEquationNumber();
            double val;
            if ( mode == VM_Total ) {
                if ( iDof->hasBc(atTime) ) { // boundary condition
                    val = iDof->giveBcValue(VM_Total, atTime);
                } else {
                    //vect = this->UnknownsField->giveSolutionVector(tStep);
                    val = vectorToStore.at(eqNum);
                }
            } else { //all other modes, e.g. VM_RhsTotal
                if ( !eqNum ) {
                    val = 0.; //assume that 0's are present in the beginning of node initiation
                } else {
                    val = vectorToStore.at(eqNum);
                }
            }

            iDof->updateUnknownsDictionary(atTime, mode, val);
        }
    }
}



void
DofDistributedPrimaryField :: advanceSolution(TimeStep *atTime)
{
    PrimaryField :: advanceSolution(atTime);
}


contextIOResultType
DofDistributedPrimaryField :: saveContext(DataStream *stream, ContextMode mode)
{
    // all the job is done by dofs alone
    return CIO_OK;
}

contextIOResultType
DofDistributedPrimaryField :: restoreContext(DataStream *stream, ContextMode mode)
{
    return CIO_OK;
}
} // end namespace oofem
