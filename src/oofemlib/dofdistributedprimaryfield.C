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

#include "dofdistributedprimaryfield.h"
#include "dofmanager.h"
#include "dof.h"
#include "domain.h"
#include "timestep.h"
#include "floatarray.h"
#include "engngm.h"

namespace oofem {
DofDistributedPrimaryField :: DofDistributedPrimaryField(EngngModel *a, int idomain,
                                                         FieldType ft, int nHist) :
    PrimaryField(a, idomain, ft, nHist)
{ }

DofDistributedPrimaryField :: ~DofDistributedPrimaryField()
{ }

double
DofDistributedPrimaryField :: giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *tStep)
{
    return dof->giveUnknown(mode, tStep);
}

FloatArray *
DofDistributedPrimaryField :: giveSolutionVector(TimeStep *tStep)
{
    return PrimaryField :: giveSolutionVector(tStep);
}

void
DofDistributedPrimaryField :: initialize(ValueModeType mode, TimeStep *tStep, FloatArray &answer, const UnknownNumberingScheme &s)
{
    Domain *domain = emodel->giveDomain(domainIndx);
    int neq =  emodel->giveNumberOfDomainEquations(domainIndx, s);
    int nnodes = domain->giveNumberOfDofManagers();

    answer.resize(neq);
    answer.zero();

    for ( int j = 1; j <= nnodes; j++ ) {
        DofManager *inode = domain->giveDofManager(j);
        for ( Dof *iDof: *inode ) {
            int eqNum = iDof->__giveEquationNumber();
            double val;
            if ( eqNum ) {
                iDof->giveUnknownsDictionaryValue(tStep, mode, val);
                answer.at(eqNum) = val;
                //answer.at(eqNum) = iDof->giveUnknown( mode, tStep);
            }
        }
    }
}

// project solutionVector to DoF unknowns dictionary
void
DofDistributedPrimaryField :: update(ValueModeType mode, TimeStep *tStep, FloatArray &vectorToStore)
{
    Domain *domain = emodel->giveDomain(domainIndx);
    int nnodes = domain->giveNumberOfDofManagers();

    for ( int j = 1; j <= nnodes; j++ ) {
        DofManager *inode = domain->giveDofManager(j);
        for ( Dof *iDof: *inode ) {
            int eqNum = iDof->__giveEquationNumber();
            double val;
            if ( mode == VM_Total ) {
                if ( iDof->hasBc(tStep) ) { // boundary condition
                    val = iDof->giveBcValue(VM_Total, tStep);
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

            iDof->updateUnknownsDictionary(tStep, mode, val);
        }
    }
}



void
DofDistributedPrimaryField :: advanceSolution(TimeStep *tStep)
{
    PrimaryField :: advanceSolution(tStep);
}


contextIOResultType
DofDistributedPrimaryField :: saveContext(DataStream &stream, ContextMode mode)
{
    // all the job is done by dofs alone
    return CIO_OK;
}

contextIOResultType
DofDistributedPrimaryField :: restoreContext(DataStream &stream, ContextMode mode)
{
    return CIO_OK;
}
} // end namespace oofem
