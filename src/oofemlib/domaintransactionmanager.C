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

#include "domaintransactionmanager.h"
#include "error.h"
#include "dofmanager.h"
#include "element.h"
#include "domain.h"

namespace oofem {
DomainTransactionManager :: DomainTransactionManager(Domain *d)
{
    this->domain = d;
}

DomainTransactionManager :: ~DomainTransactionManager()
{
    if ( !( dofmanTransactions.empty() && elementTransactions.empty() ) ) {
        OOFEM_WARNING("uncommited transactions exist");
    }
}

void
DomainTransactionManager :: initialize()
{ }

int
DomainTransactionManager :: addDofManTransaction(DomainTransactionType dtt, int label, DofManager *obj)
{
    if ( dtt == DTT_Remove ) {
        obj = NULL;
    }

    if ( dofmanTransactions.find(label) != dofmanTransactions.end() ) {
        // local enry exist
        // delete previous record
        if ( dofmanTransactions [ label ] ) {
            delete dofmanTransactions [ label ];
        }
    }

    // set new record
    dofmanTransactions [ label ] = obj;

    return 1;
}

int
DomainTransactionManager :: addElementTransaction(DomainTransactionType dtt,  int label, Element *obj)
{
    if ( dtt == DTT_Remove ) {
        obj = NULL;
    }

    if ( elementTransactions.find(label) != elementTransactions.end() ) {
        // local enry exist
        // delete previous record
        if ( elementTransactions [ label ] ) {
            delete elementTransactions [ label ];
        }
    }

    // set new record
    elementTransactions [ label ] = obj;

    return 1;
}


DofManager *DomainTransactionManager :: giveDofManager(int label)
{
    if ( dofmanTransactions.find(label) != dofmanTransactions.end() ) {
        // if modified record exist return it
        return dofmanTransactions [ label ];
    } else {
        // no modification recorded -> return NULL
        return NULL;
    }
}


Element *DomainTransactionManager :: giveElement(int label)
{
    if ( elementTransactions.find(label) != elementTransactions.end() ) {
        // if modified record exist return it
        return elementTransactions [ label ];
    } else {
        // no modification recorded -> return NULL
        return NULL;
    }
}


int
DomainTransactionManager :: commitTransactions()
{
    return domain->commitTransactions(this);
}
} // end namespace oofem
