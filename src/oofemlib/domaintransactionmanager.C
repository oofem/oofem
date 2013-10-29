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

#ifdef __PARALLEL_MODE
#include "domaintransactionmanager.h"
#include "error.h"
#include "femcmpnn.h"
#include "domain.h"

namespace oofem {
DomainTransactionManager :: DomainTransactionManager(Domain *d)
{
    this->domain = d;
}

DomainTransactionManager :: ~DomainTransactionManager()
{
    if ( !( dofmanTransactions.empty() && elementTransactions.empty() ) ) {
        OOFEM_WARNING("DomainTransactionManager::~DomainTransactionManager: uncommited transactions exist");
    }
}

void
DomainTransactionManager :: initialize()
{ }

int
DomainTransactionManager :: addTransaction(DomainTransactionType dtt, DomainComponentType dct, int label, FEMComponent *obj)
{
    std :: map< int, FEMComponent * > *rec;
    if ( dct == DCT_Element ) {
        rec = & elementTransactions;
    } else {
        rec = & dofmanTransactions;
    }

    if ( dtt == DTT_Remove ) {
        obj = NULL;
    }

    if ( rec->find(label) != rec->end() ) {
        // local enry exist
        // delete previous record
        if ( ( * rec ) [ label ] ) {
            delete(* rec) [ label ];
        }
    }

    // set new record
    ( * rec ) [ label ] = obj;

    return 1;
}


DofManager *DomainTransactionManager :: giveDofManager(int label)
{
    DofManager *answer;
    if ( dofmanTransactions.find(label) != dofmanTransactions.end() ) {
        // if modified record exist return it
        answer = ( DofManager * ) dofmanTransactions [ label ];
    } else {
        // no modification recorded -> return NULL
        answer = NULL;
    }

    return answer;
}


Element *DomainTransactionManager :: giveElement(int label)
{
    Element *answer;
    if ( elementTransactions.find(label) != elementTransactions.end() ) {
        // if modified record exist return it
        answer = ( Element * ) elementTransactions [ label ];
    } else {
        // no modification recorded -> return NULL
        answer = NULL;
    }

    return answer;
}


int
DomainTransactionManager :: commitTransactions()
{
    return domain->commitTransactions(this);
}
} // end namespace oofem
#endif
