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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
#ifndef domaintransactionmanager_h
#define domaintransactionmanager_h

#ifdef __PARALLEL_MODE
#ifndef __MAKEDEPEND
#include <map>
#endif

namespace oofem {

class FEMComponent;
class Domain;
class DofManager;
class Element;
/**
 * Class representing domain transaction mamager. The purpose of this class is to
 * make the domain modification in terms of adding and deleting components) possible.
 *
 * The changes are recorded in transaction manager and untill the are commited,
 * no change is reflected in domain itself.
 *
 * When adding componentsm new component has to be created and send as transaction parameter.
 * The allocated component will be further maintained (and dealocated) in transaction manager or domain itself.
 * When adding component a reference to already existing as well as newly added (but not yet commited)
 * component must be possible. Therefore, the concept of temporary numbering is introduced.
 * Existing domain components are assigned with temporary numbers equal to to their local numbers
 * (in parallel version to their global numbers)  and newly addded components can be assigned with
 * arbitrary temp number (but which should be different from those already used). Then references
 * between components for newly added components is specified using temporary numbering.
 * After transactions are commited, the local numbering can change.
 *
 * The user is responsible for data consistency, so for example if a node is removed,
 * but an element exist that needs such node, then (hopefully)
 * an error will be thrown since there is no mechanism for enforcing data consistency.
 *
 * Transaction manager provides an access to recorded transactions; so it is possible
 * to modify an existing transaction instead of posting a new one.
 * Hovewer, transaction manager does not provides access to original umnodified domain components,
 * since this can lead to their modification and since in the current implementation the transaction
 * manager keeps only modified records, the original records are not returned. They must be requested from original domain,
 * so it is clear that their modification is not part of transaction but it directly changes domain data.
 *
 */
class DomainTransactionManager
{
public:
    enum DomainComponentType { DCT_DofManager, DCT_Element };
    enum DomainTransactionType { DTT_ADD, DTT_Remove };
    /*
     * struct DTM_Transaction {
     * /// transaction type
     * DomainTransactionType _ttype;
     * /// component type
     * DomainComponentType _ctype;
     * /// temporary number
     * int _num;
     * /// component adress
     * FemComponent* _obj;
     * };
     */
protected:
    Domain *domain;
    std :: map< int, FEMComponent * >dofmanTransactions, elementTransactions;
    //list<DTM_Transaction> transactions;

public:
    DomainTransactionManager(Domain *d);
    ~DomainTransactionManager();

    void initialize();
    int addTransaction(DomainTransactionType, DomainComponentType, int, FEMComponent *);
    DofManager *giveDofManager(int label);
    Element *giveElement(int labal);
    int  commitTransactions();


    friend class Domain;
};

} // end namespace oofem
#endif
#endif // domaintransactionmanager_h




