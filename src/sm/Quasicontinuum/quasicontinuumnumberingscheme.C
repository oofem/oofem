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

#include "quasicontinuumnumberingscheme.h"
#include <algorithm>

namespace oofem {
QuasicontinuumNumberingscheme :: QuasicontinuumNumberingscheme() :
    UnknownNumberingScheme()
    , neq(0)
    , pres_neq(0)
    , isInitialized(false)
{ }

QuasicontinuumNumberingscheme :: ~QuasicontinuumNumberingscheme()
{
    delete equationMap;
}

void
QuasicontinuumNumberingscheme :: init(Domain *domain, std :: vector< bool >activatedNodeList, TimeStep *tStep)
{
    isInitialized = true;
    /*
     * // compute list of selected nodes numbers
     * int noActovedNodes = std::count_if(activatedNodeList.begin(), activatedNodeList.end(), [](bool i){return i;});
     * this->selectedNodes.resize(noActovedNodes);
     * selectedNodes.zero();
     * int k = 0;
     * for (int i=1; i<=(int)activatedNodeList.size(); i++) {
     *  if (activatedNodeList[i-1]) {
     *      k++;
     *      selectedNodes.at(k) = i;
     *  }
     * }
     */
    this->selectedNodes.resize( activatedNodeList.size() );
    selectedNodes.zero();
    for ( int i = 1; i <= ( int ) activatedNodeList.size(); i++ ) {
        if ( activatedNodeList [ i - 1 ] ) {
            selectedNodes.at(i) = 1;
        }
    }


    //int nnode = domain->giveNumberOfDofManagers();
    //int nnode = 0;
    DofManager *idofman;

    //   this->dofEquationNumbers.resize(nnode);

    equationMap = new std :: map< int, std :: map< int, int > >();


    for ( int inode = 1; inode <= selectedNodes.giveSize(); inode++ ) {
        if ( selectedNodes.at(inode) == 1 ) {
            std :: map< int, int > *dof2EquationMap = new std :: map< int, int >();
            idofman = domain->giveDofManager(inode);
            IntArray dofIDArray;
            idofman->giveCompleteMasterDofIDArray(dofIDArray);
            if ( !idofman->hasAnySlaveDofs() ) {
                for ( int k = 1; k <= dofIDArray.giveSize(); k++ ) {
                    if ( idofman->giveDofWithID( dofIDArray.at(k) )->hasBc(tStep) ) {
                        ( * dof2EquationMap ) [ dofIDArray.at(k) ] = --pres_neq;
                    } else {
                        ( * dof2EquationMap ) [ dofIDArray.at(k) ] = ++neq;
                    }
                }
            } else {
                /*
                 * IntArray masterDofIDArray, masterDofMans;
                 * //idofman->giveMasterDofIDArray({D_u, D_v},masterDofIDArray);
                 * idofman->giveMasterDofMans(masterDofMans);
                 * DofManager *mdofman;
                 * for (int m = 1; m <= masterDofMans.giveSize(); m++ ) {
                 *  mdofman = domain->giveDofManager(masterDofMans.at(m));
                 *  mdofman->giveCompleteMasterDofIDArray(dofIDArray);
                 *  for (int k = 1; k <= dofIDArray.giveSize(); k++) {
                 *      if(mdofman->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)) {
                 *          ( *dof2EquationMap ) [ dofIDArray.at(k) ] = --pres_neq;
                 *      } else {
                 *          ( *dof2EquationMap ) [ dofIDArray.at(k) ] = ++neq;
                 *      }
                 *  }
                 * }
                 */
            }

            ( * this->equationMap ) [ inode ] = * dof2EquationMap;
            delete dof2EquationMap;
        }
    }
}

void
QuasicontinuumNumberingscheme :: reset()
{
    neq = 0;
    pres_neq = 0;
}

int
QuasicontinuumNumberingscheme :: giveDofEquationNumber(Dof *dof) const
{
    //std::map<int,int> map = (equationMap->find(dof->giveDofManNumber())->second)->find((int)dof->giveDofID)->second;


    int dofEqNum = 0;

    if ( selectedNodes.at( dof->giveDofManNumber() ) == 1 ) {
        dofEqNum = ( equationMap->find( dof->giveDofManNumber() )->second ).find( dof->giveDofID() )->second;
    } else {
        dofEqNum = 0;
    }

    //this->nodalEquationNumbers.at( dof->giveDofManNumber() );
    if ( dofEqNum < 0 ) {
        dofEqNum = 0;
    }

    return dofEqNum;
}


int
QuasicontinuumNumberingscheme :: giveTotalNumberOfEquations() const
{
    return neq;
}


int
QuasicontinuumNumberingscheme :: giveRequiredNumberOfDomainEquation() const
{
    return neq;
}

int
QuasicontinuumNumberingscheme :: giveTotalNumberOfPrescribedEquations() const
{
    return -1 * pres_neq;
}
} // end namespace oofem
