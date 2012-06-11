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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef __MAKEDEPEND
 #include <cstdio>
 #include <cstdlib>
 #include <set>
#endif

#include "freestor.h"
#include "domain.h"
#include "element.h"
#include "dofmanager.h"
#include "conTable.h"

namespace oofem {
ConnectivityTable :: ~ConnectivityTable()
// destructor
{
    this->reset();
}

void
ConnectivityTable :: reset()
{
    if ( dofManagersConnectivity ) {
        freeInt(dofManagersConnectivity);
    }

    dofManagersConnectivity = NULL;
    nodalConnectivityFlag = 0;
}

#ifdef __OOFEG
/*****
 * void ConnectivityTable ::  allocateConnectivityTable ()
 * //
 * // creates connectivity table
 * // (finds for each node, to which elements  belongs)
 * //
 * // connectivity[inode] = nrec - contains number of elements per node
 * // which have stress corresponding to stressMode
 * //
 * {
 * if ((dofManagersConnectivity != NULL))
 * { this->resetYourself(); return;} // previously crated
 *
 * int nnodes = domain -> giveNumberOfDofManagers ();
 * dofManagersConnectivity = allocInt (nnodes);
 * dofManagersValues = allocDouble (nnodes );
 *
 * this-> resetYourself ();
 * }
 *
 *
 * void ConnectivityTable :: instanciateYourself(oofegGraphicContext& gc)
 * //
 * // instanciates yoursef - ie build node table
 * //
 * {
 * Element *elem;
 * int i,j,nelems, enodes,ires,inode;
 * double stress;
 * DrawMode mode = gc.getDrawMode();
 *
 * // ensure allocated table
 * this->allocateConnectivityTable ();
 * // start loop over elements - node part
 * nelems = domain-> giveNumberOfElements();
 * for (i = 1; i<= nelems; i++) {
 *  elem = domain -> giveElement(i);
 *  enodes = elem ->giveNumberOfNodes();
 *  // loop over each element node
 *  for (j = 1; j<= enodes; j++) {
 *    inode = elem -> giveNode (j)->giveNumber();
 *    ires  = elem -> giveInternalStateAtNode (gc,j,&stress);
 *    // store inode into connectivity table
 *    // check for free space
 *    if (ires != 0) {
 *  if (mode == stressErrorState) {
 *  // store
 *  dofManagersConnectivity[inode-1] ++;
 *  if (stress > dofManagersValues [inode-1])   dofManagersValues [inode-1] = stress;
 *  } else {
 *  // store
 *  dofManagersConnectivity[inode-1] ++;
 *  dofManagersValues [inode-1]      += stress;
 *  }
 *    } else {
 *  break;
 *    };
 *  }
 * }
 *
 * // side part
 * for (i = 1; i<= nelems; i++) {
 *  elem = domain -> giveElement(i);
 *  enodes = elem ->giveNumberOfSides();
 *  // loop over each element node
 *  for (j = 1; j<= enodes; j++) {
 *    inode = elem -> giveSide (j)->giveNumber();
 *    ires  = elem -> giveInternalStateAtSide (gc, j, &stress);
 *    // store inode into connectivity table
 *    // check for free space
 *    if (ires != 0) {
 *  if (mode == stressErrorState) {
 *  // store
 *  dofManagersConnectivity[inode-1] ++;
 *  if (stress > dofManagersValues [inode-1])   dofManagersValues [inode-1] = stress;
 *  } else {
 *  // store
 *  dofManagersConnectivity[inode-1] ++;
 *  dofManagersValues [inode-1]      += stress;
 *  }
 *    } else {
 *  break;
 *    };
 *  }
 * }
 *
 * if (mode != stressErrorState) {
 * int nman = domain -> giveNumberOfDofManagers ();
 * for (i = 0; i< nman; i++) {
 * if (dofManagersConnectivity[i] != 0)
 *  dofManagersValues[i] /= dofManagersConnectivity[i];
 * }
 * }
 * currentMode = mode;
 * }
 *
 *
 *
 * int ConnectivityTable :: giveNumberOfElementsInDofMngr (int i)
 * //
 * //returning number of Elements belonging to dof Managar
 * //
 * {
 * if (dofManagersConnectivity != NULL) error("giveNumberOfElementsInNode : not initialized");
 * return dofManagersConnectivity[i-1];
 * }
 *
 *
 * double ConnectivityTable :: giveDofMngrValue(oofegGraphicContext& gc, int i)
 * //
 * // returning j th element belonging to node i
 * //
 * {
 * if (currentMode != gc.getDrawMode() ) instanciateYourself(gc);
 * if (dofManagersConnectivity == NULL) error("giveNodeValue : not initialized");
 * return dofManagersValues [i-1];
 * }
 *
 *
 * void ConnectivityTable :: giveMinMaxVal (oofegGraphicContext& gc,double *min, double *max)
 * //
 * // returns min and max value for given mode
 * // initializes for mode also
 * //
 * {
 * double val;
 * int i;
 *
 * this->instanciateYourself(gc);
 * int nnodes = domain -> giveNumberOfDofManagers ();
 **min = +1.e30;
 **max = -1.e30;
 * for (i = 0; i< nnodes; i++) {
 * if (dofManagersConnectivity[i] != 0) {
 * val = dofManagersValues[i];
 * if (*min > val) *min = val;
 * if (*max < val) *max = val;
 * }
 * }
 *
 * return ;
 * }
 *************/
#endif


/*
 * #include "dof.h"
 * #include "structuralelement.h"
 *
 * void
 * ConnectivityTable :: instanciateReactionForceTable (TimeStep* tStep)
 * //
 * // assembles table of restrained nodes and elementts connected to this nodes
 * //
 * {
 *
 * // int ndofMan = domain -> giveNumberOfDofManagers();
 * // int nelems = domain -> giveNumberOfElements();
 * // DofManager* inode;
 * // Element* ielem;
 * // int i,j,k,kk,jnode,restr,indofs,nnodes,nChunks = 1;
 * // int *help;
 * //
 * //
 * // if (ReactionTableInit) return; // already initialized
 * //
 * // restrDofManager = allocInt (nChunks * CHUNK + 1);
 * // // first find nodes which are restrained
 * // for (i = 1; i<= ndofMan; i++) {
 * //  restr = 0;
 * //  inode = domain->giveDofManager(i);
 * //  indofs = inode -> giveNumberOfDofs();
 * //  for (j=1; j<= indofs; j++)
 * //   if (inode->giveDof(j)->hasBc(tStep)) restr = 1;
 * //  if (restr) {
 * //   numRestrManagers ++;
 * //   if (numRestrManagers > nChunks * CHUNK) {
 * //    // realocate table and copy to new place
 * //    help = allocInt((nChunks+1) * CHUNK + 1);
 * //    for (j=1 ; j<= nChunks * CHUNK; j++) help[j]=restrDofManager[j];
 * //    freeInt (restrDofManager);
 * //    restrDofManager = help;
 * //    nChunks ++;
 * //   }
 * //   restrDofManager[numRestrManagers] = i;
 * //  }
 * // }
 * //
 * // if (numRestrManagers != 0) {
 * //  // create coonectivity table for elements sharing
 * //  // restrained nodes
 * //  numOfSharedElements = new IntArray (numRestrManagers);
 * //  sharedElements = (int**) malloc (sizeof(int*)*(numRestrManagers+1));
 * //
 * //  for( kk =1; kk<=numRestrManagers; kk++)
 * //   sharedElements[kk] = (int *) allocInt ((MAX_CONNECTIVITY+1));
 * //
 * //  for (i=1 ; i<= nelems; i++) {
 * //   ielem = domain->giveElement(i);
 * //   nnodes = ielem -> giveNumberOfNodes();
 * //   for (j=1; j<= nnodes; j++) {
 * //    jnode = ielem->giveNode (j)->giveNumber();
 * //    // find if jnode is in restrNodes table
 * //    for (k=1; k<= numRestrManagers; k++)
 * //     if (jnode == restrDofManager[k]) {
 * //      if (numOfSharedElements->at(k) >= MAX_CONNECTIVITY) {
 * //       freeInt (restrDofManager);
 * //       for(kk =1; kk<=numRestrManagers; kk++) {freeInt (sharedElements[kk]);}
 * //       free  (sharedElements);
 * //       delete numOfSharedElements;
 * //       error ("printReactionForces: MAX_Connectivity reached, reactions are not computed");
 * //       return;
 * //      }
 * //      numOfSharedElements->at(k) ++;
 * //      sharedElements[k][numOfSharedElements->at(k)] = i;
 * //     }
 * //   }
 * //  }
 * //
 * //  for (i=1 ; i<= nelems; i++) {
 * //   ielem = domain->giveElement(i);
 * //   nnodes = ielem -> giveNumberOfSides();
 * //   for (j=1; j<= nnodes; j++) {
 * //    jnode = ielem->giveSide (j)->giveNumber();
 * //    // find if jnode is in restrNodes table
 * //    for (k=1; k<= numRestrManagers; k++)
 * //     if (jnode == restrDofManager[k]) {
 * //      if (numOfSharedElements->at(k) >= MAX_CONNECTIVITY) {
 * //       freeInt (restrDofManager);
 * //       for(kk =1; kk<=numRestrManagers; kk++) {freeInt (sharedElements[kk]);}
 * //       free  (sharedElements);
 * //       delete numOfSharedElements;
 * //       error ("printReactionForces: MAX_Connectivity reached, reactions are not computed");
 * //       return;
 * //      }
 * //      numOfSharedElements->at(k) ++;
 * //      sharedElements[k][numOfSharedElements->at(k)] = i;
 * //     }
 * //   }
 * //  }
 * // }
 * // ReactionTableInit = 1;
 * //
 *
 * int ndofMan = domain -> giveNumberOfDofManagers();
 * DofManager* inode;
 * //Element* ielem;
 * Dof * jdof;
 * int i,j, indofs, iRestrDof;
 *
 * this->numRestrDofs = 0;
 *
 * // determine number of restrained dofs
 * for (i=1; i<= ndofMan; i++) {
 * inode = domain->giveDofManager(i);
 * indofs = inode -> giveNumberOfDofs();
 * for (j=1; j<= indofs; j++) {
 * jdof = inode->giveDof(j);
 * if ((jdof->giveClassID() != SimpleSlaveDofClass) && (jdof-> hasBc (tStep))) this->numRestrDofs++; // skip slave dofs
 * }
 * }
 *
 * // initialize corresponding dofManagers and dofs for each restrained dof
 * this->restrDofManager.resize (this->numRestrDofs);
 * this->restrDof.resize (this->numRestrDofs);
 *
 * iRestrDof = 0;
 * for (i=1; i<= ndofMan; i++) {
 * inode = domain->giveDofManager(i);
 * indofs = inode -> giveNumberOfDofs();
 * for (j=1; j<= indofs; j++) {
 * jdof = inode->giveDof(j);
 * if ((jdof->giveClassID() != SimpleSlaveDofClass) && (jdof-> hasBc (tStep))) { // skip slave dofs
 *  iRestrDof ++;
 *  this->restrDofManager.at(iRestrDof) = i;
 *  this->restrDof.at(iRestrDof) = j;
 * }
 * }
 * }
 * return;
 * }
 *
 * int
 * ConnectivityTable::giveReactionIndx (int dofMan, int dof) const
 * {
 * //
 * //   return corresponding reaction index, zero if error
 *
 * int size = this->giveNumRestrDofs();
 * if (size == 0) return 0;
 * if (size == 1) {
 * if ((dofMan == giveReactionDofManager(1)) && (dof == giveReactionDof(1))) return 1;
 * else return 0;
 * }
 *
 * int indx = 0;
 * int startIndx, endIndx, middleIndx, middleIndxDofMan;
 * startIndx = 1;
 * endIndx = size;
 *
 * if (dofMan == giveReactionDofManager (startIndx)) {
 * indx = startIndx;
 * } else  if (dofMan == giveReactionDofManager (endIndx)) {
 * indx = endIndx;
 * } else {
 *
 * while (startIndx!= endIndx) {
 * middleIndx = (endIndx+startIndx)/2;
 * middleIndxDofMan = giveReactionDofManager (middleIndx);
 * if (dofMan == middleIndxDofMan) {
 *  indx = middleIndx;
 *  break;
 * } else if (dofMan < middleIndxDofMan) {
 *  endIndx = middleIndx;
 * } else {
 *  if (startIndx != middleIndx) startIndx = middleIndx; else startIndx ++;
 * }
 * }
 * }
 *
 * int ii;
 *
 * if (indx == 0) return 0; // reaction entry not found
 * if (dof == giveReactionDof (indx)) { // found!
 * return indx;
 * } else if (dof < giveReactionDof (indx)) {
 * ii = -1;
 * } else {
 * ii = 1;
 * }
 *
 * do {
 * indx += ii;
 * if ((indx == 0) || (indx == size+1)) return 0;
 * if (dofMan != giveReactionDofManager (indx)) return 0; // dof entry for dofManager not found!
 * } while (dof != giveReactionDof (indx));
 * return indx;
 * }
 */

/*
 * int
 * ConnectivityTable :: giveNumOfSharedElements (int inode)
 * //
 * // returns number of shared elements for restrained node inode
 * //
 * {
 * if (numOfSharedElements)
 * return numOfSharedElements->at(inode);
 * return 0;
 * }
 *
 * int
 * ConnectivityTable :: giveSharedElement (int inode, int ielem)
 * //
 * // returns ielem-th element shared by inode-th restrained node
 * //
 * {
 * if (ReactionTableInit) {
 * return sharedElements[inode][ielem];
 * }
 * return 0;
 * }
 */

void ConnectivityTable :: printYourself()
//
// print Yourself
//
{
    /*
     * printf ("\nConnectivityTable {");
     * if (dofManagersConnectivity) {
     * int nnodes = domain -> giveNumberOfDofManagers ();
     * for (int i = 0; i < nnodes; i++) {
     *  printf ("\tDoofManager %d {value: %f}\n",i+1,dofManagersValues[i]);
     * }
     * printf("}\n");
     * }
     */
    /*
     * if (ReactionTableInit) {
     * printf("\tnumberOfRestrainedDofManagers is %d\n",giveNumRestrDofMngrs());
     * int nnodes = giveNumRestrDofMngrs();
     * int nelems;
     *
     * for (int i=1; i <= nnodes; i++) {
     *  printf("\n\t%d th restrained dofManager %d,  connected elements:\n\t{",
     *     i,giveRestrDofManager(i));
     *  nelems = giveNumOfSharedElements(i);
     *  for (int j=1; j<=nelems;j++) {
     *   printf(" %d", giveSharedElement(i,j));
     *  }
     *  printf(" }");
     * }
     * }
     */
}


#ifdef __OOFEG
/****
 * void ConnectivityTable :: resetYourself ()
 * //
 * // reset tables
 * //
 * {
 * int nnodes = domain -> giveNumberOfDofManagers ();
 * for (int i= 0; i< nnodes; i++) {
 *  dofManagersConnectivity[i]=0;
 *  dofManagersValues[i]=0.;
 * }
 *
 * currentMode = unknown;
 * }
 ******/
#endif




void
ConnectivityTable :: instanciateConnectivityTable()
//
// assembles table of nodal connectivities
//
{
    int ndofMan = domain->giveNumberOfDofManagers();
    int nelems = domain->giveNumberOfElements();
    int i, j, jnode, nnodes;
    Element *ielem;
    IntArray dofManConnectivity(ndofMan);

    if ( nodalConnectivityFlag ) {
        return;                     // already initialized
    }

    OOFEM_LOG_INFO("ConnectivityTable: initializing\n");

    for ( i = 1; i <= nelems; i++ ) {
        ielem = domain->giveElement(i);
        nnodes = ielem->giveNumberOfDofManagers();
        for ( j = 1; j <= nnodes; j++ ) {
            jnode = ielem->giveDofManager(j)->giveNumber();
            dofManConnectivity.at(jnode)++;
        }
    }

    // allocate Nodal connectivity table for domain
    nodalConnectivity.growTo(ndofMan);
    for ( i = 1; i <= ndofMan; i++ ) {
        nodalConnectivity.put( i, new IntArray( dofManConnectivity.at(i) ) );
    }

    // build Nodal connectivity table for domain
    dofManConnectivity.zero();

    for ( i = 1; i <= nelems; i++ ) {
        ielem = domain->giveElement(i);
        ielem = domain->giveElement(i);
        nnodes = ielem->giveNumberOfDofManagers();
        for ( j = 1; j <= nnodes; j++ ) {
            jnode = ielem->giveDofManager(j)->giveNumber();
            nodalConnectivity.at(jnode)->at( ++dofManConnectivity.at(jnode) ) = i;
        }
    }

    nodalConnectivityFlag = 1;
}

const IntArray *
ConnectivityTable :: giveDofManConnectivityArray(int dofman)
{
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    return this->nodalConnectivity.at(dofman);
}


void
ConnectivityTable :: giveElementNeighbourList(IntArray &answer, IntArray &elemList)
{
    int i, j, k, nnode, jnode, nelems = elemList.giveSize();
    Element *ielem;
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    std :: set< int >neighbours;

    for ( i = 1; i <= nelems; i++ ) {
        ielem = domain->giveElement( elemList.at(i) );
        nnode = ielem->giveNumberOfDofManagers();
        for ( j = 1; j <= nnode; j++ ) {
            jnode = ielem->giveDofManager(j)->giveNumber();
            for ( k = 1; k <= this->nodalConnectivity.at(jnode)->giveSize(); k++ ) {
                neighbours.insert( this->nodalConnectivity.at(jnode)->at(k) );
            }
        }
    }

    answer.resize( neighbours.size() );
    std :: set< int > :: iterator pos;
    for ( pos = neighbours.begin(), i = 1; pos != neighbours.end(); ++pos, i++ ) {
        answer.at(i) = * pos;
    }
}


void
ConnectivityTable :: giveNodeNeighbourList(IntArray &answer, IntArray &nodeList)
{
    int i, k, inode, nnodes = nodeList.giveSize();
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    std :: set< int >neighbours;

    for ( i = 1; i <= nnodes; i++ ) {
        inode = nodeList.at(i);
        for ( k = 1; k <= this->nodalConnectivity.at(inode)->giveSize(); k++ ) {
            neighbours.insert( this->nodalConnectivity.at(inode)->at(k) );
        }
    }

    answer.resize( neighbours.size() );
    std :: set< int > :: iterator pos;
    for ( pos = neighbours.begin(), i = 1; pos != neighbours.end(); ++pos, i++ ) {
        answer.at(i) = * pos;
    }
}
} // end namespace oofem
