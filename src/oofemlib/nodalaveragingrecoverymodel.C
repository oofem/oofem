/* $Header: /home/cvs/bp/oofem/oofemlib/src/nodalaveragingrecoverymodel.C,v 1.2.4.1 2004/04/05 15:19:43 bp Exp $ */
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

//
// file nodalaveragingrecoverymodel.C
//

#include "nodalaveragingrecoverymodel.h"
#include "timestep.h"
#include "element.h"
#include "node.h"
#include "elementside.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

namespace oofem {

NodalAveragingRecoveryModel :: NodalAveragingRecoveryModel(Domain *d) : NodalRecoveryModel(d)
{ }

NodalAveragingRecoveryModel :: ~NodalAveragingRecoveryModel()
{ }

int
NodalAveragingRecoveryModel :: recoverValues(InternalStateType type, TimeStep *tStep)
{
    int ireg, nregions = domain->giveNumberOfRegions();
    int ielem, nelem = domain->giveNumberOfElements();
    int inode, nnodes = domain->giveNumberOfDofManagers();
    int elemNodes;
    int regionValSize;
    int elementNode, node;
    int regionDofMans;
    int i, neq, eq;
    Element *element;
    NodalAveragingRecoveryModelInterface *interface;
    IntArray skipRegionMap(nregions), regionRecSize(nregions);
    IntArray regionNodalNumbers(nnodes);
    IntArray regionDofMansConnectivity;
    FloatArray lhs, val;


    if ( ( this->valType == type ) && ( this->stateCounter == tStep->giveSolutionStateCounter() ) ) {
        return 1;
    }

#ifdef __PARALLEL_MODE
    this->initCommMaps();
#endif

    // clear nodal table
    this->clear();

    // init region table indicating regions to skip
    this->initRegionMap(skipRegionMap, regionRecSize, type);

#ifdef __PARALLEL_MODE
    // synchronize skipRegionMap over all cpus
    IntArray temp_skipRegionMap(skipRegionMap);
    MPI_Allreduce(temp_skipRegionMap.givePointer(), skipRegionMap.givePointer(), nregions, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
#endif

    // loop over regions
    for ( ireg = 1; ireg <= nregions; ireg++ ) {
        // skip regions
        if ( skipRegionMap.at(ireg) ) {
            continue;
        }

        // loop over elements and determine local region node numbering and determine and check nodal values size
        if ( this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, ireg) == 0 ) {
            break;
        }

        regionValSize = regionRecSize.at(ireg);
        neq = regionDofMans * regionValSize;
        lhs.resize(neq);
        lhs.zero();
        regionDofMansConnectivity.resize(regionDofMans);
        regionDofMansConnectivity.zero();

        // assemble element contributions
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            element = domain->giveElement(ielem);

#ifdef __PARALLEL_MODE
            if (element->giveParallelMode() != Element_local) continue;
#endif
            if ( element->giveRegionNumber() != ireg ) {
                continue;
            }

            if ( ( interface = ( NodalAveragingRecoveryModelInterface * )
                               element->giveInterface(NodalAveragingRecoveryModelInterfaceType) ) == NULL ) {
                abort();
            }

            elemNodes = element->giveNumberOfDofManagers();
            // ask element contributions
            for ( elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
                node = element->giveDofManager(elementNode)->giveNumber();
                interface->NodalAveragingRecoveryMI_computeNodalValue(val, elementNode, type, tStep);
                eq = ( regionNodalNumbers.at(node) - 1 ) * regionValSize;
                for ( i = 1; i <= regionValSize; i++ ) {
                    lhs.at(eq + i) += val.at(i);
                }

                regionDofMansConnectivity.at( regionNodalNumbers.at(node) )++;
            }
        } // end assemble element contributions

#ifdef __PARALLEL_MODE
        this->exchangeDofManValues(ireg, lhs, regionDofMansConnectivity, regionNodalNumbers, regionValSize);
#endif

        // solve for recovered values of active region
        for ( inode = 1; inode <= nnodes; inode++ ) {
            if ( regionNodalNumbers.at(inode) ) {
                eq = ( regionNodalNumbers.at(inode) - 1 ) * regionValSize;
                for ( i = 1; i <= regionValSize; i++ ) {
                    lhs.at(eq + i) /= regionDofMansConnectivity.at( regionNodalNumbers.at(inode) );
                }
            }
        }

        // update recovered values
        this->updateRegionRecoveredValues(ireg, regionNodalNumbers, regionValSize, lhs);
    } // end loop over regions

    this->valType = type;
    this->stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}

void
NodalAveragingRecoveryModel :: initRegionMap(IntArray &regionMap, IntArray &regionValSize, InternalStateType type)
{
    int nregions = domain->giveNumberOfRegions();
    int ielem, nelem = domain->giveNumberOfElements();
    int i, regionsSkipped = 0;
    Element *element;
    NodalAveragingRecoveryModelInterface *interface;

    regionMap.resize(nregions);
    regionMap.zero();
    regionValSize.resize(nregions);
    regionValSize.zero();

    // loop over elements and check if implement interface
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = domain->giveElement(ielem);

#ifdef __PARALLEL_MODE
        if (element->giveParallelMode() != Element_local) continue;
#endif

        if ( ( interface =  ( NodalAveragingRecoveryModelInterface * ) element->
                           giveInterface(NodalAveragingRecoveryModelInterfaceType) ) == NULL ) {
            /*
             * printf ("NodalAveragingRecoveryModel :: initRegionMap: Element %d does not support required interface", ielem);
             * printf ("-skipping region %d\n", element->giveRegionNumber());
             */
            regionsSkipped = 1;
            regionMap.at( element->giveRegionNumber() ) = 1;
            continue;
        } else {
            if ( regionValSize.at( element->giveRegionNumber() ) ) {
                if ( regionValSize.at( element->giveRegionNumber() ) !=
                    interface->NodalAveragingRecoveryMI_giveDofManRecordSize(type) ) {
                    regionMap.at( element->giveRegionNumber() ) = 1;
                    /*
                     *       printf ("NodalAveragingRecoveryModel :: initRegionMap: element %d has incompatible value size, skipping region\n",
                     *               ielem);
                     */
                    regionsSkipped = 1;
                }
            } else {
                regionValSize.at( element->giveRegionNumber() ) = interface->
                                                                  NodalAveragingRecoveryMI_giveDofManRecordSize(type);
            }
        }
    }

    if ( regionsSkipped ) {
        OOFEM_LOG_RELEVANT("NodalAveragingRecoveryModel :: initRegionMap: skipping regions ");
        for ( i = 1; i <= nregions; i++ ) {
            if ( regionMap.at(i) ) {
                OOFEM_LOG_RELEVANT("%d ", i);
            }
        }

        OOFEM_LOG_RELEVANT("\n");
    }
}


#ifdef __PARALLEL_MODE

void
NodalAveragingRecoveryModel :: initCommMaps ()
{
#ifdef __PARALLEL_MODE
  if (initCommMap) {
    EngngModel *emodel=domain->giveEngngModel();
    ProblemCommunicatorMode commMode = emodel->giveProblemCommMode();
    if (commMode == ProblemCommMode__NODE_CUT) {
      commBuff = new CommunicatorBuff(emodel->giveNumberOfProcesses(), CBT_dynamic);
      communicator = new ProblemCommunicator(emodel, commBuff, emodel->giveRank(),
                                             emodel->giveNumberOfProcesses(),
                                             commMode);
      communicator->setUpCommunicationMaps(domain->giveEngngModel(), true, true);
      OOFEM_LOG_INFO ("NodalAveragingRecoveryModel :: initCommMaps: initialized comm maps");
      initCommMap = false;
    } else {
      OOFEM_ERROR ("NodalAveragingRecoveryModel :: initCommMaps: unsupported comm mode");
    }
  }
#endif
}

void
NodalAveragingRecoveryModel :: exchangeDofManValues (int ireg, FloatArray& lhs, IntArray& regionDofMansConnectivity,
                                                     IntArray& regionNodalNumbers, int regionValSize) 
{
  EngngModel *emodel = domain->giveEngngModel();
  ProblemCommunicatorMode commMode = emodel->giveProblemCommMode();
  
  if (commMode == ProblemCommMode__NODE_CUT) {
    parallelStruct ls (&lhs, &regionDofMansConnectivity, &regionNodalNumbers, regionValSize);
    
    // exchange data for shared nodes
    communicator->packAllData( this, &ls, &NodalAveragingRecoveryModel::packSharedDofManData );
    communicator->initExchange(789+ireg);
    communicator->unpackAllData( this, &ls, &NodalAveragingRecoveryModel::unpackSharedDofManData );
    communicator->finishExchange();
  } else {
    OOFEM_ERROR ("NodalAveragingRecoveryModel :: exchangeDofManValues: Unsupported commMode");
  }
}

int
NodalAveragingRecoveryModel :: packSharedDofManData (parallelStruct* s, ProcessCommunicator &processComm)
{
  int result = 1, i, j, indx, eq, size;
  ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
  IntArray const *toSendMap = processComm.giveToSendMap();
  
  size = toSendMap->giveSize();
  for ( i = 1; i <= size; i++ ) {
    // toSendMap contains all shared dofmans with remote partition
    // one has to check, if particular shared node value is available for given region
    indx = s->regionNodalNumbers->at(toSendMap->at(i));
    if (indx) {
      // pack "1" to indicate that for given shared node this is a valid contribution
      result &= pcbuff->packInt (1);
      result &= pcbuff->packInt(s->regionDofMansConnectivity->at(indx));
      eq = ( indx - 1 ) * s->regionValSize;
      for (j=1; j<=s->regionValSize; j++) result &= pcbuff->packDouble (s->lhs->at(eq+j));
    } else {
      // ok shared node is not in active region (determined by s->regionNodalNumbers)
      result &= pcbuff->packInt (0);
    }
  } 
  return result;
  
}

int 
NodalAveragingRecoveryModel :: unpackSharedDofManData (parallelStruct* s, ProcessCommunicator &processComm)
{
    int result = 1;
    int i, j, eq, indx, size, flag, intValue;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    double value;

    size = toRecvMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
      indx = s->regionNodalNumbers->at(toRecvMap->at(i));
      // toRecvMap contains all shared dofmans with remote partition
      // one has to check, if particular shared node received contribution is available for given region
      result &= pcbuff->unpackInt (flag);
      if (flag) {
        // "1" to indicates that for given shared node this is a valid contribution
        result &= pcbuff->unpackInt(intValue);
        // now check if we have a valid number
        if (indx) s->regionDofMansConnectivity->at(indx) += intValue;
        eq = ( indx - 1 ) * s->regionValSize;
        for (j=1; j<=s->regionValSize; j++) {
          result &= pcbuff->unpackDouble(value);
          if (indx) s->lhs->at(eq+j) += value;
        } 
      } 
    }
    return result;
}

#endif

} // end namespace oofem
