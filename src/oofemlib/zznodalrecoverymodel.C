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

#include "zznodalrecoverymodel.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "gausspnt.h"
#include "integrationrule.h"
#include "feinterpol.h"
#include "mathfem.h"
#include "error.h"

#ifdef __PARALLEL_MODE
 #include "communicator.h"
#endif

#define ZZNRM_ZERO_VALUE 1.e-12

namespace oofem {
ZZNodalRecoveryModel :: ZZNodalRecoveryModel(Domain *d) : NodalRecoveryModel(d)
{ }

ZZNodalRecoveryModel :: ~ZZNodalRecoveryModel()
{ }

int
ZZNodalRecoveryModel :: recoverValues(InternalStateType type, TimeStep *tStep)
{
    int ireg, nregions = this->giveNumberOfVirtualRegions();
    int ielem, nelem = domain->giveNumberOfElements();
    int nnodes = domain->giveNumberOfDofManagers();
    int elemNodes;
    int regionValSize;
    int elementNode, node;
    int regionDofMans;
    int i, j, neq, eq;
    Element *element;
    ZZNodalRecoveryModelInterface *interface;
    IntArray skipRegionMap(nregions), regionRecSize(nregions);
    IntArray regionNodalNumbers(nnodes);
    // IntArray loc;
    FloatArray lhs, nn, sol;
    FloatMatrix rhs, nsig;


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
        lhs.resize(regionDofMans);
        lhs.zero();
        rhs.resize(regionDofMans, regionValSize);
        rhs.zero();
        sol.resize(neq);
        sol.zero();
        // assemble element contributions
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            element = domain->giveElement(ielem);

#ifdef __PARALLEL_MODE
            if ( element->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            if ( this->giveElementVirtualRegionNumber(ielem) != ireg ) {
                continue;
            }

            // If an element doesn't implement the interface, it is ignored.
            if ( ( interface = ( ZZNodalRecoveryModelInterface * ) element->giveInterface(ZZNodalRecoveryModelInterfaceType) ) == NULL ) {
                //abort();
                continue;
            }


            // ask element contributions
            interface->ZZNodalRecoveryMI_computeNNMatrix(nn, type);
            interface->ZZNodalRecoveryMI_computeNValProduct(nsig, type, tStep);
            // assemble contributions
            elemNodes = element->giveNumberOfDofManagers();

            //loc.resize ((elemNodes+elemSides)*regionValSize);
            eq = 1;
            for ( elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
                node = element->giveDofManager(elementNode)->giveNumber();
                lhs.at( regionNodalNumbers.at(node) ) += nn.at(eq);
                for ( i = 1; i <= regionValSize; i++ ) {
                    rhs.at(regionNodalNumbers.at(node), i) += nsig.at(eq, i);
                }

                eq++;
            }
        } // end assemble element contributions

#ifdef __PARALLEL_MODE
        this->exchangeDofManValues(ireg, lhs, rhs, regionNodalNumbers);
#endif

        bool missingDofManContribution = false;
        // solve for recovered values of active region
        for ( i = 1; i <= regionDofMans; i++ ) {
            eq = ( i - 1 ) * regionValSize;
            for ( j = 1; j <= regionValSize; j++ ) {
                // rhs will be overriden by recovered values
                if ( fabs( lhs.at(i) ) > ZZNRM_ZERO_VALUE ) {
                    sol.at(eq + j) = rhs.at(i, j) / lhs.at(i);
                } else {
                    missingDofManContribution = true;
                    sol.at(eq + j) = 0.0;
                }
            }
        }

        // update recovered values
        this->updateRegionRecoveredValues(ireg, regionNodalNumbers, regionValSize, sol);

        if ( missingDofManContribution ) {
            OOFEM_WARNING2("ZZNodalRecoveryModel::recoverValues: region %d: values of some dofmanagers undetermined", ireg);
        }
    } // end loop over regions

    this->valType = type;
    this->stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}


void
ZZNodalRecoveryModel :: initRegionMap(IntArray &regionMap, IntArray &regionValSize, InternalStateType type)
{
    int nregions = this->giveNumberOfVirtualRegions();
    int ielem, nelem = domain->giveNumberOfElements();
    int i, regionsSkipped = 0;
    Element *element;
    ZZNodalRecoveryModelInterface *interface;

    regionMap.resize(nregions);
    regionMap.zero();
    regionValSize.resize(nregions);
    regionValSize.zero();

    // loop over elements and check if implement interface
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = domain->giveElement(ielem);
#ifdef __PARALLEL_MODE
        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

#endif
        if ( ( interface =  ( ZZNodalRecoveryModelInterface * ) element->giveInterface(ZZNodalRecoveryModelInterfaceType) ) == NULL ) {
            // If an element doesn't implement the interface, it is ignored.

            //regionsSkipped = 1;
            //regionMap.at( element->giveRegionNumber() ) = 1;
            continue;
        } else {
            int elementVR = this->giveElementVirtualRegionNumber(ielem);
            if ( regionValSize.at(elementVR) ) {
                if ( regionValSize.at(elementVR) != interface->ZZNodalRecoveryMI_giveDofManRecordSize(type) ) {
                    // This indicates a size mis-match between different elements, no choice but to skip the region.
                    regionMap.at(elementVR) = 1;
                    regionsSkipped = 1;
                }
            } else {
                regionValSize.at(elementVR) = interface->ZZNodalRecoveryMI_giveDofManRecordSize(type);
                if ( regionValSize.at(elementVR) == 0 ) {
                    regionMap.at(elementVR) = 1;
                    regionsSkipped = 1;
                    OOFEM_LOG_RELEVANT( "ZZNodalRecoveryModel :: initRegionMap: unknown size of InternalStateType %s\n", __InternalStateTypeToString(type) );
                }
            }
        }
    }

    if ( regionsSkipped ) {
        OOFEM_LOG_RELEVANT( "ZZNodalRecoveryModel :: initRegionMap: skipping regions for InternalStateType %s\n", __InternalStateTypeToString(type) );
        for ( i = 1; i <= nregions; i++ ) {
            if ( regionMap.at(i) ) {
                OOFEM_LOG_RELEVANT("%d ", i);
            }
        }

        OOFEM_LOG_RELEVANT("\n");
    }
}


void
ZZNodalRecoveryModelInterface :: ZZNodalRecoveryMI_computeNValProduct(FloatMatrix &answer, InternalStateType type,
                                                                      TimeStep *tStep)
{  // evaluates N^T sigma over element volume
   // N(nsigma, nsigma*nnodes)
   // Definition : sigmaVector = N * nodalSigmaVector
    int i, j, k, size;
    double dV;
    FloatArray stressVector, help, n;
    Element *elem  = this->ZZNodalRecoveryMI_giveElement();
    IntegrationRule *iRule = elem->giveDefaultIntegrationRulePtr();
    GaussPoint *gp;

    size = ZZNodalRecoveryMI_giveDofManRecordSize(type);
    answer.resize(elem->giveNumberOfDofManagers(), size);

    answer.zero();
    for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        dV  = elem->computeVolumeAround(gp);
        //this-> computeStressVector(stressVector, gp, stepN);
        elem->giveIPValue(stressVector, gp, type, tStep);

        this->ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(n, gp, type);
        for ( j = 1; j <= elem->giveNumberOfDofManagers(); j++ ) {
            for ( k = 1; k <= size; k++ ) {
                answer.at(j, k) += n.at(j) * stressVector.at(k) * dV;
            }
        }

        //  help.beTProductOf(n,stressVector);
        //  answer.add(help.times(dV));
    }
}

void
ZZNodalRecoveryModelInterface :: ZZNodalRecoveryMI_computeNNMatrix(FloatArray &answer, InternalStateType type)
{
    //
    // Returns NTN matrix (lumped) for Zienkiewicz-Zhu
    // The size of N mtrx is (nstresses, nnodes*nstreses)
    // Definition : sigmaVector = N * nodalSigmaVector
    //
    int i, size;
    double sum, dV, volume = 0.0;
    FloatMatrix fullAnswer;
    FloatArray n;
    Element *elem  = this->ZZNodalRecoveryMI_giveElement();
    IntegrationRule *iRule = elem->giveDefaultIntegrationRulePtr();
    GaussPoint *gp;

    size = elem->giveNumberOfDofManagers(); //ZZNodalRecoveryMI_giveDofManRecordSize (type);
    fullAnswer.resize(size, size);
    fullAnswer.zero();
    double pok = 0.0;

    for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        dV  = elem->computeVolumeAround(gp);
        this->ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(n, gp, type);
        fullAnswer.plusDyadSymmUpper(n, n, dV);
        pok += ( n.at(1) * dV ); ///@todo What is this? Completely unused.
        volume += dV;
    }


    fullAnswer.symmetrized();
    answer.resize(size);
    for ( i = 1; i <= size; i++ ) {
        sum = 0.0;
        for ( int j = 1; j <= size; j++ ) {
            sum += fullAnswer.at(i, j);
        }

        answer.at(i) = sum;
    }
}

void
ZZNodalRecoveryModelInterface :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *gp, InternalStateType type)
{
    // evaluates N matrix (interpolation estimated stress matrix)
    // according to Zienkiewicz & Zhu paper
    // N(nsigma, nsigma*nnodes)
    // Definition : sigmaVector = N * nodalSigmaVector

    Element* elem = ZZNodalRecoveryMI_giveElement();
    FEInterpolation* interpol = elem->giveInterpolation();

    // test if underlying element provides interpolation
    if (interpol) {
        if ( !elem->giveIPValueSize(type, gp) ) {
            OOFEM_ERROR3("ZZNodalRecoveryMI_computeNNMatrix: Element %d not supporting type %d", elem->giveNumber(), type);
            return;
        }

        interpol->evalN(answer, *gp->giveCoordinates(), FEIElementGeometryWrapper(elem));
    } else {
        // ok default implementation can not work, as element is not providing valid interpolation
        // to resolve this, one can overload this method for element implementing ZZNodalRecoveryModelInterface
        // or element should provide interpolation.
        OOFEM_ERROR2("ZZNodalRecoveryMI_computeNNMatrix: Element %d not providing valid interpolation", elem->giveNumber());
    }
}

int
ZZNodalRecoveryModelInterface :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    IntegrationRule* iRule = ZZNodalRecoveryMI_giveElement()->giveDefaultIntegrationRulePtr();
    if (iRule) {
        return ZZNodalRecoveryMI_giveElement()->giveIPValueSize(type, iRule->getIntegrationPoint(0));
    } else {
        return 0;
    }
}



#ifdef __PARALLEL_MODE

void
ZZNodalRecoveryModel :: initCommMaps()
{
 #ifdef __PARALLEL_MODE
    if ( initCommMap ) {
        EngngModel *emodel = domain->giveEngngModel();
        ProblemCommunicatorMode commMode = emodel->giveProblemCommMode();
        if ( commMode == ProblemCommMode__NODE_CUT ) {
            commBuff = new CommunicatorBuff(emodel->giveNumberOfProcesses(), CBT_dynamic);
            communicator = new ProblemCommunicator(emodel, commBuff, emodel->giveRank(),
                                                   emodel->giveNumberOfProcesses(),
                                                   commMode);
            communicator->setUpCommunicationMaps(domain->giveEngngModel(), true, true);
            OOFEM_LOG_INFO("ZZNodalRecoveryModel :: initCommMaps: initialized comm maps");
            initCommMap = false;
        } else {
            OOFEM_ERROR("ZZNodalRecoveryModel :: initCommMaps: unsupported comm mode");
        }
    }

 #endif
}

void
ZZNodalRecoveryModel :: exchangeDofManValues(int ireg, FloatArray &lhs, FloatMatrix &rhs, IntArray &rn)
{
    EngngModel *emodel = domain->giveEngngModel();
    ProblemCommunicatorMode commMode = emodel->giveProblemCommMode();

    if ( commMode == ProblemCommMode__NODE_CUT ) {
        parallelStruct ls( &lhs, &rhs, &rn);

        // exchange data for shared nodes
        communicator->packAllData(this, & ls, & ZZNodalRecoveryModel :: packSharedDofManData);
        communicator->initExchange(789 + ireg);
        communicator->unpackAllData(this, & ls, & ZZNodalRecoveryModel :: unpackSharedDofManData);
        communicator->finishExchange();
    } else {
        OOFEM_ERROR("ZZNodalRecoveryModel :: exchangeDofManValues: Unsupported commMode");
    }
}

int
ZZNodalRecoveryModel :: packSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm)
{
    int result = 1, i, j, indx, nc, size;
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    IntArray const *toSendMap = processComm.giveToSendMap();
    nc = s->rhs->giveNumberOfColumns();

    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        // toSendMap contains all shared dofmans with remote partition
        // one has to check, if particular shared node value is available for given region
        indx = s->regionNodalNumbers->at( toSendMap->at(i) );
        if ( indx ) {
            // pack "1" to indicate that for given shared node this is a valid contribution
            result &= pcbuff->packInt(1);
            result &= pcbuff->packDouble( s->lhs->at(indx) );
            for ( j = 1; j <= nc; j++ ) {
                result &= pcbuff->packDouble( s->rhs->at(indx, j) );
            }

            //printf("[%d] ZZ: Sending data for shred node %d[%d]\n", domain->giveEngngModel()->giveRank(),
            //       toSendMap->at(i), domain->giveDofManager(toSendMap->at(i))->giveGlobalNumber());
        } else {
            // ok shared node is not in active region (determined by s->regionNodalNumbers)
            result &= pcbuff->packInt(0);
        }
    }

    return result;
}

int
ZZNodalRecoveryModel :: unpackSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm)
{
    int result = 1;
    int i, j, nc, indx, size, flag;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    double value;
    nc = s->rhs->giveNumberOfColumns();

    size = toRecvMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        indx = s->regionNodalNumbers->at( toRecvMap->at(i) );
        // toRecvMap contains all shared dofmans with remote partition
        // one has to check, if particular shared node received contribution is available for given region
        result &= pcbuff->unpackInt(flag);
        if ( flag ) {
            // "1" to indicates that for given shared node this is a valid contribution
            result &= pcbuff->unpackDouble(value);
            // now check if we have a valid number
            if ( indx ) {
                s->lhs->at(indx) += value;
            }

            for ( j = 1; j <= nc; j++ ) {
                result &= pcbuff->unpackDouble(value);
                if ( indx ) {
                    s->rhs->at(indx, j) += value;
                }
            }

            //if (indx) printf("[%d] ZZ: Receiving data for shred node %d[%d]\n", domain->giveEngngModel()->giveRank(),
            //                 toRecvMap->at(i), domain->giveDofManager(toRecvMap->at(i))->giveGlobalNumber());
        }
    }

    return result;
}

#endif
} // end namespace oofem
