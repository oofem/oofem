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

#include "zznodalrecoverymodel.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "gausspoint.h"
#include "integrationrule.h"
#include "feinterpol.h"
#include "mathfem.h"
#include "feinterpol.h"
#include "error.h"
#include "engngm.h"

#include <sstream>
#include <set>

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "communicator.h"
#endif

#define ZZNRM_ZERO_VALUE 1.e-12

namespace oofem {
ZZNodalRecoveryModel :: ZZNodalRecoveryModel(Domain *d) : NodalRecoveryModel(d)
{ }

ZZNodalRecoveryModel :: ~ZZNodalRecoveryModel()
{ }

int
ZZNodalRecoveryModel :: recoverValues(Set elementSet, InternalStateType type, TimeStep *tStep)
{
    int nnodes = domain->giveNumberOfDofManagers();
    IntArray regionNodalNumbers(nnodes);
    // following variable is for better error reporting only
    std :: set< int >unresolvedDofMans;
    // IntArray loc;
    FloatArray lhs, nn, sol;
    FloatMatrix rhs, nsig;


    if ( this->valType == type && this->stateCounter == tStep->giveSolutionStateCounter() ) {
        return 1;
    }

#ifdef __PARALLEL_MODE
    if ( this->domain->giveEngngModel()->isParallel() ) {
        this->initCommMaps();
    }
#endif

    // clear nodal table
    this->clear();

    int elemNodes;
    int regionValSize;
    int regionDofMans;

    // loop over elements and determine local region node numbering and determine and check nodal values size
    if ( this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, elementSet) == 0 ) {
        return 0;
    }

    regionValSize = 0;
    lhs.resize(regionDofMans);
    lhs.zero();
    IntArray elements = elementSet.giveElementList();
    // assemble element contributions
    for ( int i = 1; i <= elements.giveSize(); i++ ) {
        int ielem = elements.at(i);
        ZZNodalRecoveryModelInterface *interface;
        Element *element = domain->giveElement(ielem);

        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

        // If an element doesn't implement the interface, it is ignored.
        if ( ( interface = static_cast< ZZNodalRecoveryModelInterface * >( element->giveInterface(ZZNodalRecoveryModelInterfaceType) ) ) == NULL ) {
            //abort();
            continue;
        }


        // ask element contributions
        if (!interface->ZZNodalRecoveryMI_computeNValProduct(nsig, type, tStep)) {
          // skip element contribution if value type recognized by element
          continue;
        }
        interface->ZZNodalRecoveryMI_computeNNMatrix(nn, type);

        // assemble contributions
        elemNodes = element->giveNumberOfDofManagers();

        if ( regionValSize == 0 ) {
            regionValSize = nsig.giveNumberOfColumns();
            rhs.resize(regionDofMans, regionValSize);
            rhs.zero();
            if ( regionValSize == 0 ) {
                OOFEM_LOG_RELEVANT( "ZZNodalRecoveryModel :: unknown size of InternalStateType %s\n", __InternalStateTypeToString(type) );
            }
        } else if ( regionValSize != nsig.giveNumberOfColumns() ) {
            nsig.resize(regionDofMans, regionValSize);
            nsig.zero();
            OOFEM_LOG_RELEVANT( "ZZNodalRecoveryModel :: changing size of for InternalStateType %s. New sized results ignored (this shouldn't happen).\n", __InternalStateTypeToString(type) );
        }

        //loc.resize ((elemNodes+elemSides)*regionValSize);
        int eq = 1;
        for ( int elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
            int node = element->giveDofManager(elementNode)->giveNumber();
            lhs.at( regionNodalNumbers.at(node) ) += nn.at(eq);
            for ( int j = 1; j <= regionValSize; j++ ) {
                rhs.at(regionNodalNumbers.at(node), j) += nsig.at(eq, j);
            }

            eq++;
        }
    } // end assemble element contributions

#ifdef __PARALLEL_MODE
    if ( this->domain->giveEngngModel()->isParallel() ) {
        this->exchangeDofManValues(lhs, rhs, regionNodalNumbers);
    }
#endif

    sol.resize(regionDofMans * regionValSize);
    sol.zero();

    bool missingDofManContribution = false;
    unresolvedDofMans.clear();
    // solve for recovered values of active region
    for ( int i = 1; i <= regionDofMans; i++ ) {
        int eq = ( i - 1 ) * regionValSize;
        for ( int j = 1; j <= regionValSize; j++ ) {
            // rhs will be overriden by recovered values
            if ( fabs( lhs.at(i) ) > ZZNRM_ZERO_VALUE ) {
                sol.at(eq + j) = rhs.at(i, j) / lhs.at(i);
            } else {
                missingDofManContribution = true;
                unresolvedDofMans.insert( regionNodalNumbers.at(i) );
                sol.at(eq + j) = 0.0;
            }
        }
    }

    // update recovered values
    this->updateRegionRecoveredValues(regionNodalNumbers, regionValSize, sol);

    if ( missingDofManContribution ) {
        std :: ostringstream msg;
        int i = 0;
        for ( int dman: unresolvedDofMans ) {
            msg << this->domain->giveDofManager(dman)->giveLabel() << ' ';
            if ( ++i > 20 ) {
                break;
            }
        }
        if ( i > 20 ) {
            msg << "...";
        }
        OOFEM_WARNING("some values of some dofmanagers undetermined (in global numbers) \n[%s]", msg.str().c_str() );
    }


    this->valType = type;
    this->stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}


bool
ZZNodalRecoveryModelInterface :: ZZNodalRecoveryMI_computeNValProduct(FloatMatrix &answer, InternalStateType type,
                                                                      TimeStep *tStep)
{  // evaluates N^T sigma over element volume
   // N(nsigma, nsigma*nnodes)
   // Definition : sigmaVector = N * nodalSigmaVector
    FloatArray stressVector, n;
    FEInterpolation *interpol = element->giveInterpolation();
    IntegrationRule *iRule = element->giveDefaultIntegrationRulePtr();
    bool success = true;

    answer.clear();
    for ( GaussPoint *gp: *iRule ) {
        double dV = element->computeVolumeAround(gp);
        //this-> computeStressVector(stressVector, gp, tStep);
        if ( !element->giveIPValue(stressVector, gp, type, tStep) ) {
          success = false;
          continue;
        }

        interpol->evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(element) );
        answer.plusDyadUnsym(n, stressVector, dV);

        //  help.beTProductOf(n,stressVector);
        //  answer.add(help.times(dV));
    }
    return success;
}

void
ZZNodalRecoveryModelInterface :: ZZNodalRecoveryMI_computeNNMatrix(FloatArray &answer, InternalStateType type)
{
    //
    // Returns NTN matrix (lumped) for Zienkiewicz-Zhu
    // The size of N mtrx is (nstresses, nnodes*nstreses)
    // Definition : sigmaVector = N * nodalSigmaVector
    //
    double volume = 0.0;
    FloatMatrix fullAnswer;
    FloatArray n;
    FEInterpolation *interpol = element->giveInterpolation();
    IntegrationRule *iRule = element->giveDefaultIntegrationRulePtr();

    if ( !interpol ) {
        OOFEM_ERROR( "Element %d not providing interpolation", element->giveNumber() );
    }

    int size = element->giveNumberOfDofManagers();
    fullAnswer.resize(size, size);
    fullAnswer.zero();
    double pok = 0.0;

    for ( GaussPoint *gp: *iRule ) {
        double dV = element->computeVolumeAround(gp);
        interpol->evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(element) );
        fullAnswer.plusDyadSymmUpper(n, dV);
        pok += ( n.at(1) * dV ); ///@todo What is this? Completely unused.
        volume += dV;
    }


    fullAnswer.symmetrized();
    answer.resize(size);
    for ( int i = 1; i <= size; i++ ) {
        double sum = 0.0;
        for ( int j = 1; j <= size; j++ ) {
            sum += fullAnswer.at(i, j);
        }

        answer.at(i) = sum;
    }
}


#ifdef __PARALLEL_MODE

void
ZZNodalRecoveryModel :: initCommMaps()
{
 #ifdef __PARALLEL_MODE
    if ( initCommMap ) {
        EngngModel *emodel = domain->giveEngngModel();
        commBuff = new CommunicatorBuff(emodel->giveNumberOfProcesses(), CBT_dynamic);
        communicator = new NodeCommunicator(emodel, commBuff, emodel->giveRank(),
                                            emodel->giveNumberOfProcesses());
        communicator->setUpCommunicationMaps(domain->giveEngngModel(), true, true);
        OOFEM_LOG_INFO("ZZNodalRecoveryModel :: initCommMaps: initialized comm maps");
        initCommMap = false;
    }

 #endif
}

void
ZZNodalRecoveryModel :: exchangeDofManValues(FloatArray &lhs, FloatMatrix &rhs, IntArray &rn)
{
    parallelStruct ls( &lhs, &rhs, &rn);

    // exchange data for shared nodes
    communicator->packAllData(this, & ls, & ZZNodalRecoveryModel :: packSharedDofManData);
    communicator->initExchange(788);
    communicator->unpackAllData(this, & ls, & ZZNodalRecoveryModel :: unpackSharedDofManData);
    communicator->finishExchange();
}

int
ZZNodalRecoveryModel :: packSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm)
{
    int result = 1, indx, nc, size;
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    IntArray const *toSendMap = processComm.giveToSendMap();
    nc = s->rhs->giveNumberOfColumns();

    size = toSendMap->giveSize();
    for ( int i = 1; i <= size; i++ ) {
        // toSendMap contains all shared dofmans with remote partition
        // one has to check, if particular shared node value is available for given region
        indx = s->regionNodalNumbers->at( toSendMap->at(i) );
        if ( indx ) {
            // pack "1" to indicate that for given shared node this is a valid contribution
            result &= pcbuff->write(1);
            result &= pcbuff->write( s->lhs->at(indx) );
            for ( int j = 1; j <= nc; j++ ) {
                result &= pcbuff->write( s->rhs->at(indx, j) );
            }

            //printf("[%d] ZZ: Sending data for shred node %d[%d]\n", domain->giveEngngModel()->giveRank(),
            //       toSendMap->at(i), domain->giveDofManager(toSendMap->at(i))->giveGlobalNumber());
        } else {
            // ok shared node is not in active region (determined by s->regionNodalNumbers)
            result &= pcbuff->write(0);
        }
    }

    return result;
}

int
ZZNodalRecoveryModel :: unpackSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm)
{
    int result = 1;
    int nc, indx, size, flag;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    double value;
    nc = s->rhs->giveNumberOfColumns();

    size = toRecvMap->giveSize();
    for ( int i = 1; i <= size; i++ ) {
        indx = s->regionNodalNumbers->at( toRecvMap->at(i) );
        // toRecvMap contains all shared dofmans with remote partition
        // one has to check, if particular shared node received contribution is available for given region
        result &= pcbuff->read(flag);
        if ( flag ) {
            // "1" to indicates that for given shared node this is a valid contribution
            result &= pcbuff->read(value);
            // now check if we have a valid number
            if ( indx ) {
                s->lhs->at(indx) += value;
            }

            for ( int j = 1; j <= nc; j++ ) {
                result &= pcbuff->read(value);
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
