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
ZZNodalRecoveryModel :: recoverValues(InternalStateType type, TimeStep *tStep)
{
    int nregions = this->giveNumberOfVirtualRegions();
    int nelem = domain->giveNumberOfElements();
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

    // loop over regions
    for ( int ireg = 1; ireg <= nregions; ireg++ ) {
        int elemNodes;
        int regionValSize;
        int regionDofMans;

        // loop over elements and determine local region node numbering and determine and check nodal values size
        if ( this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, ireg) == 0 ) {
            break;
        }

        regionValSize = 0;
        lhs.resize(regionDofMans);
        lhs.zero();
        // assemble element contributions
        for ( int ielem = 1; ielem <= nelem; ielem++ ) {
            ZZNodalRecoveryModelInterface *interface;
            Element *element = domain->giveElement(ielem);

#ifdef __PARALLEL_MODE
            if ( element->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            if ( this->giveElementVirtualRegionNumber(ielem) != ireg ) {
                continue;
            }

            // If an element doesn't implement the interface, it is ignored.
            if ( ( interface = static_cast< ZZNodalRecoveryModelInterface * >( element->giveInterface(ZZNodalRecoveryModelInterfaceType) ) ) == NULL ) {
                //abort();
                continue;
            }


            // ask element contributions
            interface->ZZNodalRecoveryMI_computeNNMatrix(nn, type);
            interface->ZZNodalRecoveryMI_computeNValProduct(nsig, type, tStep);
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
                for ( int i = 1; i <= regionValSize; i++ ) {
                    rhs.at(regionNodalNumbers.at(node), i) += nsig.at(eq, i);
                }

                eq++;
            }
        } // end assemble element contributions

#ifdef __PARALLEL_MODE
        if ( this->domain->giveEngngModel()->isParallel() ) {
            this->exchangeDofManValues(ireg, lhs, rhs, regionNodalNumbers);
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
        this->updateRegionRecoveredValues(ireg, regionNodalNumbers, regionValSize, sol);

        if ( missingDofManContribution ) {
            std :: ostringstream msg;
            int i = 0;
            for ( std :: set< int > :: const_iterator sit = unresolvedDofMans.begin(); sit != unresolvedDofMans.end(); ++sit ) {
                msg << * sit << ' ';
                if ( ++i > 20 ) {
                    break;
                }
            }
            if ( i > 20 ) {
                msg << "...";
            }
            OOFEM_WARNING3( "ZZNodalRecoveryModel::recoverValues: region %d: values of some dofmanagers undetermined\n[%s]", ireg, msg.str().c_str() );
        }
    } // end loop over regions

    this->valType = type;
    this->stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}


void
ZZNodalRecoveryModelInterface :: ZZNodalRecoveryMI_computeNValProduct(FloatMatrix &answer, InternalStateType type,
                                                                      TimeStep *tStep)
{  // evaluates N^T sigma over element volume
   // N(nsigma, nsigma*nnodes)
   // Definition : sigmaVector = N * nodalSigmaVector
    FloatArray stressVector, n;
    Element *elem  = this->ZZNodalRecoveryMI_giveElement();
    FEInterpolation *interpol = elem->giveInterpolation();
    IntegrationRule *iRule = elem->giveDefaultIntegrationRulePtr();

    answer.clear();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        double dV = elem->computeVolumeAround(gp);
        //this-> computeStressVector(stressVector, gp, tStep);
        if ( !elem->giveIPValue(stressVector, gp, type, tStep) ) {
            continue;
        }

        interpol->evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(elem) );
        answer.plusDyadUnsym(n, stressVector, dV);

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
    double volume = 0.0;
    FloatMatrix fullAnswer;
    FloatArray n;
    Element *elem  = this->ZZNodalRecoveryMI_giveElement();
    FEInterpolation *interpol = elem->giveInterpolation();
    IntegrationRule *iRule = elem->giveDefaultIntegrationRulePtr();

    if ( !interpol ) {
        OOFEM_ERROR2( "ZZNodalRecoveryMI_computeNNMatrix: Element %d not providing interpolation", elem->giveNumber() );
    }

    int size = elem->giveNumberOfDofManagers();
    fullAnswer.resize(size, size);
    fullAnswer.zero();
    double pok = 0.0;

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        double dV = elem->computeVolumeAround(gp);
        interpol->evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(elem) );
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
            result &= pcbuff->packInt(1);
            result &= pcbuff->packDouble( s->lhs->at(indx) );
            for ( int j = 1; j <= nc; j++ ) {
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
        result &= pcbuff->unpackInt(flag);
        if ( flag ) {
            // "1" to indicates that for given shared node this is a valid contribution
            result &= pcbuff->unpackDouble(value);
            // now check if we have a valid number
            if ( indx ) {
                s->lhs->at(indx) += value;
            }

            for ( int j = 1; j <= nc; j++ ) {
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
