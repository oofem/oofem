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

#include "sprnodalrecoverymodel.h"
#include "timestep.h"
#include "element.h"
#include "node.h"
#include "connectivitytable.h"
#include "integrationrule.h"
#include "gausspoint.h"
#include "engngm.h"
#include "classfactory.h"

#ifdef __PARALLEL_MODE
 #include "processcomm.h"
 #include "problemcomm.h"
#endif

#include <cstdlib>
#include <list>

namespace oofem {
REGISTER_NodalRecoveryModel(SPRNodalRecoveryModel, NodalRecoveryModel :: NRM_SPR);

SPRNodalRecoveryModel :: SPRNodalRecoveryModel(Domain *d) : NodalRecoveryModel(d)
{ }

SPRNodalRecoveryModel :: ~SPRNodalRecoveryModel()
{ }

int
SPRNodalRecoveryModel :: recoverValues(Set elementSet, InternalStateType type, TimeStep *tStep)
{
    int nnodes = domain->giveNumberOfDofManagers();
    IntArray regionNodalNumbers(nnodes);
    IntArray patchElems, dofManToDetermine, pap;
    FloatMatrix a;
    FloatArray dofManValues;
    IntArray dofManPatchCount;

    if ( ( this->valType == type ) && ( this->stateCounter == tStep->giveSolutionStateCounter() ) ) {
        return 1;
    }

#ifdef __PARALLEL_MODE
    this->initCommMaps();
#endif

    // clear nodal table
    this->clear();

    int regionValSize;
    int regionDofMans;


    regionValSize = 0;
    // loop over elements and determine local region node numbering and determine and check nodal values size
    if ( this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, elementSet) == 0 ) {
        return 0;
    }

    SPRPatchType regType = this->determinePatchType(elementSet);

    dofManPatchCount.resize(regionDofMans);
    dofManPatchCount.zero();


    //pap = patch assembly points
    this->determinePatchAssemblyPoints(pap, regType, elementSet);

    int npap = pap.giveSize();
    for ( int ipap = 1; ipap <= npap; ipap++ ) {
        int papNumber = pap.at(ipap);
        int oldSize = regionValSize;

        this->initPatch(patchElems, dofManToDetermine, pap, papNumber, elementSet);
        this->computePatch(a, patchElems, regionValSize, regType, type, tStep);
        if ( oldSize == 0 ) {
            dofManValues.resize(regionDofMans * regionValSize);
            dofManValues.zero();
        }
        this->determineValuesFromPatch(dofManValues, dofManPatchCount, regionNodalNumbers,
                                       dofManToDetermine, a, regType);
    }

#ifdef __PARALLEL_MODE
    this->exchangeDofManValues(dofManValues, dofManPatchCount, regionNodalNumbers, regionValSize);
#endif

    // average  recovered values of active region
    bool abortFlag = false;
    for ( int i = 1; i <= nnodes; i++ ) {
        if ( regionNodalNumbers.at(i) &&
            ( ( domain->giveDofManager(i)->giveParallelMode() == DofManager_local ) ||
             ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) ) ) {
            int eq = ( regionNodalNumbers.at(i) - 1 ) * regionValSize;
            if ( dofManPatchCount.at( regionNodalNumbers.at(i) ) ) {
                for ( int j = 1; j <= regionValSize; j++ ) {
                    dofManValues.at(eq + j) /= dofManPatchCount.at( regionNodalNumbers.at(i) );
                }
            } else {
                OOFEM_WARNING("values of %s in dofmanager %d undetermined", __InternalStateTypeToString(type), i);

                for ( int j = 1; j <= regionValSize; j++ ) {
                    dofManValues.at(eq + j) = 0.0;
                }
                //abortFlag = true;
            }
        }

        if ( abortFlag ) {
            abort();
        }

        // update recovered values
        this->updateRegionRecoveredValues(regionNodalNumbers, regionValSize, dofManValues);
    }

    this->valType = type;
    this->stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}

void
SPRNodalRecoveryModel :: determinePatchAssemblyPoints(IntArray &pap, SPRPatchType regType, Set &elementSet)
{
    int idofMan, ndofMan = domain->giveNumberOfDofManagers();
    int ielem;
    int npap, ipap, count, neq, nip;
    IntArray dofManFlags(ndofMan);
    IntArray elemPap;
    SPRNodalRecoveryModelInterface *interface;
    Element *element;
    const IntArray *papDofManConnectivity;
    enum { papStatus_noPap, papStatus_regular, papStatus_boundary, papStatus_littleNIPs };

    // init all dof man statuses
    for ( idofMan = 1; idofMan <= ndofMan; idofMan++ ) {
        dofManFlags.at(idofMan) = papStatus_noPap;
    }

    IntArray elements = elementSet.giveElementList();
    // assign all possible paps with corresponding count
    for ( int i = 1; i <= elements.giveSize(); i++ ) {
        ielem = elements.at(i);
        element = domain->giveElement(ielem);

        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

        if ( ( interface = static_cast< SPRNodalRecoveryModelInterface * >( element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) ) {
            interface->SPRNodalRecoveryMI_giveSPRAssemblyPoints(elemPap);
            npap = elemPap.giveSize();
            for ( ipap = 1; ipap <= npap; ipap++ ) {
                dofManFlags.at( elemPap.at(ipap) ) = papStatus_regular;
            }
        }
    }

    // after loop all possible paps (patch assembly points) will have papStatus_regular flag

    // but we now have to skip those pap reported by elements, which have not enough integration points
    // to determine the least square fit of patch
    // and also we mark those dofManagers which are on boundary

    neq = this->giveNumberOfUnknownPolynomialCoefficients(regType);
    for ( idofMan = 1; idofMan <= ndofMan; idofMan++ ) {
        // mark boundary dofManagers
        if ( domain->giveDofManager(idofMan)->isBoundary() ) {
            dofManFlags.at(idofMan) = papStatus_boundary;
        }

        nip = 0;
        if ( dofManFlags.at(idofMan) != papStatus_noPap ) {
            papDofManConnectivity = domain->giveConnectivityTable()->giveDofManConnectivityArray(idofMan);
            for ( ielem = 1; ielem <= papDofManConnectivity->giveSize(); ielem++ ) {
                element = domain->giveElement( papDofManConnectivity->at(ielem) );

                if ( element->giveParallelMode() != Element_local ) {
                    continue;
                }

                if ( elementSet.hasElement( element->giveNumber() ) ) {
                    if ( ( interface = static_cast< SPRNodalRecoveryModelInterface * >( element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) ) {
                        nip += interface->SPRNodalRecoveryMI_giveNumberOfIP();
                    }
                }
            }

            if ( nip < neq ) {
                // this pap has not enough integration points to determine patch polynomial
                // reset its count to zero
                dofManFlags.at(idofMan) =  papStatus_littleNIPs;
            }
        }
    }

    //
    // generally boundary pap can be removed from pap list
    // if their value can be determined from other paps
    // and if they are  not the last resort to determine other dofManagers values (for example those with little nips).
    //
    //
    // here only test if paps with papStatus_littleNIPs can be determined using regular paps (papStatus_regular)
    // or the boundary paps must be employed. In such case these boundary paps are marked as regular ones
    // to force paches to be assembled.
    //
    bool foundRegularPap, foundBoundaryPap, abort_flag = false;
    // loop over boundary candidates to remove and try to confirm whether they can be removed
    for ( idofMan = 1; idofMan <= ndofMan; idofMan++ ) {
        foundRegularPap = foundBoundaryPap = false;
        if ( dofManFlags.at(idofMan) == papStatus_littleNIPs ) {
            papDofManConnectivity = domain->giveConnectivityTable()->giveDofManConnectivityArray(idofMan);
            for ( ielem = 1; ielem <= papDofManConnectivity->giveSize(); ielem++ ) {
                // try to determine if they can be determined from surronuding elements paps
                element = domain->giveElement( papDofManConnectivity->at(ielem) );

                if ( element->giveParallelMode() != Element_local ) {
                    continue;
                }

                if ( !elementSet.hasElement( element->giveNumber() ) ) {
                    continue;
                }

                if ( ( interface = static_cast< SPRNodalRecoveryModelInterface * >( element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) ) {
                    interface->SPRNodalRecoveryMI_giveSPRAssemblyPoints(elemPap);
                    npap = elemPap.giveSize();
                    for ( ipap = 1; ipap <= npap; ipap++ ) {
                        // skip other dofMans with littleNIPs
                        if ( dofManFlags.at( elemPap.at(ipap) ) == papStatus_littleNIPs ) {
                            continue;
                        }

                        if ( dofManFlags.at( elemPap.at(ipap) ) == papStatus_regular ) {
                            foundRegularPap = true;
                        } else if ( dofManFlags.at( elemPap.at(ipap) ) == papStatus_boundary ) {
                            foundBoundaryPap = true;
                        }
                    }
                }
            }

            if ( foundRegularPap ) {
                continue;         // can be determined from regular pap - ok
            }

            // boundary dof man can be removed <= can be determined
            if ( foundBoundaryPap ) {
                // try the last possibility-> determine its value from boundary patches
                // mark boundaryPap as regulars -> they can be used to assemble patch (they have enough nips)
                for ( ielem = 1; ielem <= papDofManConnectivity->giveSize(); ielem++ ) {
                    element = domain->giveElement( papDofManConnectivity->at(ielem) );
                    if ( !elementSet.hasElement( element->giveNumber() ) ) {
                        continue;
                    }

                    if ( ( interface = static_cast< SPRNodalRecoveryModelInterface * >( element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) ) {
                        interface->SPRNodalRecoveryMI_giveSPRAssemblyPoints(elemPap);
                        npap = elemPap.giveSize();
                        for ( ipap = 1; ipap <= npap; ipap++ ) {
                            if ( dofManFlags.at( elemPap.at(ipap) ) == papStatus_boundary ) {
                                // change status to regular pap
                                dofManFlags.at( elemPap.at(ipap) ) = papStatus_regular;
                            }
                        }
                    }
                }
            } else {
                // if the pap with papStatus_littleNIPs status found, which values could not be determined using
                // regular pap or boundary pap then we are unable to determine such value
                if ( dofManFlags.at(idofMan) == papStatus_littleNIPs ) {
                    OOFEM_WARNING("unable to determine dofMan %d\n", idofMan);
                    //abort_flag = true;
                }
            }
        }
    }

    if ( abort_flag ) {
        abort();
    }


    count = 0;
    // count regular paps - those for which patch will be assembled
    for ( idofMan = 1; idofMan <= ndofMan; idofMan++ ) {
        if ( dofManFlags.at(idofMan) ==  papStatus_regular ) {
            count++;
        }
    }

    pap.resize(count);
    count = 0;
    for ( idofMan = 1; idofMan <= ndofMan; idofMan++ ) {
        if ( dofManFlags.at(idofMan) ==  papStatus_regular ) {
            pap.at(++count) = idofMan;
        }
    }
}


void
SPRNodalRecoveryModel :: initPatch(IntArray &patchElems, IntArray &dofManToDetermine,
                                   IntArray &pap, int papNumber, Set &elementSet)
{
    int nelem, ndofman, count, patchElements, j, includes, npap, ipap;
    const IntArray *papDofManConnectivity = domain->giveConnectivityTable()->giveDofManConnectivityArray(papNumber);
    std :: list< int >dofManToDetermineList;
    SPRNodalRecoveryModelInterface *interface;
    IntArray toDetermine, toDetermine2, elemPap, papInv;
    Element *element;

    IntArray regionelements = elementSet.giveElementList();

    // loop over elements sharing dofManager with papNumber and
    // determine those in region in ireg
    //
    nelem = papDofManConnectivity->giveSize();
    count = 0;
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        if ( domain->giveElement( papDofManConnectivity->at(ielem) )->giveParallelMode() != Element_local ) {
            continue;
        }

        if ( elementSet.hasElement(papDofManConnectivity->at(ielem)) ) {
            count++;
        }
    }

    patchElems.resize(count);
    patchElements = 0;
    //for ( int i = 1; i <= nelem; i++ ) {
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        //int ielem = regionelements.at(i);

        if ( domain->giveElement( papDofManConnectivity->at(ielem) )->giveParallelMode() != Element_local ) {
            continue;
        }

        if ( elementSet.hasElement(papDofManConnectivity->at(ielem)) ) {
            patchElems.at(++patchElements) = papDofManConnectivity->at(ielem);
        }
    }

    // Invert the pap array for faster access later
    ndofman = this->domain->giveNumberOfDofManagers();
    papInv.resize(ndofman);
    papInv.zero();
    for ( int i = 1; i <= pap.giveSize(); ++i ) {
        papInv.at( pap.at(i) ) = 1;
    }

    // determine dofManagers which values will be determined by this patch
    // first add those required by elements participating in patch
    dofManToDetermine.clear();
    for (int ielem = 1; ielem <= patchElements; ielem++ ) {
        element = domain->giveElement( patchElems.at(ielem) );

        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

        if ( ( interface = static_cast< SPRNodalRecoveryModelInterface * >( element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) ) {
            // add element reported dofMans for pap dofMan
            interface->SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(toDetermine, papNumber);
            for ( int i = 1; i <= toDetermine.giveSize(); i++ ) {
                includes = 0;
                // test if INCLUDED
                for ( int dman: dofManToDetermineList ) {
                    if ( dman == toDetermine.at(i) ) {
                        includes = 1;
                    }
                }

                if ( !includes ) {
                    dofManToDetermineList.push_back( toDetermine.at(i) );
                }

                // determine those dofManagers which are not reported by elements,
                // but their values shoud be determined from this patch
                // Example include pap DofMans with connectivity 1
                interface->SPRNodalRecoveryMI_giveSPRAssemblyPoints(elemPap);
                npap = elemPap.giveSize();
                for ( ipap = 1; ipap <= npap; ipap++ ) {
                    // test if element reported SPRAssembly point is not global assembly point
                    // then determine this point from this patch
                    if ( papInv.at( elemPap.at(ipap) ) ==  0 ) {
                        includes = 0;
                        // test if INCLUDED
                        for ( int dman: dofManToDetermineList ) {
                            if ( dman == elemPap.at(ipap) ) {
                                includes = 1;
                            }
                        }

                        if ( !includes ) {
                            dofManToDetermineList.push_back( elemPap.at(ipap) );
                        }

                        // add also all dofManagers which are reported by element for this Assembly node
                        interface->SPRNodalRecoveryMI_giveDofMansDeterminedByPatch( toDetermine2, elemPap.at(ipap) );
                        for ( j = 1; j <= toDetermine2.giveSize(); j++ ) {
                            includes = 0;
                            // test if INCLUDED
                            for ( int dman: dofManToDetermineList ) {
                                if ( dman == toDetermine2.at(j) ) {
                                    includes = 1;
                                }
                            }

                            if ( !includes ) {
                                dofManToDetermineList.push_back( toDetermine2.at(j) );
                            }
                        }
                    }
                }
            }
        }
    } // end loop over patch elements

    // transform set to dofManToDetermine array
    count = (int)dofManToDetermineList.size();

    dofManToDetermine.resize(count);

    count = 0;
    for ( int dman: dofManToDetermineList ) {
        dofManToDetermine.at(++count) = dman;
    }
}



void
SPRNodalRecoveryModel :: computePatch(FloatMatrix &a, IntArray &patchElems, int &regionValSize,
                                      SPRPatchType regType, InternalStateType type, TimeStep *tStep)
{
    int nelem, neq;
    FloatArray ipVal, coords, P;
    FloatMatrix A, rhs;

    neq = this->giveNumberOfUnknownPolynomialCoefficients(regType);
    rhs.resize(neq, regionValSize);
    rhs.zero();
    A.resize(neq, neq);
    A.zero();

    // loop over elements in patch
    nelem = patchElems.giveSize();
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *element = domain->giveElement( patchElems.at(ielem) );
        if ( element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) {
            IntegrationRule *iRule = element->giveDefaultIntegrationRulePtr();
            for ( GaussPoint *gp: *iRule ) {
                int hasVal = element->giveIPValue(ipVal, gp, type, tStep);
                if ( !hasVal ) {
                    ipVal.resize(regionValSize);
                    ipVal.zero();
                } else if ( regionValSize == 0 ) {
                    regionValSize = ipVal.giveSize();
                    rhs.resize(neq, regionValSize);
                    rhs.zero();
                }

                element->computeGlobalCoordinates( coords, gp->giveSubPatchCoordinates() );
                // compute ip contribution
                this->computePolynomialTerms(P, coords, regType);
                for ( int j = 1; j <= neq; j++ ) {
                    for ( int k = 1; k <= regionValSize; k++ ) {
                        rhs.at(j, k) += P.at(j) * ipVal.at(k);
                    }

                    for ( int k = 1; k <= neq; k++ ) {
                        A.at(j, k) += P.at(j) * P.at(k);
                    }
                }
            } // end loop over nip
        }
    } // end loop over elements

    A.solveForRhs(rhs, a);
}

void
SPRNodalRecoveryModel :: determineValuesFromPatch(FloatArray &dofManValues, IntArray &dofManCount,
                                                  IntArray &regionNodalNumbers, IntArray &dofManToDetermine,
                                                  FloatMatrix &a, SPRPatchType type)
{
    int ndofMan = dofManToDetermine.giveSize();
    FloatArray P, vals;

    for ( int dofMan = 1; dofMan <= ndofMan; dofMan++ ) {
        const auto &coords = domain->giveNode( dofManToDetermine.at(dofMan) )->giveCoordinates();
        this->computePolynomialTerms(P, coords, type);
        vals.beTProductOf(a, P);

        // assemble values
        int eq = ( regionNodalNumbers.at( dofManToDetermine.at(dofMan) ) - 1 ) * vals.giveSize();
        for ( int i = 1; i <= vals.giveSize(); i++ ) {
            dofManValues.at(eq + i) += vals.at(i);
        }

        dofManCount.at( regionNodalNumbers.at( dofManToDetermine.at(dofMan) ) )++;
    }
}

void
SPRNodalRecoveryModel :: computePolynomialTerms(FloatArray &P, const FloatArray &coords, SPRPatchType type)
{
    if ( type == SPRPatchType_2dxy ) {
        P.resize(3);
        P.at(1) = 1.0;
        P.at(2) = coords.at(1);
        P.at(3) = coords.at(2);
    } else if ( type == SPRPatchType_3dBiLin ) {
        P.resize(4);
        P.at(1) = 1.0;
        P.at(2) = coords.at(1);
        P.at(3) = coords.at(2);
        P.at(4) = coords.at(3);
    } else if ( type == SPRPatchType_2dquadratic ) {
        P.resize(6);
        P.at(1) = 1.0;
        P.at(2) = coords.at(1);
        P.at(3) = coords.at(2);
        P.at(4) = coords.at(1) * coords.at(2);
        P.at(5) = coords.at(1) * coords.at(1);
        P.at(6) = coords.at(2) * coords.at(2);
    } else if ( type == SPRPatchType_3dBiQuadratic ) {
        P.resize(10);
        P.at(1) = 1.0;
        P.at(2) = coords.at(1);
        P.at(3) = coords.at(2);
        P.at(4) = coords.at(3);
        P.at(5) = coords.at(1) * coords.at(1);
        P.at(6) = coords.at(1) * coords.at(2);
        P.at(7) = coords.at(1) * coords.at(3);
        P.at(8) = coords.at(2) * coords.at(2);
        P.at(9) = coords.at(2) * coords.at(3);
        P.at(10) = coords.at(3) * coords.at(3);
    } else {
        OOFEM_ERROR("unknown regionType");
    }
}

int
SPRNodalRecoveryModel :: giveNumberOfUnknownPolynomialCoefficients(SPRPatchType regType)
{
    if ( regType == SPRPatchType_2dxy ) {
        return 3;
    } else if ( regType == SPRPatchType_3dBiLin ) {
        return 4;
    } else if ( regType == SPRPatchType_2dquadratic ) {
        return 6;
    } else if ( regType == SPRPatchType_3dBiQuadratic ) {
        return 10;
    } else {
        return 0;
    }
}


SPRPatchType SPRNodalRecoveryModel :: determinePatchType(Set &elementList)
{
    SPRNodalRecoveryModelInterface *interface;
    if ( elementList.giveElementList().giveSize() ) {
        Element *e = domain->giveElement( elementList.giveElementList().at(1) );
        if ( ( interface = static_cast< SPRNodalRecoveryModelInterface * >( e->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) ) {
            return interface->SPRNodalRecoveryMI_givePatchType();
        } else {
            OOFEM_ERROR("unable to determine region patchtype");
        }
    } else {
        OOFEM_ERROR("empty region set");
    }
    return SPRPatchType_none; // to make compiler happy
}



#ifdef __PARALLEL_MODE

void
SPRNodalRecoveryModel :: initCommMaps()
{
 #ifdef __PARALLEL_MODE
    if ( initCommMap ) {
        EngngModel *emodel = domain->giveEngngModel();
        commBuff = new CommunicatorBuff(emodel->giveNumberOfProcesses(), CBT_dynamic);
        communicator = new NodeCommunicator(emodel, commBuff, emodel->giveRank(),
                                            emodel->giveNumberOfProcesses());
        communicator->setUpCommunicationMaps(domain->giveEngngModel(), true, true);
        OOFEM_LOG_INFO("SPRNodalRecoveryModel :: initCommMaps: initialized comm maps");
        initCommMap = false;
    }

 #endif
}

void
SPRNodalRecoveryModel :: exchangeDofManValues(FloatArray &dofManValues, IntArray &dofManPatchCount,
                                              IntArray &regionNodalNumbers, int regionValSize)
{
    parallelStruct ls( &dofManValues, &dofManPatchCount, &regionNodalNumbers, regionValSize);

    // exchange data for shared nodes
    communicator->packAllData(this, & ls, & SPRNodalRecoveryModel :: packSharedDofManData);
    communicator->initExchange(789);
    communicator->unpackAllData(this, & ls, & SPRNodalRecoveryModel :: unpackSharedDofManData);
    communicator->finishExchange();
}

int
SPRNodalRecoveryModel :: packSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm)
{
    int result = 1;
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    const IntArray &toSendMap = processComm.giveToSendMap();

    for ( int inode : toSendMap ) {
        // toSendMap contains all shared dofmans with remote partition
        // one has to check, if particular shared node value is available for given region
        int indx = s->regionNodalNumbers->at( inode );
        if ( indx && s->dofManPatchCount->at(indx) ) {
            // pack "1" to indicate that for given shared node this is a valid contribution
            result &= pcbuff->write(1);
            int eq = ( indx - 1 ) * s->regionValSize;
            for ( int j = 1; j <= s->regionValSize; j++ ) {
                result &= pcbuff->write( s->dofManValues->at(eq + j) );
            }
        } else {
            // ok shared node is not in active region (determined by s->regionNodalNumbers)
            result &= pcbuff->write(0);
        }
    }

    return result;
}

int
SPRNodalRecoveryModel :: unpackSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm)
{
    int result = 1;
    int flag;
    const IntArray &toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();

    for ( int inode : toRecvMap ) {
        int indx = s->regionNodalNumbers->at( inode );
        // toRecvMap contains all shared dofmans with remote partition
        // one has to check, if particular shared node received contribution is available for given region
        result &= pcbuff->read(flag);
        if ( flag ) {
            // "1" to indicates that for given shared node this is a valid contribution
            int eq = ( indx - 1 ) * s->regionValSize;
            for ( int j = 1; j <= s->regionValSize; j++ ) {
                double value;
                result &= pcbuff->read(value);
                if ( indx ) {
                    s->dofManValues->at(eq + j) += value;
                }
            }

            if ( indx ) {
                s->dofManPatchCount->at(indx)++;
            }
        }
    }

    return result;
}

#endif
} // end namespace oofem
