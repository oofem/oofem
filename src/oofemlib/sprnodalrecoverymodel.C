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

#include "sprnodalrecoverymodel.h"
#include "timestep.h"
#include "element.h"
#include "node.h"
#include "conTable.h"
#include "integrationrule.h"
#include "gausspnt.h"

namespace oofem {
SPRNodalRecoveryModel :: SPRNodalRecoveryModel(Domain *d) : NodalRecoveryModel(d)
{ }

SPRNodalRecoveryModel :: ~SPRNodalRecoveryModel()
{ }

int
SPRNodalRecoveryModel :: recoverValues(InternalStateType type, TimeStep *tStep)
{
    int ireg, nregions = this->giveNumberOfVirtualRegions();
    int nnodes = domain->giveNumberOfDofManagers();
    int regionValSize;
    int regionDofMans;
    int i, j, neq, eq, npap, ipap, papNumber;
    IntArray skipRegionMap(nregions), regionRecSize(nregions);
    IntArray regionNodalNumbers(nnodes);
    IntArray patchElems, dofManToDetermine, pap, regionTypes;
    SPRPatchType regType;
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

    // init region table indicating regions to skip
    this->initRegionMap(skipRegionMap, regionRecSize, regionTypes, type);

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

        regType = ( SPRPatchType ) regionTypes.at(ireg);
        regionValSize = regionRecSize.at(ireg);
        // loop over elements and determine local region node numbering and determine and check nodal values size
        if ( this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, ireg) == 0 ) {
            continue;
        }

        neq = regionDofMans * regionValSize;
        dofManValues.resize(neq);
        dofManValues.zero();
        dofManPatchCount.resize(regionDofMans);
        dofManPatchCount.zero();

        //pap = patch assembly points
        this->determinePatchAssemblyPoints(pap, ireg, regType);

        npap = pap.giveSize();
        for ( ipap = 1; ipap <= npap; ipap++ ) {
            papNumber = pap.at(ipap);

            this->initPatch(patchElems, dofManToDetermine, pap, papNumber, ireg);
            this->computePatch(a, patchElems, papNumber, regionValSize, regType, type, tStep);
            this->determineValuesFromPatch(dofManValues, dofManPatchCount, regionNodalNumbers,
                                           dofManToDetermine, a, papNumber, regionValSize, regType);
        }

#ifdef __PARALLEL_MODE
        this->exchangeDofManValues(ireg, dofManValues, dofManPatchCount, regionNodalNumbers, regionValSize);
#endif

        // average  recovered values of active region
        bool abortFlag = false;
        for ( i = 1; i <= nnodes; i++ ) {
#ifndef __PARALLEL_MODE
            if ( regionNodalNumbers.at(i) ) {
#else
            if ( regionNodalNumbers.at(i) &&
                ( ( domain->giveDofManager(i)->giveParallelMode() == DofManager_local ) ||
                 ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) ) ) {
#endif
                eq = ( regionNodalNumbers.at(i) - 1 ) * regionValSize;
                if ( dofManPatchCount.at( regionNodalNumbers.at(i) ) ) {
                    for ( j = 1; j <= regionValSize; j++ ) {
                        dofManValues.at(eq + j) /= dofManPatchCount.at( regionNodalNumbers.at(i) );
                    }
                } else {
#ifndef __PARALLEL_MODE
                    OOFEM_WARNING3("SPRNodalRecoveryModel::recoverValues : values of %s in dofmanager %d undetermined", __InternalStateTypeToString(type), i);
#else
                    OOFEM_WARNING4("[%d] SPRNodalRecoveryModel::recoverValues : values of %s in dofmanager %d undetermined",
                                   domain->giveEngngModel()->giveRank(), __InternalStateTypeToString(type), i);
#endif
                    dofManValues.at(eq + j) = 0.0;
                    //abortFlag = true;
                }
            }
        }

        if ( abortFlag ) {
            abort();
        }

        // update recovered values
        this->updateRegionRecoveredValues(ireg, regionNodalNumbers, regionValSize, dofManValues);
    }

    this->valType = type;
    this->stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}

void
SPRNodalRecoveryModel :: initRegionMap(IntArray &regionMap, IntArray &regionValSize,
                                       IntArray &regionTypes, InternalStateType type)
{
    int nregions = this->giveNumberOfVirtualRegions();
    int ielem, nelem = domain->giveNumberOfElements();
    int i, regionsSkipped = 0;
    int elemVR;
    Element *element;
    SPRNodalRecoveryModelInterface *interface;

    regionMap.resize(nregions);
    regionMap.zero();
    regionValSize.resize(nregions);
    regionValSize.zero();
    regionTypes.resize(nregions);
    regionTypes.zero();

    // loop over elements and check if implement interface
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = domain->giveElement(ielem);
#ifdef __PARALLEL_MODE
        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

#endif
        if ( ( interface =  ( SPRNodalRecoveryModelInterface * ) element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) == NULL ) {
            // If an element doesn't implement the interface, it is ignored.
            //regionsSkipped = 1;
            //regionMap.at( element->giveRegionNumber() ) = 1;
            continue;
        } else {
            if ( ( elemVR = this->giveElementVirtualRegionNumber(ielem) ) ) { // test if elementVR is nonzero
                if ( regionValSize.at(elemVR) ) {
                    if ( regionValSize.at(elemVR) != interface->SPRNodalRecoveryMI_giveDofManRecordSize(type) ) {
                        // This indicates a size mis-match between different elements, no choice but to skip the region.
                        regionMap.at(elemVR) = 1;
                        regionsSkipped = 1;
                    }

                    if ( regionTypes.at(elemVR) != ( int ) interface->SPRNodalRecoveryMI_givePatchType() ) {
                        regionMap.at(elemVR) = 1;
                        /*
                         *   printf ("NodalRecoveryModel :: initRegionMap: element %d has incompatible Patch type, skipping region\n",ielem);
                         */
                        regionsSkipped = 1;
                    }
                } else {
                    regionValSize.at(elemVR) = interface->SPRNodalRecoveryMI_giveDofManRecordSize(type);
                    regionTypes.at(elemVR) = ( int ) interface->SPRNodalRecoveryMI_givePatchType();
                    if ( regionValSize.at(elemVR) == 0 ) {
                        regionMap.at(elemVR) = 1;
                        regionsSkipped = 1;
                    }
                }
            }
        }
    }

    if ( regionsSkipped ) {
        OOFEM_LOG_RELEVANT( "SPRNodalRecoveryModel :: initRegionMap: skipping regions for InternalStateType %s\n", __InternalStateTypeToString(type) );
        for ( i = 1; i <= nregions; i++ ) {
            if ( regionMap.at(i) ) {
                OOFEM_LOG_RELEVANT("%d ", i);
            }
        }

        OOFEM_LOG_RELEVANT("\n");
    }
}

void
SPRNodalRecoveryModel :: determinePatchAssemblyPoints(IntArray &pap, int ireg, SPRPatchType regType)
{
    int idofMan, ndofMan = domain->giveNumberOfDofManagers();
    int ielem, nelem = domain->giveNumberOfElements();
    int npap, ipap, count, neq, nip;
    IntArray dofManFlags(ndofMan), elemPap;
    SPRNodalRecoveryModelInterface *interface;
    Element *element;
    const IntArray *papDofManConnectivity;
    enum { papStatus_noPap, papStatus_regular, papStatus_boundary, papStatus_littleNIPs };

    // init all dof man statuses
    for ( idofMan = 1; idofMan <= ndofMan; idofMan++ ) {
        dofManFlags.at(idofMan) = papStatus_noPap;
    }

    // assign all possible paps with corresponding count
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

        if ( ( interface = ( SPRNodalRecoveryModelInterface * ) element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) {
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
#ifdef __PARALLEL_MODE
                if ( element->giveParallelMode() != Element_local ) {
                    continue;
                }

#endif
                if ( this->giveElementVirtualRegionNumber( element->giveNumber() ) == ireg ) {
                    if ( ( interface = ( SPRNodalRecoveryModelInterface * ) element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) {
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
#ifdef __PARALLEL_MODE
                if ( element->giveParallelMode() != Element_local ) {
                    continue;
                }

#endif
                if ( this->giveElementVirtualRegionNumber( element->giveNumber() ) != ireg ) {
                    continue;
                }

                if ( ( interface = ( SPRNodalRecoveryModelInterface * ) element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) {
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
                    if ( this->giveElementVirtualRegionNumber( element->giveNumber() ) != ireg ) {
                        continue;
                    }

                    if ( ( interface = ( SPRNodalRecoveryModelInterface * ) element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) {
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
                    OOFEM_WARNING2("SPRNodalRecoveryModel::determinePatchAssemblyPoints - unable to determine dofMan %d\n", idofMan);
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
                                   IntArray &pap, int papNumber, int ireg)
{
    int nelem, ndofman, ielem, count, patchElements, i, j, includes, npap, ipap;
    const IntArray *papDofManConnectivity = domain->giveConnectivityTable()->giveDofManConnectivityArray(papNumber);
    dynaList< int >dofManToDetermineList;
    dynaList< int > :: iterator dofManToDetermineListIter;
    SPRNodalRecoveryModelInterface *interface;
    IntArray toDetermine, toDetermine2, elemPap, papInv;
    Element *element;

    // looop over elements sharing dofManager with papNumber and
    // determine those in region in ireg
    //
    nelem = papDofManConnectivity->giveSize();
    count = 0;
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
#ifdef __PARALLEL_MODE
        if ( domain->giveElement( papDofManConnectivity->at(ielem) )->giveParallelMode() != Element_local ) {
            continue;
        }

#endif
        if ( this->giveElementVirtualRegionNumber( papDofManConnectivity->at(ielem) ) == ireg ) {
            count++;
        }
    }

    patchElems.resize(count);
    patchElements = 0;
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
#ifdef __PARALLEL_MODE
        if ( domain->giveElement( papDofManConnectivity->at(ielem) )->giveParallelMode() != Element_local ) {
            continue;
        }

#endif
        if ( this->giveElementVirtualRegionNumber( papDofManConnectivity->at(ielem) ) == ireg ) {
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
    dofManToDetermine.resize(0);
    for ( ielem = 1; ielem <= patchElements; ielem++ ) {
        element = domain->giveElement( patchElems.at(ielem) );

#ifdef __PARALLEL_MODE
        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

#endif
        if ( ( interface = ( SPRNodalRecoveryModelInterface * ) element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) {
            // add element reported dofMans for pap dofMan
            interface->SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(toDetermine, papNumber);
            for ( i = 1; i <= toDetermine.giveSize(); i++ ) {
                includes = 0;
                // test if INCLUDED
                for ( dofManToDetermineListIter = dofManToDetermineList.begin();
                      dofManToDetermineListIter != dofManToDetermineList.end();
                      ++dofManToDetermineListIter ) {
                    if ( * ( dofManToDetermineListIter ) == toDetermine.at(i) ) {
                        includes = 1;
                    }
                }

                if ( !includes ) {
                    dofManToDetermineList.pushBack( toDetermine.at(i) );
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
                        for ( dofManToDetermineListIter = dofManToDetermineList.begin();
                              dofManToDetermineListIter != dofManToDetermineList.end();
                              ++dofManToDetermineListIter ) {
                            if ( ( * dofManToDetermineListIter ) == elemPap.at(ipap) ) {
                                includes = 1;
                            }
                        }

                        if ( !includes ) {
                            dofManToDetermineList.pushBack( elemPap.at(ipap) );
                        }

                        // add also all dofManagers which are reported by element for this Assembly node
                        interface->SPRNodalRecoveryMI_giveDofMansDeterminedByPatch( toDetermine2, elemPap.at(ipap) );
                        for ( j = 1; j <= toDetermine2.giveSize(); j++ ) {
                            includes = 0;
                            // test if INCLUDED
                            for ( dofManToDetermineListIter = dofManToDetermineList.begin();
                                  dofManToDetermineListIter != dofManToDetermineList.end();
                                  ++dofManToDetermineListIter ) {
                                if ( * dofManToDetermineListIter == toDetermine2.at(j) ) {
                                    includes = 1;
                                }
                            }

                            if ( !includes ) {
                                dofManToDetermineList.pushBack( toDetermine2.at(j) );
                            }
                        }
                    }
                }
            }
        }
    } // end loop over patch elements

    // transform set to dofManToDetermine array
    count = 0;
    for ( dofManToDetermineListIter = dofManToDetermineList.begin();
          dofManToDetermineListIter != dofManToDetermineList.end();
          ++dofManToDetermineListIter ) {
        count++;
    }

    dofManToDetermine.resize(count);

    count = 0;
    for ( dofManToDetermineListIter = dofManToDetermineList.begin();
          dofManToDetermineListIter != dofManToDetermineList.end();
          ++dofManToDetermineListIter ) {
        dofManToDetermine.at(++count) = * dofManToDetermineListIter;
    }
}



void
SPRNodalRecoveryModel :: computePatch(FloatMatrix &a, IntArray &patchElems, int papNumber, int regionValSize,
                                      SPRPatchType regType, InternalStateType type, TimeStep *tStep)
{
    int i, j, k, nelem, ielem, nip, neq;
    Element *element;
    FloatArray ipVal, coords, P;
    FloatMatrix A, rhs;
    SPRNodalRecoveryModelInterface *interface;
    IntegrationRule *iRule;
    GaussPoint *gp;

    neq = this->giveNumberOfUnknownPolynomialCoefficients(regType);
    a.resize(neq, regionValSize);
    a.zero();
    rhs.resize(neq, regionValSize);
    rhs.zero();
    A.resize(neq, neq);
    A.zero();

    // loop over elements in patch
    nelem = patchElems.giveSize();
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = domain->giveElement( patchElems.at(ielem) );
        if ( ( interface = ( SPRNodalRecoveryModelInterface * ) element->giveInterface(SPRNodalRecoveryModelInterfaceType) ) ) {
            iRule = element->giveDefaultIntegrationRulePtr();
            nip = iRule->getNumberOfIntegrationPoints();
            for ( i = 0; i < nip; i++ ) {
                gp  = iRule->getIntegrationPoint(i);
                if ( !element->giveIPValue(ipVal, gp, type, tStep) ) {
                    ipVal.resize(regionValSize);
                    ipVal.zero();
                }

                element->computeGlobalCoordinates( coords, * gp->giveLocalCoordinates() );
                // compute ip contribution
                this->computePolynomialTerms(P, coords, regType);
                for ( j = 1; j <= neq; j++ ) {
                    for ( k = 1; k <= regionValSize; k++ ) {
                        rhs.at(j, k) += P.at(j) * ipVal.at(k);
                    }

                    for ( k = 1; k <= neq; k++ ) {
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
                                                  FloatMatrix &a, int papNumber, int regionValSize,
                                                  SPRPatchType type)
{
    int i, j, eq, dofMan, ndofMan = dofManToDetermine.giveSize();
    FloatArray P, *coords, vals(regionValSize);
    int lneq = this->giveNumberOfUnknownPolynomialCoefficients(type);

    for ( dofMan = 1; dofMan <= ndofMan; dofMan++ ) {
        vals.zero();
        coords = domain->giveNode( dofManToDetermine.at(dofMan) )->giveCoordinates();
        this->computePolynomialTerms(P, * coords, type);
        for ( i = 1; i <= regionValSize; i++ ) {
            for ( j = 1; j <= lneq; j++ ) {
                vals.at(i) += P.at(j) * a.at(j, i);
            }
        }

        // assemble values

        eq = ( regionNodalNumbers.at( dofManToDetermine.at(dofMan) ) - 1 ) * regionValSize;
        for ( i = 1; i <= regionValSize; i++ ) {
            dofManValues.at(eq + i) += vals.at(i);
        }

        dofManCount.at( regionNodalNumbers.at( dofManToDetermine.at(dofMan) ) )++;
    }
}

void
SPRNodalRecoveryModel :: computePolynomialTerms(FloatArray &P, FloatArray &coords, SPRPatchType type)
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
        OOFEM_ERROR("SPRNodalRecoveryModel::computePolynomialTerms - unknown regionType");
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



#ifdef __PARALLEL_MODE

void
SPRNodalRecoveryModel :: initCommMaps()
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
            OOFEM_LOG_INFO("SPRNodalRecoveryModel :: initCommMaps: initialized comm maps");
            initCommMap = false;
        } else {
            OOFEM_ERROR("SPRNodalRecoveryModel :: initCommMaps: unsupported comm mode");
        }
    }

 #endif
}

void
SPRNodalRecoveryModel :: exchangeDofManValues(int ireg, FloatArray &dofManValues, IntArray &dofManPatchCount,
                                              IntArray &regionNodalNumbers, int regionValSize)
{
    EngngModel *emodel = domain->giveEngngModel();
    ProblemCommunicatorMode commMode = emodel->giveProblemCommMode();

    if ( commMode == ProblemCommMode__NODE_CUT ) {
        parallelStruct ls( &dofManValues, &dofManPatchCount, &regionNodalNumbers, regionValSize);

        // exchange data for shared nodes
        communicator->packAllData(this, & ls, & SPRNodalRecoveryModel :: packSharedDofManData);
        communicator->initExchange(789 + ireg);
        communicator->unpackAllData(this, & ls, & SPRNodalRecoveryModel :: unpackSharedDofManData);
        communicator->finishExchange();
    } else {
        OOFEM_ERROR("SPRNodalRecoveryModel :: exchangeDofManValues: Unsupported commMode");
    }
}

int
SPRNodalRecoveryModel :: packSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm)
{
    int result = 1, i, j, indx, eq, size;
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    IntArray const *toSendMap = processComm.giveToSendMap();

    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        // toSendMap contains all shared dofmans with remote partition
        // one has to check, if particular shared node value is available for given region
        indx = s->regionNodalNumbers->at( toSendMap->at(i) );
        if ( indx && s->dofManPatchCount->at(indx) ) {
            // pack "1" to indicate that for given shared node this is a valid contribution
            result &= pcbuff->packInt(1);
            eq = ( indx - 1 ) * s->regionValSize;
            for ( j = 1; j <= s->regionValSize; j++ ) {
                result &= pcbuff->packDouble( s->dofManValues->at(eq + j) );
            }
        } else {
            // ok shared node is not in active region (determined by s->regionNodalNumbers)
            result &= pcbuff->packInt(0);
        }
    }

    return result;
}

int
SPRNodalRecoveryModel :: unpackSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm)
{
    int result = 1;
    int i, j, eq, indx, size, flag;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    double value;
    bool accept;

    size = toRecvMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        indx = s->regionNodalNumbers->at( toRecvMap->at(i) );
        accept = indx && s->dofManPatchCount->at(indx);
        // toRecvMap contains all shared dofmans with remote partition
        // one has to check, if particular shared node received contribution is available for given region
        result &= pcbuff->unpackInt(flag);
        if ( flag ) {
            // "1" to indicates that for given shared node this is a valid contribution
            eq = ( indx - 1 ) * s->regionValSize;
            for ( j = 1; j <= s->regionValSize; j++ ) {
                result &= pcbuff->unpackDouble(value);
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
