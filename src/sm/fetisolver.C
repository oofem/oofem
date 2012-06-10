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

#ifdef __PARALLEL_MODE

#include "mathfem.h"
#include "fetisolver.h"
#include "compiler.h"
#include "skyline.h"
#include "verbose.h"
#include "dofmanager.h"
#include "engngm.h"
#include "feticommunicator.h"

namespace oofem {
FETISolver :: FETISolver(int i, Domain *d, EngngModel *m) : SparseLinearSystemNM(i, d, m), pcbuff(CBT_static), processCommunicator(m, & pcbuff, 0)
{
    err    = 1.e-6;
    ni     = 20;
    energyNorm_comput_flag = 0;
}


FETISolver :: ~FETISolver()
{
    if ( engngModel->giveRank() == 0 ) {
        delete commBuff;
        delete masterCommunicator;
    }
}



int
FETISolver :: estimateMaxPackSize(IntArray &map, CommunicationBuffer &buff, int &packUnpackType)
{
    int rank = domain->giveEngngModel()->giveRank();
    int mapSize = map.giveSize();
    int i, j, eqNum, ndofs, count = 0;
    IntArray locationArray;
    EModelDefaultEquationNumbering dn;

    if ( rank == 0 ) {
        // master comm maps contain boundary dof managers
        for ( i = 1; i <= mapSize; i++ ) {
            count += masterCommunicator->giveDofManager( map.at(i) )->giveNumberOfDofs();
        }
    } else {
        for ( i = 1; i <= mapSize; i++ ) {
            domain->giveDofManager( map.at(i) )->giveCompleteLocationArray(locationArray, dn);
            ndofs = locationArray.giveSize();
            for ( j = 1; j <= ndofs; j++ ) {
                if ( ( eqNum = locationArray.at(j) ) ) {
                    count++;
                }
            }
        }
    }

    return ( buff.givePackSize(MPI_DOUBLE, 1) * count * FETISOLVER_MAX_RBM );
}


IRResultType
FETISolver :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                            // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, ni, IFT_FETISolver_maxiter, "maxiter"); // Macro
    IR_GIVE_FIELD(ir, err, IFT_FETISolver_maxerr, "maxerr"); // Macro
    IR_GIVE_FIELD(ir, limit, IFT_FETISolver_limit, "limit"); // Macro
    energyNorm_comput_flag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, energyNorm_comput_flag, IFT_FETISolver_energynormflag, "energynormflag"); // Macro

    if ( fabs(limit) < 1.e-20 ) {
        limit = 1.e-20;
    }

    if ( err < 1.e-20 ) {
        err = 1.e-20;
    }

    if ( err > 0.1 ) {
        err = 0.1;
    }

    return IRRT_OK;
}

void FETISolver :: setUpCommunicationMaps()
{
    int nnodes = domain->giveNumberOfDofManagers();
    int boundaryDofManNum = 0;
    int i, j, indx = 1, neq;
    StaticCommunicationBuffer commBuff(MPI_COMM_WORLD);
    IntArray commMap;
    EModelDefaultEquationNumbering dn;

    // determine the total number of boundary dof managers
    for ( i = 1; i <= nnodes; i++ ) {
        if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
            boundaryDofManNum++;
        }
    }

    if ( domain->giveEngngModel()->giveRank() != 0 ) {
        commBuff.resize( commBuff.givePackSize(MPI_INT, 1) );
        commBuff.packInt(boundaryDofManNum);

#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Sending data to partition 0 (send %d)\n",
                        this->giveEngngModel()->giveRank(),
                        "FETISolver :: setUpCommunicationMaps : send number of boundary dofMans", boundaryDofManNum);
#endif
        commBuff.iSend(0, FETICommunicator :: NumberOfBoundaryDofManagersMsg);

        MPI_Barrier(MPI_COMM_WORLD);
    }

#ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Barrier reached\n",
                    this->giveEngngModel()->giveRank(), "FETISolver :: setUpCommunicationMaps");
#endif




    commMap.resize(boundaryDofManNum);

    commBuff.resize(commBuff.givePackSize(MPI_INT, 1) * 2 * boundaryDofManNum);
    commBuff.init();
    // determine total number of DOFs per each boundary node
    // and build communication map simultaneously
    indx = 1;
    IntArray locNum;
    if ( domain->giveEngngModel()->giveRank() != 0 ) {
        for ( i = 1; i <= nnodes; i++ ) {
            if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
                // remember comm map entry
                commMap.at(indx++) = i;
                // determine number of DOFs
                domain->giveDofManager(i)->giveCompleteLocationArray(locNum, dn);
                neq = 0;
                for ( j = 1; j <= locNum.giveSize(); j++ ) {
                    if ( locNum.at(j) ) {
                        neq++;
                    }
                }

                commBuff.packInt( domain->giveDofManager(i)->giveGlobalNumber() );
                commBuff.packInt(neq);
            }
        }
    } else {
        for ( i = 1; i <= nnodes; i++ ) {
            if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
                // remember comm map entry
                commMap.at(indx++) = i;
            }
        }
    }

    if ( domain->giveEngngModel()->giveRank() != 0 ) {
        processCommunicator.setToSendArry(this, commMap, 0);
        processCommunicator.setToRecvArry(this, commMap, 0);
        commBuff.iSend(0, FETICommunicator :: BoundaryDofManagersRecMsg);
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Barrier reached\n",
                        this->giveEngngModel()->giveRank(), "FETISolver :: setUpCommunicationMaps");
#endif
        MPI_Barrier(MPI_COMM_WORLD);
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Barrier reached\n",
                        this->giveEngngModel()->giveRank(), "FETISolver :: setUpCommunicationMaps");
#endif
        MPI_Barrier(MPI_COMM_WORLD);
    } else {
        // remember comm map for direct mapping at master
        this->masterCommMap = commMap;
    }
}





/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

void
FETISolver :: projection(FloatArray &v, FloatMatrix &l, FloatMatrix &l1)
/*
 * funkce provadi projekci v modifikovane metode sdruzenych gradientu
 *
 * vstupy
 * v - vektor pred projekci
 * l - matice L(m,n)
 * l1 - inverzni matice k matici L^T L, rozmery (n,n)
 *
 * vystup
 * v - vektor po projekci
 *
 * 7.12.1998
 */
{
    FloatArray help1, help2;

    help1.beTProductOf(l, v);
    help2.beProductOf(l1, help1);
    help1.beProductOf(l, help2);

    v.subtract(help1);
}




int
FETISolver :: packRBM(ProcessCommunicator &processComm)
{
    int result = 1;
    int i, ir, size;
    int j, ndofs, eqNum;
    IntArray const *toSendMap = processComm.giveToSendMap();
    CommunicationBuffer *send_buff = processComm.giveProcessCommunicatorBuff()->giveSendBuff();
    IntArray locationArray;
    EModelDefaultEquationNumbering dn;

    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        domain->giveDofManager( toSendMap->at(i) )->giveCompleteLocationArray(locationArray, dn);
        ndofs = locationArray.giveSize();
        for ( j = 1; j <= ndofs; j++ ) {
            if ( ( eqNum = locationArray.at(j) ) ) {
                for ( ir = 1; ir <= nse; ir++ ) {
                    result &= send_buff->packDouble( rbm.at(eqNum, ir) );
                }

                if ( nse == 0 ) {
                    result &= send_buff->packDouble(0.0);
                }
            }
        }
    }

    return result;
}


int
FETISolver :: masterUnpackRBM(ProcessCommunicator &processComm)
{
    int to, result = 1;
    int i, irbm, idof, size, receivedRank;
    int j, nshared, part, eqNum;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    CommunicationBuffer *recv_buff = processComm.giveProcessCommunicatorBuff()->giveRecvBuff();
    IntArray locationArray;
    double value;

    receivedRank = processComm.giveRank();
    size = toRecvMap->giveSize();
    if ( receivedRank != 0 ) {
        //  for (irbm = 1; irbm <= nsem.at(receivedRank+1); irbm++) {
        for ( i = 1; i <= size; i++ ) {
            to = toRecvMap->at(i);
            //
            // loop over all dofs
            for ( idof = 1; idof <= masterCommunicator->giveDofManager(to)->giveNumberOfDofs(); idof++ ) {
                for ( irbm = 1; irbm <= nsem.at(receivedRank + 1); irbm++ ) {
                    // unpack contribution
                    result &= recv_buff->unpackDouble(value);
                    if ( masterCommunicator->giveDofManager(to)->giveReferencePratition() == receivedRank ) { // contribution from reference partition localizes to all
                                                                                                              // corresponding DOFs in boundary node

                        // localize to corresponding places
                        nshared = masterCommunicator->giveDofManager(to)->giveNumberOfSharedPartitions();
                        for ( j = 1; j <= nshared; j++ ) {
                            part = masterCommunicator->giveDofManager(to)->giveSharedPartition(j);
                            if ( part == processComm.giveRank() ) {
                                continue;
                            }

                            eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(part, idof);
                            l.at(eqNum, rbmAddr.at(receivedRank + 1) + irbm - 1) = value;
                        }
                    } else { // no reference partition
                        eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(receivedRank, idof);
                        l.at(eqNum, rbmAddr.at(receivedRank + 1) + irbm - 1) = ( -1.0 ) * value;
                    }
                }
            }
        }
    }

    return result;
}

int
FETISolver :: masterMapRBM()
{
    int i, to, from, idof, irbm, receivedRank = 0, nshared, j, part, eqNum, result;
    IntArray locationArray;
    double value;
    int size = masterCommMap.giveSize();
    // master will map its own values directly
    int locpos;
    EModelDefaultEquationNumbering dn;

    for ( irbm = 1; irbm <= nsem.at(1); irbm++ ) {
        for ( i = 1; i <= size; i++ ) {
            to = masterCommunicator->giveMasterCommMapPtr()->at(i);
            // use receive map, send map is empty to prevent master to send
            // itself any data. Note, however, that send and receive maps are same.
            from = masterCommMap.at(i);

            domain->giveDofManager(from)->giveCompleteLocationArray(locationArray, dn);
            locpos = 1;
            //
            //   ndofs = locationArray.giveSize(); // including supported
            //   for (j=1; j<=ndofs; j++) {
            //   if (eqNum = locationArray.at(j)) {

            //
            // loop over all dofs
            for ( idof = 1; idof <= masterCommunicator->giveDofManager(to)->giveNumberOfDofs(); idof++, locpos++ ) {
                // extract source value
                while ( locationArray.at(locpos) == 0 ) {
                    locpos++;
                    // test if position of nonzero dof is allowed
                    if ( locpos > locationArray.giveSize() ) {
                        _error("Consistency dof error");
                    }
                }

                value = rbm.at(locationArray.at(locpos), irbm);
                if ( masterCommunicator->giveDofManager(to)->giveReferencePratition() == receivedRank ) { // contribution from reference partition localizes to all
                                                                                                          // corresponding DOFs in boundary node

                    // localize to corresponding places
                    nshared = masterCommunicator->giveDofManager(to)->giveNumberOfSharedPartitions();
                    for ( j = 1; j <= nshared; j++ ) {
                        part = masterCommunicator->giveDofManager(to)->giveSharedPartition(j);
                        if ( part == 0 ) {
                            continue;
                        }

                        eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(part, idof);
                        l.at(eqNum, rbmAddr.at(receivedRank + 1) + irbm - 1) = value;
                    }
                } else { // no reference partition
                    eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(receivedRank, idof);
                    l.at(eqNum, rbmAddr.at(receivedRank + 1) + irbm - 1) = ( -1.0 ) * value;
                }
            }
        }
    }

    result = 1;
    return result;
}


int
FETISolver :: packQQProducts(ProcessCommunicator &processComm)
{
    int result = 1;
    int i, size;
    IntArray const *toSendMap = processComm.giveToSendMap();
    CommunicationBuffer *send_buff = processComm.giveProcessCommunicatorBuff()->giveSendBuff();
    IntArray locationArray;


    size = toSendMap->giveSize();
    for ( i = 1; i <= nse; i++ ) {
        result &= send_buff->packDouble( qq.at(i) );
    }

    return result;
}


int
FETISolver :: masterUnpackQQProduct(ProcessCommunicator &processComm)
{
    int result = 1;
    int i, receivedRank;
    //IntArray const* toRecvMap = processComm.giveToRecvMap();
    CommunicationBuffer *recv_buff = processComm.giveProcessCommunicatorBuff()->giveRecvBuff();
    IntArray locationArray;
    //double value;

    receivedRank = processComm.giveRank();

    if ( receivedRank != 0 ) {
        for ( i = 1; i <= nsem.at(receivedRank + 1); i++ ) {
            result &= recv_buff->unpackDouble( q.at(rbmAddr.at(receivedRank + 1) + i - 1) );
        }
    }

    return result;
}


int
FETISolver :: masterMapQQProduct()
{
    int i, result = 0;

    for ( i = 1; i <= nsem.at(1); i++ ) {
        q.at(rbmAddr.at(1) + i - 1) = qq.at(i);
        result = 1;
    }

    return result;
}


int
FETISolver :: packSolution(ProcessCommunicator &processComm)
{
    // master

    int result = 1;
    int i, size;
    int j, k, ndofs, eqNum, nshared, part, from;
    double val;
    IntArray const *toSendMap = processComm.giveToSendMap();
    CommunicationBuffer *send_buff = processComm.giveProcessCommunicatorBuff()->giveSendBuff();
    IntArray locationArray;

    int rank = processComm.giveRank();
    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        from = toSendMap->at(i);
        ndofs = masterCommunicator->giveDofManager(from)->giveNumberOfDofs();
        if ( rank == masterCommunicator->giveDofManager(from)->giveReferencePratition() ) {
            // summ corresponding values (multipliers)
            nshared = masterCommunicator->giveDofManager(from)->giveNumberOfSharedPartitions();
            for ( k = 1; k <= ndofs; k++ ) {
                val = 0.0;
                for ( j = 1; j <= nshared; j++ ) {
                    part = masterCommunicator->giveDofManager(from)->giveSharedPartition(j);
                    if ( part == processComm.giveRank() ) {
                        continue;
                    }

                    eqNum = masterCommunicator->giveDofManager(from)->giveCodeNumber(part, k);
                    val += w.at(eqNum);
                }

                result &= send_buff->packDouble(val);
            }
        } else {
            masterCommunicator->giveDofManager(from)->giveCompleteLocationArray(rank, locationArray);
            for ( j = 1; j <= ndofs; j++ ) {
                if ( ( eqNum = locationArray.at(j) ) ) {
                    result &= send_buff->packDouble( ( -1.0 ) * w.at(eqNum) );
                }
            }
        }
    }

    return result;
}


int
FETISolver :: unpackSolution(ProcessCommunicator &processComm)
{
    // slaves unpack their slotion contributions
    int receivedRank, result = 1, to;
    int i, size;
    int j, ndofs, eqNum;
    double value;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    CommunicationBuffer *recv_buff = processComm.giveProcessCommunicatorBuff()->giveRecvBuff();
    IntArray locationArray;
    EModelDefaultEquationNumbering dn;

    receivedRank = processComm.giveRank();

    size = toRecvMap->giveSize();
    // if (receivedRank != 0) {
    for ( i = 1; i <= size; i++ ) {
        to = toRecvMap->at(i);
        domain->giveDofManager(to)->giveCompleteLocationArray(locationArray, dn);
        ndofs = locationArray.giveSize();
        for ( j = 1; j <= ndofs; j++ ) {
            if ( ( eqNum = locationArray.at(j) ) ) {
                result &= recv_buff->unpackDouble(value);
                dd.at(eqNum) = value;

#ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Unpacking solution value %f at %d\n",
                                this->giveEngngModel()->giveRank(), "FETISolver :: solveYourselfAt", value, eqNum);
#endif
            }
        }

        //  }
    }

    return result;
}

int
FETISolver :: masterMapSolution()
{
    int i, to, from, idof, receivedRank = 0, nshared, j, part, eqNum, result, locpos;
    IntArray locationArray;
    double value;
    int size = masterCommMap.giveSize();
    EModelDefaultEquationNumbering dn;

    for ( i = 1; i <= size; i++ ) {
        from = masterCommunicator->giveMasterCommMapPtr()->at(i);
        // use receive map, send map is empty to prevent master to send
        // itself any data. Note, however, that send and receive maps are same.
        to = masterCommMap.at(i);

        domain->giveDofManager(to)->giveCompleteLocationArray(locationArray, dn);
        locpos = 1;
        //
        // loop over all dofs
        for ( idof = 1; idof <= masterCommunicator->giveDofManager(from)->giveNumberOfDofs(); idof++, locpos++ ) {
            // extract source value
            while ( locationArray.at(locpos) == 0 ) {
                locpos++;
                // test if position of nonzero dof is allowed
                if ( locpos > locationArray.giveSize() ) {
                    _error("Consistency dof error");
                }
            }

            value = 0.0;
            if ( masterCommunicator->giveDofManager(from)->giveReferencePratition() == receivedRank ) { // contribution from reference partition localizes to all
                                                                                                        // corresponding DOFs in boundary node

                // localize to corresponding places
                nshared = masterCommunicator->giveDofManager(from)->giveNumberOfSharedPartitions();
                for ( j = 1; j <= nshared; j++ ) {
                    part = masterCommunicator->giveDofManager(from)->giveSharedPartition(j);
                    if ( part == 0 ) {
                        continue;
                    }

                    eqNum = masterCommunicator->giveDofManager(from)->giveCodeNumber(part, idof);
                    value += w.at(eqNum);
                }
            } else { // no reference partition
                eqNum = masterCommunicator->giveDofManager(from)->giveCodeNumber(receivedRank, idof);
                value = ( -1.0 ) * w.at(eqNum);
            }

            // unpack to locaL CORRESPONDING EQUATION
            dd.at( locationArray.at(locpos) ) = value;
        }
    }

    result = 1;
    return result;
}



int
FETISolver :: packResiduals(ProcessCommunicator &processComm)
{
    // slaves

    int result = 1;
    int i, size;
    int j, ndofs, eqNum;
    IntArray const *toSendMap = processComm.giveToSendMap();
    CommunicationBuffer *send_buff = processComm.giveProcessCommunicatorBuff()->giveSendBuff();
    IntArray locationArray;
    EModelDefaultEquationNumbering dn;

    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        domain->giveDofManager( toSendMap->at(i) )->giveCompleteLocationArray(locationArray, dn);
        ndofs = locationArray.giveSize();
        for ( j = 1; j <= ndofs; j++ ) {
            if ( ( eqNum = locationArray.at(j) ) ) {
                result &= send_buff->packDouble( pp.at(eqNum) );
            }
        }
    }

    return result;
}


int
FETISolver :: unpackResiduals(ProcessCommunicator &processComm)
{
    // master

    int result = 1;
    int i, size, receivedRank, to;
    int j, idof, nshared, part, eqNum;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    CommunicationBuffer *recv_buff = processComm.giveProcessCommunicatorBuff()->giveRecvBuff();
    IntArray locationArray;
    double value;

    receivedRank = processComm.giveRank();

    size = toRecvMap->giveSize();
    if ( receivedRank != 0 ) {
        for ( i = 1; i <= size; i++ ) {
            to = toRecvMap->at(i);
            //
            // loop over all dofs
            for ( idof = 1; idof <= masterCommunicator->giveDofManager(to)->giveNumberOfDofs(); idof++ ) {
                // unpack contribution
                result &= recv_buff->unpackDouble(value);
                if ( masterCommunicator->giveDofManager(to)->giveReferencePratition() == receivedRank ) { // contribution from reference partition localizes to all
                                                                                                          // corresponding DOFs in boundary node

                    // localize to corresponding places
                    nshared = masterCommunicator->giveDofManager(to)->giveNumberOfSharedPartitions();
                    for ( j = 1; j <= nshared; j++ ) {
                        part = masterCommunicator->giveDofManager(to)->giveSharedPartition(j);
                        if ( part == processComm.giveRank() ) {
                            continue;
                        }

                        eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(part, idof);
                        g.at(eqNum) += value;
                    }
                } else { // no reference partition
                    eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(receivedRank, idof);
                    g.at(eqNum) += ( -1.0 ) * value;
                }
            }
        }
    }

    return result;
}

int
FETISolver :: masterMapResiduals()
{ // master will map its own values directly
    int i, to, from, idof, receivedRank = 0, nshared, j, part, eqNum, result, locpos;
    IntArray locationArray;
    double value;
    int size = masterCommMap.giveSize();
    EModelDefaultEquationNumbering dn;

    for ( i = 1; i <= size; i++ ) {
        to = masterCommunicator->giveMasterCommMapPtr()->at(i);
        // use receive map, send map is empty to prevent master to send
        // itself any data. Note, however, that send and receive maps are same.
        from = masterCommMap.at(i);

        domain->giveDofManager(from)->giveCompleteLocationArray(locationArray, dn);
        locpos = 1;
        //
        //   ndofs = locationArray.giveSize(); // including supported
        //   for (j=1; j<=ndofs; j++) {
        //   if (eqNum = locationArray.at(j)) {

        //
        // loop over all dofs
        for ( idof = 1; idof <= masterCommunicator->giveDofManager(to)->giveNumberOfDofs(); idof++, locpos++ ) {
            // extract source value
            while ( locationArray.at(locpos) == 0 ) {
                locpos++;
                // test if position of nonzero dof is allowed
                if ( locpos > locationArray.giveSize() ) {
                    _error("Consistency dof error");
                }
            }

            value = pp.at( locationArray.at(locpos) );
            if ( masterCommunicator->giveDofManager(to)->giveReferencePratition() == receivedRank ) { // contribution from reference partition localizes to all
                                                                                                      // corresponding DOFs in boundary node

                // localize to corresponding places
                nshared = masterCommunicator->giveDofManager(to)->giveNumberOfSharedPartitions();
                for ( j = 1; j <= nshared; j++ ) {
                    part = masterCommunicator->giveDofManager(to)->giveSharedPartition(j);
                    if ( part == 0 ) {
                        continue;
                    }

                    eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(part, idof);
                    g.at(eqNum) += value;
                }
            } else { // no reference partition
                eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(receivedRank, idof);
                g.at(eqNum) += ( -1.0 ) * value;
            }
        }
    }

    result = 1;
    return result;
}




int
FETISolver :: packDirectionVector(ProcessCommunicator &processComm)
{
    // master

    int result = 1;
    int i, size;
    int j, k, ndofs, eqNum, nshared, part, from;
    double val;
    IntArray const *toSendMap = processComm.giveToSendMap();
    CommunicationBuffer *send_buff = processComm.giveProcessCommunicatorBuff()->giveSendBuff();
    IntArray locationArray;

    int rank = processComm.giveRank();
    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        from = toSendMap->at(i);

        ndofs = masterCommunicator->giveDofManager(from)->giveNumberOfDofs();
        if ( rank == masterCommunicator->giveDofManager(from)->giveReferencePratition() ) {
            // summ corresponding values (multipliers)
            nshared = masterCommunicator->giveDofManager(from)->giveNumberOfSharedPartitions();
            for ( k = 1; k <= ndofs; k++ ) {
                val = 0.0;
                for ( j = 1; j <= nshared; j++ ) {
                    part = masterCommunicator->giveDofManager(from)->giveSharedPartition(j);
                    if ( part == processComm.giveRank() ) {
                        continue;
                    }

                    eqNum = masterCommunicator->giveDofManager(from)->giveCodeNumber(part, k);
                    val += d.at(eqNum);
                }

                result &= send_buff->packDouble(val);
            }
        } else {
            masterCommunicator->giveDofManager(from)->giveCompleteLocationArray(rank, locationArray);
            for ( j = 1; j <= ndofs; j++ ) {
                if ( ( eqNum = locationArray.at(j) ) ) {
                    result &= send_buff->packDouble( ( -1.0 ) * d.at(eqNum) );
                }
            }
        }
    }

    return result;
}


int
FETISolver :: unpackDirectionVector(ProcessCommunicator &processComm)
{
    // slaves unpack their slotion contributions
    int result = 1;
    int i, size;
    int j, ndofs, eqNum;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    CommunicationBuffer *recv_buff = processComm.giveProcessCommunicatorBuff()->giveRecvBuff();
    IntArray locationArray;
    // int receivedRank = domainComm.giveRank();
    EModelDefaultEquationNumbering dn;

    size = toRecvMap->giveSize();
    // if (receivedRank != 0) {
    for ( i = 1; i <= size; i++ ) {
        domain->giveDofManager( toRecvMap->at(i) )->giveCompleteLocationArray(locationArray, dn);
        ndofs = locationArray.giveSize();
        for ( j = 1; j <= ndofs; j++ ) {
            if ( ( eqNum = locationArray.at(j) ) ) {
                result &= recv_buff->unpackDouble( dd.at(eqNum) );
            }
        }
    }


    // }
    return result;
}

int
FETISolver :: masterMapDirectionVector()
{
    int i, to, from, idof, receivedRank = 0, nshared, j, part, eqNum, result, locpos;
    IntArray locationArray;
    double value;
    int size = masterCommMap.giveSize();
    EModelDefaultEquationNumbering dn;

    for ( i = 1; i <= size; i++ ) {
        from = masterCommunicator->giveMasterCommMapPtr()->at(i);
        // use receive map, send map is empty to prevent master to send
        // itself any data. Note, however, that send and receive maps are same.
        to = masterCommMap.at(i);

        domain->giveDofManager(to)->giveCompleteLocationArray(locationArray, dn);
        locpos = 1;
        //
        // loop over all dofs
        for ( idof = 1; idof <= masterCommunicator->giveDofManager(from)->giveNumberOfDofs(); idof++, locpos++ ) {
            // extract source value
            while ( locationArray.at(locpos) == 0 ) {
                locpos++;
                // test if position of nonzero dof is allowed
                if ( locpos > locationArray.giveSize() ) {
                    _error("Consistency dof error");
                }
            }

            value = 0.0;
            if ( masterCommunicator->giveDofManager(from)->giveReferencePratition() == receivedRank ) { // contribution from reference partition localizes to all
                                                                                                        // corresponding DOFs in boundary node

                // localize to corresponding places
                nshared = masterCommunicator->giveDofManager(from)->giveNumberOfSharedPartitions();
                for ( j = 1; j <= nshared; j++ ) {
                    part = masterCommunicator->giveDofManager(from)->giveSharedPartition(j);
                    if ( part == 0 ) {
                        continue;
                    }

                    eqNum = masterCommunicator->giveDofManager(from)->giveCodeNumber(part, idof);
                    value += d.at(eqNum);
                }
            } else { // no reference partition
                eqNum = masterCommunicator->giveDofManager(from)->giveCodeNumber(receivedRank, idof);
                value = ( -1.0 ) * d.at(eqNum);
            }

            // unpack to locaL CORRESPONDING EQUATION
            dd.at( locationArray.at(locpos) ) = value;
        }
    }

    result = 1;
    return result;
}


int
FETISolver :: packPPVector(ProcessCommunicator &processComm)
{
    // slaves

    int result = 1;
    int i, size;
    int j, ndofs, eqNum;
    IntArray const *toSendMap = processComm.giveToSendMap();
    CommunicationBuffer *send_buff = processComm.giveProcessCommunicatorBuff()->giveSendBuff();
    IntArray locationArray;
    EModelDefaultEquationNumbering dn;

    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        domain->giveDofManager( toSendMap->at(i) )->giveCompleteLocationArray(locationArray, dn);
        ndofs = locationArray.giveSize();
        for ( j = 1; j <= ndofs; j++ ) {
            if ( ( eqNum = locationArray.at(j) ) ) {
                result &= send_buff->packDouble( pp.at(eqNum) );
            }
        }
    }

    return result;
}


int
FETISolver :: unpackPPVector(ProcessCommunicator &processComm)
{
    // master

    int result = 1;
    int i, size, receivedRank, to;
    int j, idof, nshared, part, eqNum;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    CommunicationBuffer *recv_buff = processComm.giveProcessCommunicatorBuff()->giveRecvBuff();
    IntArray locationArray;
    double value;

    receivedRank = processComm.giveRank();

    size = toRecvMap->giveSize();
    if ( receivedRank != 0 ) {
        for ( i = 1; i <= size; i++ ) {
            to = toRecvMap->at(i);
            //
            // loop over all dofs
            for ( idof = 1; idof <= masterCommunicator->giveDofManager(to)->giveNumberOfDofs(); idof++ ) {
                // unpack contribution
                result &= recv_buff->unpackDouble(value);
                if ( masterCommunicator->giveDofManager(to)->giveReferencePratition() == receivedRank ) { // contribution from reference partition localizes to all
                                                                                                          // corresponding DOFs in boundary node

                    // localize to corresponding places
                    nshared = masterCommunicator->giveDofManager(to)->giveNumberOfSharedPartitions();
                    for ( j = 1; j <= nshared; j++ ) {
                        part = masterCommunicator->giveDofManager(to)->giveSharedPartition(j);
                        if ( part == processComm.giveRank() ) {
                            continue;
                        }

                        eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(part, idof);
                        p.at(eqNum) += value;
                    }
                } else { // no reference partition
                    eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(receivedRank, idof);
                    p.at(eqNum) += ( -1.0 ) * value;
                }
            }
        }
    }

    return result;
}

int
FETISolver :: masterMapPPVector()
{ // master will map its own values directly
    int i, to, from, idof, receivedRank = 0, nshared, j, part, eqNum, result, locpos;
    IntArray locationArray;
    double value;
    int size = masterCommMap.giveSize();
    EModelDefaultEquationNumbering dn;

    for ( i = 1; i <= size; i++ ) {
        to = masterCommunicator->giveMasterCommMapPtr()->at(i);
        // use receive map, send map is empty to prevent master to send
        // itself any data. Note, however, that send and receive maps are same.
        from = masterCommMap.at(i);

        domain->giveDofManager(from)->giveCompleteLocationArray(locationArray, dn);
        locpos = 1;
        //
        //   ndofs = locationArray.giveSize(); // including supported
        //   for (j=1; j<=ndofs; j++) {
        //   if (eqNum = locationArray.at(j)) {

        //
        // loop over all dofs
        for ( idof = 1; idof <= masterCommunicator->giveDofManager(to)->giveNumberOfDofs(); idof++, locpos++ ) {
            // extract source value
            while ( locationArray.at(locpos) == 0 ) {
                locpos++;
                // test if position of nonzero dof is allowed
                if ( locpos > locationArray.giveSize() ) {
                    _error("Consistency dof error");
                }
            }

            value = pp.at( locationArray.at(locpos) );
            if ( masterCommunicator->giveDofManager(to)->giveReferencePratition() == receivedRank ) { // contribution from reference partition localizes to all
                                                                                                      // corresponding DOFs in boundary node

                // localize to corresponding places
                nshared = masterCommunicator->giveDofManager(to)->giveNumberOfSharedPartitions();
                for ( j = 1; j <= nshared; j++ ) {
                    part = masterCommunicator->giveDofManager(to)->giveSharedPartition(j);
                    if ( part == 0 ) {
                        continue;
                    }

                    eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(part, idof);
                    p.at(eqNum) += value;
                }
            } else { // no reference partition
                eqNum = masterCommunicator->giveDofManager(to)->giveCodeNumber(receivedRank, idof);
                p.at(eqNum) += ( -1.0 ) * value;
            }
        }
    }

    result = 1;
    return result;
}


int
FETISolver :: packGammas(ProcessCommunicator &processComm)
{
    // master

    int irbm, result = 1;
    int rank = processComm.giveRank();
    CommunicationBuffer *send_buff = processComm.giveProcessCommunicatorBuff()->giveSendBuff();

    for ( irbm = 1; irbm <= nsem.at(rank + 1); irbm++ ) {
        result &= send_buff->packDouble( gamma.at(rbmAddr.at(rank + 1) + irbm - 1) );
    }

    return result;
}


int
FETISolver :: unpackGammas(ProcessCommunicator &processComm)
{
    // slaves
    int irbm, result = 1;
    CommunicationBuffer *recv_buff = processComm.giveProcessCommunicatorBuff()->giveRecvBuff();

    localGammas.resize(nse);
    for ( irbm = 1; irbm <= nse; irbm++ ) {
        result &= recv_buff->unpackDouble( localGammas.at(irbm) );
    }

    return result;
}

int
FETISolver :: masterMapGammas()
{
    // slaves
    int irbm, result = 1;

    localGammas.resize(nse);
    for ( irbm = 1; irbm <= nse; irbm++ ) {
        localGammas.at(irbm) = gamma.at(rbmAddr.at(1) + irbm - 1);
    }

    return result;
}


NM_Status
FETISolver :: solve(SparseMtrx *A, FloatArray *partitionLoad, FloatArray *partitionSolution)
{
    int i, j, tnse = 0, rank = domain->giveEngngModel()->giveRank();
    int ani;
    int source, tag;
    int masterLoopStatus;
    double nom = 0.0, denom, alpha, beta, ares, energyNorm = 0.0;
    FloatMatrix l1;
    StaticCommunicationBuffer commBuff(MPI_COMM_WORLD);
    Skyline *partitionStiffness;

    if ( A->giveType() != SMT_Skyline ) {
        _error("solve: unsuported sparse matrix type");
    }

    partitionStiffness = ( Skyline * ) A;

    if ( ( partitionSolution->giveSize() ) != partitionLoad->giveSize() ) {
        _error("solveYourselfAt: size mismatch");
    }

    int neq = partitionStiffness->giveNumberOfRows();
    int size = domain->giveEngngModel()->giveNumberOfProcesses();
    nsem.resize(size);

    // initialize communication maps
    this->setUpCommunicationMaps();
    if ( rank == 0 ) {
        this->commBuff = new CommunicatorBuff(size);
        masterCommunicator = new FETICommunicator(engngModel, this->commBuff, engngModel->giveRank(), size);
        masterCommunicator->setUpCommunicationMaps(engngModel);
    }

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Barrier reached\n",
                    this->giveEngngModel()->giveRank(), "FETISolver :: solveYourselfAt");
#endif

    /*********************************/
    /*  rozklad matice               */
    /*  vypocet rigid body modes     */
    /*********************************/
    /*  pole indexu zavislych rovnic  */
    se.resize(6);
    se.zero();
    /*  pole posunu jako tuheho telesa  */
    // rbm.resize (neq, 6); rbm.zero();

#ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: rbmodes computation startup, lneq is %d\n",
                    rank, "FETISolver :: solveYourselfAt", neq);
#endif

    partitionStiffness->rbmodes(rbm, nse, se, limit, 3);

#ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Number of Rigid body modes %3d\n",
                    rank, "FETISolver :: solveYourselfAt", nse);

    //fprintf(stdout, "\n[process rank %3d]: %-30s: RBM dump",
    //    rank,"FETISolver :: solveYourselfAt");
    //rbm.printYourself();

#endif

    /*****************************/
    /*  vypocet soucinu R^T . f  */
    /*****************************/

    if ( nse != 0 ) {
        qq.resize(nse);
        qq.zero();

        qq.beTProductOf(rbm, * partitionLoad);
        qq.times(-1.0);
    }

    /***************************************************************************/
    /*  zjisteni rozmeru matice L, vektoru q, poctu neznamych na podoblastech  */
    /***************************************************************************/

    //  commBuff.init();
    if ( rank == 0 ) {
        commBuff.resize( commBuff.givePackSize(MPI_INT, 1) );
        //
        // receive data
        //
        for ( i = 1; i < size; i++ ) {
            commBuff.iRecv(MPI_ANY_SOURCE, FETISolver :: NumberOfRBMMsg);
            while ( !commBuff.testCompletion(source, tag) ) {
                ;
            }

            // unpack data
#ifdef __VERBOSE_PARALLEL
            OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Received data from partition %3d\n",
                            rank, "FETICommunicator :: setUpCommunicationMaps : received number of partition rbm", source);
#endif
            commBuff.init();
            commBuff.unpackInt( nsem.at(source + 1) );
            tnse += nsem.at(source + 1);
        }

        nsem.at(1) = nse;
        tnse += nse;

        OOFEM_LOG_INFO("Number of RBM per partion\npart. rbm\n-------------------------------\n");
        for ( i = 1; i <= size; i++ ) {
            OOFEM_LOG_INFO( "%-4d %8d\n", i - 1, nsem.at(i) );
        }


        // initialize partition start indexes into rbm array
        rbmAddr.resize(size);
        rbmAddr.zero();
        rbmAddr.at(1) = 1;
        for ( i = 2; i <= size; i++ ) {
            rbmAddr.at(i) = rbmAddr.at(i - 1) + nsem.at(i - 1);
        }
    } else { // slave code
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Sending number of Rigid body modes %3d\n",
                        rank, "FETISolver :: solveYourselfAt", nse);
#endif

        commBuff.resize( commBuff.givePackSize(MPI_INT, 1) );
        commBuff.packInt(nse);
        commBuff.iSend(0, FETISolver :: NumberOfRBMMsg);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*****************************/
    /*  sestaveni cele matice L  */
    /*****************************/

    if ( rank == 0 ) {
        // assemble l matrix containing rbm localized from partition contibutions
        if ( tnse ) {
            l.resize(masterCommunicator->giveNumberOfEquations(), tnse);
            l.zero();
        }

        masterCommunicator->initReceive(FETISolver :: RBMMessage);
        masterCommunicator->unpackAllData(this, & FETISolver :: masterUnpackRBM);
        this->masterMapRBM();

        if ( tnse ) {
            l.times(-1.0);
        }
    } else {
        // pack RBMs
        processCommunicator.packData(this, & FETISolver :: packRBM);
        processCommunicator.initSend(FETISolver :: RBMMessage);
        MPI_Barrier(MPI_COMM_WORLD);
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: RBMMessage Barrier finished",
                        rank, "FETISolver :: solveYourselfAt");
#endif
    }

    if ( rank == 0 ) {
        // receive products of rbm * LoadVector
        if ( tnse ) {
            q.resize(tnse);
            q.zero();
        }

        masterCommunicator->initReceive(FETISolver :: QQMessage);
        masterCommunicator->unpackAllData(this, & FETISolver :: masterUnpackQQProduct);
        this->masterMapQQProduct();
    } else {
        // pack QQ products
        processCommunicator.packData(this, & FETISolver :: packQQProducts);
        processCommunicator.initSend(FETISolver :: QQMessage);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /*******************************************/
    /*******************************************/
    /*  sestaveni matice L^T.L a jeji inverze  */
    /*******************************************/
    /*******************************************/

    if ( ( rank == 0 ) && ( tnse != 0 ) ) {
        FloatMatrix l2;
        // compute l1, the inverse of (l^Tl)
        l2.beTProductOf(l, l);
        l1.beInverseOf(l2);
    }

    /****************************************************************************/
    /****************************************************************************/
    /*  vypocet pocatecni aproximace neznamych do modifikovane met. sdr. grad.  */
    /****************************************************************************/
    /****************************************************************************/

    /****************************/
    /*  alokace nekterych poli  */
    /****************************/
    /*  vektor smeru na jedne podoblasti  */
    dd.resize(neq);
    /*  pomocny vektor na jedne podoblasti  */
    pp.resize(neq);


    if ( rank == 0 ) {
        /*  vektor smeru  */
        d.resize(masterCommunicator->giveNumberOfEquations() + 1);
        /*  pomocny vektor  */
        p.resize( masterCommunicator->giveNumberOfEquations() );
        /*  vektor gradientu  */
        g.resize( masterCommunicator->giveNumberOfEquations() );
        /*  vektor neznamych  */
        w.resize( masterCommunicator->giveNumberOfEquations() );

        /*  soucin (L^T.L)^{-1}.q  */
        if ( tnse ) {
            FloatArray help(tnse);
            help.beProductOf(l1, q);
            w.beProductOf(l, help);
        } else {
            w.zero();
        }

        // send aproximation of solution to slaves
        masterCommunicator->packAllData(this, & FETISolver :: packSolution);
        masterCommunicator->initSend(FETISolver :: SolutionMessage);
        this->masterMapSolution();
    } else {
        // receive send aproximation of solution
        processCommunicator.initReceive(FETISolver :: SolutionMessage);
        while ( !processCommunicator.receiveCompleted() ) {
            ;
        }

        processCommunicator.unpackData(this, & FETISolver :: unpackSolution);
    }

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Solution Approx Barrier finished\n",
                    rank, "FETISolver :: solveYourselfAt");
#endif

    dd.subtract(*partitionLoad);
    partitionStiffness->ldl_feti_sky(pp, dd, nse, limit, se);

    if ( rank == 0 ) {
        // receive contributions (to resudals) from slaves

        masterCommunicator->initReceive(FETISolver :: ResidualMessage);
        masterCommunicator->unpackAllData(this, & FETISolver :: unpackResiduals);
        this->masterMapResiduals();
    } else {
        // send contributions (to resudals) to master
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Residual contribution packing initiated\n",
                        rank, "FETISolver :: solveYourselfAt");
#endif
        processCommunicator.packData(this, & FETISolver :: packResiduals);
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Residual contribution send initiated\n",
                        rank, "FETISolver :: solveYourselfAt");
#endif
        processCommunicator.initSend(FETISolver :: ResidualMessage);

        MPI_Barrier(MPI_COMM_WORLD);
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Residual contribution Barrier finished\n",
                        rank, "FETISolver :: solveYourselfAt");
#endif
    }


    /*  v tuto chvili je vypocten nulty vektor gradientu  */

    /***********************************/
    /*  vypocet nulteho vektoru smeru  */
    /***********************************/
    if ( rank == 0 ) {
        if ( tnse ) {
            this->projection(g, l, l1);
        }

        d = g;
        d.negated();

        /*  vypocet nulteho citatele  */
        nom = g.computeSquaredNorm();

        //#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("\nIteration process\n");
        if ( energyNorm_comput_flag ) {
            OOFEM_LOG_DEBUG("iteration      gradient vector norm      energy norm\n=====================================================================\n");
        } else {
            OOFEM_LOG_DEBUG("iteration      gradient vector norm\n================================================\n");
        }

        //#endif
    }

    /*  cyklus vyjadrujici pocet iteraci  */
    for ( i = 0; i < ni; i++ ) {
        dd.zero();
        if ( rank == 0 ) {
            /***************************************/
            /*  vypocet pomocneho vektoru K.d = p  */
            /***************************************/

            // send direction vector to slaves
            masterCommunicator->packAllData(this, & FETISolver :: packDirectionVector);
            masterCommunicator->initSend(FETISolver :: DirectionVectorMessage);
            this->masterMapDirectionVector();
        } else {
            // receive direction vector on slaves
            processCommunicator.initReceive(FETISolver :: DirectionVectorMessage);
            while ( !processCommunicator.receiveCompleted() ) {
                ;
            }

            processCommunicator.unpackData(this, & FETISolver :: unpackDirectionVector);
        }

        MPI_Barrier(MPI_COMM_WORLD);


        pp.zero();
        partitionStiffness->ldl_feti_sky(pp, dd, nse, limit, se);

        if ( rank == 0 ) {
            p.zero();

            // receive contributions (to resudals) from slaves

            masterCommunicator->initReceive(FETISolver :: PPVectorMessage);
            masterCommunicator->unpackAllData(this, & FETISolver :: unpackPPVector);
            this->masterMapPPVector();
        } else {
            // pack pp on all partitions
            processCommunicator.packData(this, & FETISolver :: packPPVector);
            processCommunicator.initSend(FETISolver :: PPVectorMessage);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        if ( rank == 0 ) {
            /*  calculation of the denominator of the fraction defining Alpha */
            denom = d.dotProduct(p);

            if ( fabs(denom) < FETISOLVER_ZERONUM ) {
                OOFEM_LOG_RELEVANT("FETISolver::solve :  v modifikovane metode sdruzenych gradientu je nulovy jmenovatel u soucinitele alpha\n");
                commBuff.resize( commBuff.givePackSize(MPI_INT, 1) );
                commBuff.init();
                commBuff.packInt(FETISolverIterationBreak);
                commBuff.bcast(0);
                break;
            }

            /*******************/
            /*  vypocet alpha  */
            /*******************/
            alpha = nom / denom;

            /**************************************************************/
            /*  vypocet noveho gradientu g a nove aproximace neznamych x  */
            /**************************************************************/
            for ( j = 1; j <= masterCommunicator->giveNumberOfEquations(); j++ ) {
                w.at(j) += alpha * d.at(j);
                g.at(j) += alpha * p.at(j);
            }



            if ( tnse ) {
                this->projection(g, l, l1);
            }

            denom = nom;
            if ( fabs(denom) < FETISOLVER_ZERONUM ) {
                OOFEM_LOG_RELEVANT("FETISolver::solve : v modifikovane metode sdruzenych gradientu je nulovy jmenovatel u soucinitele beta\n");
                commBuff.resize( commBuff.givePackSize(MPI_INT, 1) );
                commBuff.init();
                commBuff.packInt(FETISolverIterationBreak);
                commBuff.bcast(0);
                break;
            }

            nom = g.computeSquaredNorm();

            if ( nom < err ) {
                commBuff.resize( commBuff.givePackSize(MPI_INT, 1) );
                commBuff.init();
                commBuff.packInt(FETISolverIterationBreak);
                commBuff.bcast(0);
                break;
            }


            /******************/
            /*  vypocet beta  */
            /******************/
            beta = nom / denom;

            /****************************/
            /*  vypocet noveho smeru d  */
            /****************************/
            for ( j = 1; j <= masterCommunicator->giveNumberOfEquations(); j++ ) {
                d.at(j) = beta * d.at(j) - g.at(j);
            }
        }


        if ( rank == 0 ) {
            commBuff.resize( commBuff.givePackSize(MPI_INT, 1) );
            commBuff.init();
            commBuff.packInt(FETISolverIterationContinue);
            commBuff.bcast(0);
        } else {
            commBuff.resize( commBuff.givePackSize(MPI_INT, 1) );
            commBuff.init();
            commBuff.bcast(0);
            commBuff.unpackInt(masterLoopStatus);
            if ( masterLoopStatus == FETISolverIterationBreak ) {
#ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Received Loop break signal from master\n",
                                rank, "FETISolver :: solveYourselfAt");
#endif

                break;
            }
        }


        // computation of energy norm
        if ( energyNorm_comput_flag ) {
            dd.zero();
            if ( rank == 0 ) {
                // send aproximation of solution to slaves
                masterCommunicator->packAllData(this, & FETISolver :: packSolution);
                masterCommunicator->initSend(FETISolver :: SolutionMessage);
                this->masterMapSolution();
            } else {
                // receive send aproximation of solution
                processCommunicator.initReceive(FETISolver :: SolutionMessage);
                while ( !processCommunicator.receiveCompleted() ) {
                    ;
                }

                processCommunicator.unpackData(this, & FETISolver :: unpackSolution);
            }

            MPI_Barrier(MPI_COMM_WORLD);

            pp.zero();
            partitionStiffness->ldl_feti_sky(pp, dd, nse, limit, se);

            if ( rank == 0 ) {
                p.zero();

                // receive contributions (to resudals) from slaves

                masterCommunicator->initReceive(FETISolver :: PPVectorMessage);
                masterCommunicator->unpackAllData(this, & FETISolver :: unpackPPVector);
                this->masterMapPPVector();
            } else {
                // pack pp on all partitions
                processCommunicator.packData(this, & FETISolver :: packPPVector);
                processCommunicator.initSend(FETISolver :: PPVectorMessage);
                MPI_Barrier(MPI_COMM_WORLD);
            }

            if ( rank == 0 ) {
                energyNorm = p.dotProduct(w);
            }
        }

        // end of energy norm computation

        //#ifdef __VERBOSE_PARALLEL
        if ( rank == 0 ) {
            if ( energyNorm_comput_flag ) {
                OOFEM_LOG_DEBUG("%-9d%15e %15e\n", i, nom, energyNorm);
            } else {
                OOFEM_LOG_DEBUG("%-9d%15e\n", i, nom);
            }
        }

        //#endif
    } // end iterative loop

    if ( rank == 0 ) {
        ani = i;
        ares = nom;
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("Konec metody sdruzenych gradientu, nite %d, err %e\n", ani, ares);
#endif
    }

    /*  prelom (break) */

    /****************************************************/
    /*  konec modifikovane metody sdruzenych gradientu  */
    /****************************************************/

    if ( rank == 0 ) {
        OOFEM_LOG_INFO("End of iteration, reached norm %15e\n", nom);
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("\nVysledne Lagrangeovy multiplikatory\n");
        for ( i = 1; i <= masterCommunicator->giveNumberOfEquations(); i++ ) {
            OOFEM_LOG_DEBUG( "lambda %4d          %f\n", i, w.at(i) );
        }

#endif
    }

    dd.zero();
    if ( rank == 0 ) {
        // send solution to slaves
        masterCommunicator->packAllData(this, & FETISolver :: packSolution);
        masterCommunicator->initSend(FETISolver :: SolutionMessage);
        this->masterMapSolution();
    } else {
        // receive solution
        processCommunicator.initReceive(FETISolver :: SolutionMessage);
        while ( !processCommunicator.receiveCompleted() ) {
            ;
        }

        processCommunicator.unpackData(this, & FETISolver :: unpackSolution);
    }

    MPI_Barrier(MPI_COMM_WORLD);


#ifdef __VERBOSE_PARALLEL
    //fprintf(stdout, "\n[process rank %3d]: %-30s: Partition lagrange multipliers dump:\n",
    //    rank,"FETISolver :: solveYourselfAt");
    //dd.printYourself();
#endif

    partitionLoad->subtract(dd);
#ifdef __VERBOSE_PARALLEL
    // fprintf(stdout, "\n[process rank %3d]: %-30s: Partition load dump:\n",
    //    rank,"FETISolver :: solveYourselfAt");
    //partitionLoad->printYourself();
#endif
    partitionStiffness->ldl_feti_sky(* partitionSolution, * partitionLoad, nse, limit, se);
    pp = * partitionSolution;

#ifdef __VERBOSE_PARALLEL
    //fprintf(stdout, "\n[process rank %3d]: %-30s: Partition solution dump:\n",
    //    rank,"FETISolver :: solveYourselfAt");
    //pp.printYourself();
#endif

    if ( rank == 0 ) {
        g.zero();
        // receive contributions (to resudals) from slaves

        masterCommunicator->initReceive(FETISolver :: ResidualMessage);
        masterCommunicator->unpackAllData(this, & FETISolver :: unpackResiduals);

#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Received residuals from slaves\n",
                        rank, "FETISolver :: solveYourselfAt");
#endif

        this->masterMapResiduals();
    } else {
        // send contributions (to resudals) to master

#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Sending residuals to master\n",
                        rank, "FETISolver :: solveYourselfAt");
#endif

        processCommunicator.packData(this, & FETISolver :: packResiduals);
        processCommunicator.initSend(FETISolver :: ResidualMessage);
        MPI_Barrier(MPI_COMM_WORLD);
    }

#ifdef __VERBOSE_PARALLEL
    if ( rank == 0 ) {
        //fprintf(stdout, "\n[process rank %3d]: %-30s: Partition residuals dump:\n",
        //    rank,"FETISolver :: solveYourselfAt");
        //g.printYourself();
    }

#endif

    if ( rank == 0 ) {
        if ( tnse ) {
            FloatArray help1;
            help1.beTProductOf(l, g);
            // coefficient of linear combinations
            gamma.beProductOf(l1, help1);
        } else {
            gamma.zero();
        }

        // send corresponding gammas to slaves
        masterCommunicator->packAllData(this, & FETISolver :: packGammas);
        masterCommunicator->initSend(FETISolver :: GammasMessage);
        this->masterMapGammas();
    } else {
        // receive gammas
        processCommunicator.initReceive(FETISolver :: GammasMessage);
        while ( !processCommunicator.receiveCompleted() ) {
            ;
        }

        processCommunicator.unpackData(this, & FETISolver :: unpackGammas);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if ( nse != 0 ) {
        FloatArray help;
        help.beProductOf(rbm, localGammas);

        for ( j = 1; j <= neq; j++ ) {
            partitionSolution->at(j) += help.at(j);
        }
    }

    return NM_Success;
}
} // end namespace oofem
#endif
