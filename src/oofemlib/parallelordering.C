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

#include "parallelordering.h"
#include "engngm.h"
#include "domain.h"
#include "unknownnumberingscheme.h"
#include "combuff.h"
#include "mathfem.h"

namespace oofem {
bool
ParallelOrdering :: isLocal(DofManager *dman)
{
    int myrank = dman->giveDomain()->giveEngngModel()->giveRank();
    if ( dman->giveParallelMode() == DofManager_local ) {
        return true;
    }

    if ( dman->giveParallelMode() == DofManager_shared ) {
        // determine if problem is the lowest one sharing the dofman; if yes the receiver is responsible to
        // deliver number
        return myrank <= dman->givePartitionList()->minimum();
    }

    return false;
}

bool
ParallelOrdering :: isShared(DofManager *dman)
{
    if ( dman->giveParallelMode() == DofManager_shared ) {
        return true;
    } else {
        return false;
    }
}

Natural2GlobalOrdering :: Natural2GlobalOrdering() : ParallelOrdering(), locGlobMap(), globLocMap()
{
    l_neqs = g_neqs = 0;
}


void
Natural2GlobalOrdering :: init(EngngModel *emodel, int di, const UnknownNumberingScheme &n)
{
    Domain *d = emodel->giveDomain(di);
    int ndofman = d->giveNumberOfDofManagers();
    int myrank = emodel->giveRank();
    DofManager *dman;
    // determine number of local eqs + number of those shared DOFs which are numbered by receiver
    // shared dofman is numbered on partition with lovest rank number

#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("Natural2GlobalOrdering :: init", "initializing N2G ordering", myrank);
#endif

    l_neqs = 0;
    for ( int i = 1; i <= ndofman; i++ ) {
        dman = d->giveDofManager(i);
#if 0
        if ( dman->giveParallelMode() == DofManager_local ) {  // count all dofman eqs
            for ( Dof *dof: *dman ) {
                if ( dof->isPrimaryDof() ) {
                    if ( dof->giveEquationNumber() ) {
                        l_neqs++;
                    }
                }
            }
        } else if ( dman->giveParallelMode() == DofManager_shared ) {
            // determine if problem is the lowest one sharing the dofman; if yes the receiver is responsible to
            // deliver number
            IntArray *plist = dman->givePartitionList();
            int n = plist->giveSize();
            int minrank = myrank;
            for ( j = 1; j <= n; j++ ) {
                minrank = min( minrank, plist->at(j) );
            }
            if ( minrank == myrank ) { // count eqs
                for ( Dof *dof: *dman ) {
                    if ( dof->isPrimaryDof() ) {
                        if ( dof->giveEquationNumber() ) {
                            l_neqs++;
                        }
                    }
                }
            }
        }   // end shared dman
#endif
        if ( isLocal(dman) ) {
            for ( Dof *dof: *dman ) {
                if ( dof->isPrimaryDof() ) {
                    if ( dof->giveEquationNumber(n) ) {
                        l_neqs++;
                    }
                }
            }
        }
    }

    // exchange with other procs the number of eqs numbered on particular procs
    int *leqs = new int [ emodel->giveNumberOfProcesses() ];
    MPI_Allgather(& l_neqs, 1, MPI_INT, leqs, 1, MPI_INT, MPI_COMM_WORLD);
    // compute local offset
    int offset = 0;
    for ( int j = 0; j < myrank; j++ ) {
        offset += leqs [ j ];
    }

    // count global number of eqs
    g_neqs = 0;
    for ( int j = 0; j < emodel->giveNumberOfProcesses(); j++ ) {
        g_neqs += leqs [ j ];
    }

    // send numbered shared ones
    locGlobMap.resize( emodel->giveNumberOfDomainEquations(di, n) );

    // determine shared dofs
    int psize, nproc = emodel->giveNumberOfProcesses();
    IntArray sizeToSend(nproc), sizeToRecv(nproc), nrecToReceive(nproc);
#ifdef __VERBOSE_PARALLEL
    IntArray nrecToSend(nproc);
#endif
    const IntArray *plist;
    for ( int i = 1; i <= ndofman; i++ ) {
        // if (domain->giveDofManager(i)->giveParallelMode() == DofManager_shared) {
        if ( isShared( d->giveDofManager(i) ) ) {
            int n = d->giveDofManager(i)->giveNumberOfDofs();
            plist = d->giveDofManager(i)->givePartitionList();
            psize = plist->giveSize();
            int minrank = myrank;
            for ( int j = 1; j <= psize; j++ ) {
                minrank = min( minrank, plist->at(j) );
            }

            if ( minrank == myrank ) { // count to send
                for ( int j = 1; j <= psize; j++ ) {
#ifdef __VERBOSE_PARALLEL
                    nrecToSend( plist->at(j) )++;
#endif
                    sizeToSend( plist->at(j) ) += ( 1 + n );  // ndofs+dofman number
                }
            } else {
                nrecToReceive(minrank)++;
                sizeToRecv(minrank) += ( 1 + n );      // ndofs+dofman number
            }
        }
    }

#ifdef __VERBOSE_PARALLEL
    for ( int i = 0; i < nproc; i++ ) {
        OOFEM_LOG_INFO("[%d] Record Statistics: Sending %d Receiving %d to %d\n",
                       myrank, nrecToSend(i), nrecToReceive(i), i);
    }

#endif


    std :: map< int, int >globloc; //  global->local mapping for shared
    // number local guys
    int globeq = offset;
    for ( int i = 1; i <= ndofman; i++ ) {
        dman = d->giveDofManager(i);
        //if (dman->giveParallelMode() == DofManager_shared) {
        if ( isShared(dman) ) {
            globloc [ dman->giveGlobalNumber() ] = i; // build global->local mapping for shared

            plist = dman->givePartitionList();
            psize = plist->giveSize();
            int minrank = myrank;
            for ( int j = 1; j <= psize; j++ ) {
                minrank = min( minrank, plist->at(j) );
            }

            if ( minrank == myrank ) { // local
                for ( Dof *dof: *dman ) {
                    if ( dof->isPrimaryDof() ) {
                        int eq = dof->giveEquationNumber(n);

                        if ( eq ) {
                            locGlobMap.at(eq) = globeq++;
                        }
                    }
                }
            }

            //} else if (dman->giveParallelMode() == DofManager_local) {
        } else {
            for ( Dof *dof: *dman ) {
                if ( dof->isPrimaryDof() ) {
                    int eq = dof->giveEquationNumber(n);

                    if ( eq ) {
                        locGlobMap.at(eq) = globeq++;
                    }
                }
            }
        }
    }


    /*
     * fprintf (stderr, "[%d] locGlobMap: ", myrank);
     * for (i=1; i<=locGlobMap.giveSize(); i++)
     * fprintf (stderr, "%d ",locGlobMap.at(i));
     */

    // pack data for remote procs
    CommunicationBuffer **buffs = new CommunicationBuffer * [ nproc ];
    for ( int p = 0; p < nproc; p++ ) {
        buffs [ p ] = new StaticCommunicationBuffer(MPI_COMM_WORLD, 0);
        buffs [ p ]->resize( buffs [ p ]->givePackSizeOfInt(1) * sizeToSend(p) );

#if 0
        OOFEM_LOG_INFO( "[%d]Natural2GlobalOrdering :: init: Send buffer[%d] size %d\n",
                       myrank, p, sizeToSend(p) );
#endif
    }


    for ( int i = 1; i <= ndofman; i++ ) {
        if ( isShared( d->giveDofManager(i) ) ) {
            dman = d->giveDofManager(i);
            plist = dman->givePartitionList();
            psize = plist->giveSize();
            int minrank = myrank;
            for ( int j = 1; j <= psize; j++ ) {
                minrank = min( minrank, plist->at(j) );
            }

            if ( minrank == myrank ) { // do send
                for ( int j = 1; j <= psize; j++ ) {
                    int p = plist->at(j);
                    if ( p == myrank ) {
                        continue;
                    }

#if 0
                    OOFEM_LOG_INFO("[%d]Natural2GlobalOrdering :: init: Sending localShared node %d[%d] to proc %d\n",
                                   myrank, i, dman->giveGlobalNumber(), p);
#endif
                    buffs [ p ]->write( dman->giveGlobalNumber() );
                    for ( Dof *dof: *dman ) {
                        if ( dof->isPrimaryDof() ) {
                            int eq = dof->giveEquationNumber(n);

                            if ( eq ) {
                                buffs [ p ]->write( locGlobMap.at(eq) );
                            }
                        }
                    }
                }
            }
        }
    }


    //fprintf (stderr, "[%d] Sending glob nums ...", myrank);
    // send buffers
    for ( int p = 0; p < nproc; p++ ) {
        if ( p != myrank ) {
            buffs [ p ]->iSend(p, 999);
        }
    }


#if 0
    for ( int p = 0; p < nproc; p++ ) {
        if ( p == myrank ) {
            continue;
        }
        for ( int i = 1; i <= ndofman; i++ ) {
            //if (domain->giveDofManager(i)->giveParallelMode() == DofManager_shared) {
            if ( isShared( d->giveDofManager(i) ) ) {
                dman = d->giveDofManager(i);
                plist = dman->givePartitionList();
                psize = plist->giveSize();
                int minrank = myrank;
                for ( j = 1; j <= psize; j++ ) {
                    minrank = min( minrank, plist->at(j) );
                }
                if ( minrank == myrank ) { // do send
                    buffs [ p ]->write( dman->giveGlobalNumber() );
                    for ( Dof *dof: *dman ) {
                        if ( dof->isPrimaryDof() ) {
                            buffs [ p ]->write( locGlobMap.at( dof->giveEquationNumber() ) );
                        }
                    }
                }
            }
        }
        // send buffer
        buffs [ p ]->iSend(p, 999);
    }
#endif

    // receive remote eqs and complete global numbering
    CommunicationBuffer **rbuffs = new CommunicationBuffer * [ nproc ];
    for ( int p = 0; p < nproc; p++ ) {
        rbuffs [ p ] = new StaticCommunicationBuffer(MPI_COMM_WORLD, 0);
        rbuffs [ p ]->resize( rbuffs [ p ]->givePackSizeOfInt(1) * sizeToRecv(p) );
#if 0
        OOFEM_LOG_INFO( "[%d]Natural2GlobalOrdering :: init: Receive buffer[%d] size %d\n",
                       myrank, p, sizeToRecv(p) );
#endif
    }


    //fprintf (stderr, "[%d] Receiving glob nums ...", myrank);
    for ( int p = 0; p < nproc; p++ ) {
        if ( p != myrank ) {
            rbuffs [ p ]->iRecv(p, 999);
        }
    }


    IntArray finished(nproc);
    finished.zero();
    int fin = 1;
    finished.at(emodel->giveRank() + 1) = 1;
    do {
        for ( int p = 0; p < nproc; p++ ) {
            if ( finished.at(p + 1) == 0 ) {
                if ( rbuffs [ p ]->testCompletion() ) {
                    // data are here
                    // unpack them
                    int nite = nrecToReceive(p);
                    int shdm, ldm;
                    for ( int i = 1; i <= nite; i++ ) {
                        rbuffs [ p ]->read(shdm);

#if 0
                        OOFEM_LOG_INFO("[%d]Natural2GlobalOrdering :: init: Received shared node [%d] from proc %d\n",
                                       myrank, shdm, p);
#endif
                        //
                        // find local guy coorecponding to shdm
                        if ( globloc.find(shdm) != globloc.end() ) {
                            ldm = globloc [ shdm ];
                        } else {
                          continue;
                            OOFEM_ERROR("[%d] invalid shared dofman received, globnum %d\n", myrank, shdm);
                            ldm = 0;
                        }

                        dman = d->giveDofManager(ldm);
                        for ( Dof *dof: *dman ) {
                            if ( dof->isPrimaryDof() ) {
                                int eq = dof->giveEquationNumber(n);

                                if ( eq ) {
                                    int val;
                                    rbuffs [ p ]->read(val);
                                    locGlobMap.at(eq) = val;
                                }
                            }
                        }
                    }

                    finished.at(p + 1) = 1;
                    fin++;
                }
            }
        }
    } while ( fin < nproc );


    /*
     * fprintf (stderr, "[%d] Finished receiving glob nums ...", myrank);
     *
     * fprintf (stderr, "[%d] locGlobMap:", myrank);
     * for (i=1; i<=locGlobMap.giveSize(); i++)
     * fprintf (stderr, "%d ",locGlobMap.at(i));
     */

#ifdef  __VERBOSE_PARALLEL
    int _eq;
    char *ptr;
    char locname[] = "local", shname[] = "shared", unkname[] = "unknown";
    for ( int i = 1; i <= ndofman; i++ ) {
      dman = d->giveDofManager(i);
      if ( dman->giveParallelMode() == DofManager_local ) {
        ptr = locname;
      } else if ( dman->giveParallelMode() == DofManager_shared ) {
        ptr = shname;
      } else {
        ptr = unkname;
      }
      
      for ( Dof *dof: *dman ) {
        DofIDItem id = dof->giveDofID();
        if ( ( _eq = dof->giveEquationNumber(n) ) ) {
          fprintf( stderr, "[%d] n:%6s %d[%d] (%d), leq = %d, geq = %d\n", emodel->giveRank(), ptr, i, dman->giveGlobalNumber(), id, _eq, locGlobMap.at(_eq) );
        } else {
          fprintf(stderr, "[%d] n:%6s %d[%d] (%d), leq = %d, geq = %d\n", emodel->giveRank(), ptr, i, dman->giveGlobalNumber(), id, _eq, 0);
        }
      }
    }

#endif


    // build reverse map
    int lneq = emodel->giveNumberOfDomainEquations(di, n);

    globLocMap.clear();
    for ( int i = 1; i <= lneq; i++ ) {
        globLocMap [ locGlobMap.at(i) ] = i;
    }

    for ( int p = 0; p < nproc; p++ ) {
        delete rbuffs [ p ];
        delete buffs [ p ];
    }

    delete[] rbuffs;
    delete[] buffs;
    delete[] leqs;

    MPI_Barrier(MPI_COMM_WORLD);
#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("Natural2GlobalOrdering :: init", "done", myrank);
#endif
}


int
Natural2GlobalOrdering :: giveNewEq(int leq)
{
    return locGlobMap.at(leq);
}

int
Natural2GlobalOrdering :: giveOldEq(int eq)
{
    std :: map< int, int > :: iterator i = globLocMap.find(eq);
    if ( i != globLocMap.end() ) {
        return i->second;
    } else {
        return 0;
    }
}

void
Natural2GlobalOrdering :: map2New(IntArray &answer, const IntArray &src, int baseOffset)
{
    int indx, n = src.giveSize();
    answer.resize(n);

    //for (i=1; i<=n; i++) answer.at(i)=locGlobMap.at(src.at(i))-baseOffset;
    for ( int i = 1; i <= n; i++ ) {
        if ( ( indx = src.at(i) ) ) {
            answer.at(i) = locGlobMap.at(indx) + baseOffset;
        } else {
            answer.at(i) = ( -1 ) + baseOffset;
        }
    }
}


void
Natural2GlobalOrdering :: map2Old(IntArray &answer, const IntArray &src, int baseOffset)
{
    int n = src.giveSize();
    int offset = baseOffset - 1;
    answer.resize(n);
    for ( int i = 1; i <= n; i++ ) {
        answer.at(i) = giveOldEq( src.at(i) ) + offset;
    }
}




Natural2LocalOrdering :: Natural2LocalOrdering() : ParallelOrdering(), n2l() { }

void
Natural2LocalOrdering :: init(EngngModel *emodel, int di, const UnknownNumberingScheme &n)
{
    Domain *d = emodel->giveDomain(di);
    int ndofman = d->giveNumberOfDofManagers(), loc_eq = 1;

    // determine number of local eqs + number of those shared DOFs which are numbered by receiver
    // shared dofman is numbered on partition with lovest rank number

    n2l.resize( emodel->giveNumberOfDomainEquations(di, n) );

    for ( int i = 1; i <= ndofman; i++ ) {
        DofManager *dman = d->giveDofManager(i);
        bool lFlag = isLocal(dman);
        for ( Dof *dof: *dman ) {
            if ( dof->isPrimaryDof() ) {
                int n_eq = dof->giveEquationNumber(n);

                if ( n_eq == 0 ) {
                    continue;
                }

                if ( lFlag ) {
                    n2l.at(n_eq) = loc_eq++;
                } else {
                    n2l.at(n_eq) = 0;
                }
            }
        }
    }
}

int
Natural2LocalOrdering :: giveNewEq(int leq)
{
    return n2l.at(leq);
}

int
Natural2LocalOrdering :: giveOldEq(int eq)
{
    // not really efficient
    // it is assumed that queries in oposite directuion take only place
    // this is here only for completness, but performance is really bad.
    return ( n2l.findFirstIndexOf(eq) );
}

void
Natural2LocalOrdering :: map2New(IntArray &answer, const IntArray &src, int baseOffset)
{
    int indx, n = src.giveSize();
    answer.resize(n);
    for ( int i = 1; i <= n; i++ ) {
        if ( ( indx = src.at(i) ) ) {
            answer.at(i) = n2l.at(indx) - baseOffset;
        } else {
            answer.at(i) = 0 - baseOffset;
        }
    }
}

void
Natural2LocalOrdering :: map2Old(IntArray &answer, const IntArray &src, int baseOffset)
{
    int n = src.giveSize();
    answer.resize(n);
    for ( int i = 1; i <= n; i++ ) {
        answer.at(i) = giveOldEq( src.at(i) ) - baseOffset;
    }
}
} // end namespace oofem
