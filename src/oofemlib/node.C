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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#include "node.h"
#include "dof.h"
#include "unknowntype.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "nodload.h"
#include "timestep.h"

#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "verbose.h"
#include "datastream.h"
#include "contextioerr.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "enrichmentitem.h"
#endif

namespace oofem {
Node :: Node(int n, Domain *aDomain) :
    DofManager(n, aDomain), coordinates()
{
    localCoordinateSystem = NULL;
}


Node :: ~Node()
{
    delete localCoordinateSystem;
}


double
Node :: giveCoordinate(int i)
// Returns the i-th coordinate of the receiver.
{
    if ( i > coordinates.giveSize() ) {
        return 0.;
    }

    return coordinates.at(i);
}


IRResultType Node :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int j, size;
    FloatArray triplets;

#  ifdef VERBOSE
    // VERBOSE_PRINT1("Instanciating node ",number)
#  endif

    DofManager :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, coordinates, IFT_Node_coords, "coords"); // Macro

    //
    // scaling of coordinates if necessary
    //
    if ( domain->giveEngngModel()->giveEquationScalingFlag() ) {
        double lscale = domain->giveEngngModel()->giveVariableScale(VST_Length);
        this->coordinates.times(1. / lscale);
    }


    // Read if available local coordinate system in this node
    triplets.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, triplets, IFT_Node_lcs, "lcs"); // Macro
    size = triplets.giveSize();
    if ( !( ( size == 0 ) || ( size == 6 ) ) ) {
        _warning2( "initializeFrom: lcs in node %d is not properly defined, will be ignored", this->giveNumber() );
    }

    if ( size == 6 ) {
        double n1 = 0.0, n2 = 0.0;
        localCoordinateSystem = new FloatMatrix(3, 3);

        for ( j = 1; j <= 3; j++ ) {
            localCoordinateSystem->at(1, j) = triplets.at(j);
            n1 += triplets.at(j) * triplets.at(j);
            localCoordinateSystem->at(2, j) = triplets.at(j + 3);
            n2 += triplets.at(j + 3) * triplets.at(j + 3);
        }

        n1 = sqrt(n1);
        n2 = sqrt(n2);
        if ( ( n1 <= 1.e-6 ) || ( n2 <= 1.e-6 ) ) {
            _error("instanciateFrom : lcs input error");
        }

        for ( j = 1; j <= 3; j++ ) { // normalize e1' e2'
            localCoordinateSystem->at(1, j) /= n1;
            localCoordinateSystem->at(2, j) /= n2;
        }

        // vector e3' computed from vector product of e1', e2'
        localCoordinateSystem->at(3, 1) =
            localCoordinateSystem->at(1, 2) * localCoordinateSystem->at(2, 3) -
            localCoordinateSystem->at(1, 3) * localCoordinateSystem->at(2, 2);
        localCoordinateSystem->at(3, 2) =
            localCoordinateSystem->at(1, 3) * localCoordinateSystem->at(2, 1) -
            localCoordinateSystem->at(1, 1) * localCoordinateSystem->at(2, 3);
        localCoordinateSystem->at(3, 3) =
            localCoordinateSystem->at(1, 1) * localCoordinateSystem->at(2, 2) -
            localCoordinateSystem->at(1, 2) * localCoordinateSystem->at(2, 1);
    }

    return IRRT_OK;
}

void
Node :: computeLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the vector of the nodal loads of the receiver.
{
    int i, n, nLoads;
    NodalLoad *loadN;
    FloatArray contribution;
    IntArray dofIDarry(0);
    FloatMatrix L2G;

    if ( this->giveLoadArray()->isEmpty() ) {
        answer.resize(0);
        return;
    } else {
        answer.resize(0);
        nLoads = loadArray.giveSize();       // the node may be subjected
        for ( i = 1; i <= nLoads; i++ ) {     // to more than one load
            n     = loadArray.at(i);
            loadN = dynamic_cast< NodalLoad * >( domain->giveLoad(n) );
            if ( !loadN ) {
                _error("computeLoadVectorAt: incompatible load");
            }

            if ( loadN->giveBCGeoType() != NodalLoadBGT ) {
                _error("computeLoadVectorAt: incompatible load type applied");
            }

            loadN->computeComponentArrayAt(contribution, stepN, mode); // can be NULL
            // Transform from Global to Local c.s.
            if ( loadN->giveCoordSystMode() == NodalLoad :: BL_GlobalMode ) {
                if ( this->computeL2GTransformation(L2G, dofIDarry) ) {
                    contribution.rotatedWith(L2G, 't');
                }
            }
            answer.add(contribution);
        }
    }
}


void
Node :: printYourself()
// Prints the receiver on screen.
{
    int i;
    double x, y;

    x = this->giveCoordinate(1);
    y = this->giveCoordinate(2);
    printf("Node %d    coord : x %f  y %f\n", number, x, y);
    for ( i = 0; i < numberOfDofs; i++ ) {
        if ( dofArray [ i ] ) {
            dofArray [ i ]->printYourself();
        } else {
            printf("dof %d is nil \n", i + 1);
        }
    }

    loadArray.printYourself();
    printf("\n");
}


void
Node :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    int i, ic;

    fMode mode = domain->giveEngngModel()->giveFormulation();

    double dt = tStep->giveTimeIncrement();

    for ( i = 1; i <= numberOfDofs; i++ ) {
        if ( mode == AL ) { // updated Lagrange
            ic = domain->giveCorrespondingCoordinateIndex(i);
            if ( ic != 0 ) {
                Dof *d = this->giveDof(i);
                DofIDItem id = d->giveDofID();
                if ( id == D_u || id == D_v || id == D_w ) {
                    coordinates.at(ic) += d->giveUnknown(EID_MomentumBalance, VM_Incremental, tStep);
                } else if ( id == V_u || id == V_v || id == V_w ) {
                    coordinates.at(ic) += d->giveUnknown(EID_MomentumBalance, VM_Total, tStep) * dt;
                }
            }
        }
    }
}


double
Node :: giveUpdatedCoordinate(int ic, TimeStep *tStep, EquationID type, double scale)
//
// returns coordinate + scale * displacement
// displacement is of updMode (UpdateMode) type
//
{
    int i, j;
    FloatMatrix *T;

    if ( ( ic < 1 ) || ( ic > 3 ) ) {
        _error("giveUpdatedCoordinate: Can't return non-existing coordinate (index not in range 1..3)");
        return 0.;
    }

    if ( tStep->isTheCurrentTimeStep() ) {
        double coordinate = this->giveCoordinate(ic);
        if ( !this->hasLocalCS() ) {
            // this has no local cs.
            for ( i = 1; i <= numberOfDofs; i++ ) {
                j = domain->giveCorrespondingCoordinateIndex(i);
                if ( ( j != 0 ) && ( j == ic ) ) {
                    coordinate +=
                        scale * this->giveDof(i)->giveUnknown(type, VM_Total, tStep);
                    break;
                }
            }
        } else {
            //
            // this has local cs.
            // We must perform transformation of displacements DOFs
            // in to global c.s and then to add them to global coordinates.
            //
            T = this->giveLocalCoordinateTriplet();
            FloatArray displacements(3);
            for ( i = 1; i <= 3; i++ ) {
                displacements.at(i) = 0.;
            }

            for ( i = 1; i <= numberOfDofs; i++ ) {
                j = domain->giveCorrespondingCoordinateIndex(i);
                if ( j ) { // && (this->giveDof(i)->giveUnknownType()==DisplacementVector))
                    displacements.at(j) = scale * this->giveDof(i)->
                                          giveUnknown(type, VM_Total, tStep);
                }
            }

            // perform transformation for desired displacement
            for ( i = 1; i <= 3; i++ ) {
                coordinate += displacements.at(i) * T->at(i, ic);
            }
        }

        return coordinate;
    } else {
        _error("Can't return updatedCoordinate for non-current timestep");
    }

    return 0.;
}

int
Node :: checkConsistency()
{
    /*
     * Checks internal data consistency in node.
     * Current implementation checks (when receiver has slave dofs) if receiver has the same
     * coordinate system as master dofManager of slave dof.
     */
    int result = 1;
    int ndofs = this->giveNumberOfDofs();
    int i, nslaves = 0;

    result = result && DofManager :: checkConsistency();

    for ( i = 1; i <= ndofs; i++ ) {
        if ( this->giveDof(i)->giveClassID() == SimpleSlaveDofClass ) {
            nslaves++;
        }
    }

    if ( nslaves == 0 ) {
        return result;             // return o.k. if no slaves exists
    }

    IntArray masterDofManagers(nslaves);
    int numberOfMDM = 0; // counter of diferent master dofManagers
    int j, master, alreadyFound = 0;
    Dof *idof;
    Node *masterNode;

    for ( i = 1; i <= ndofs; i++ ) {
        if ( ( idof = this->giveDof(i) )->giveClassID() == SimpleSlaveDofClass ) {
            alreadyFound  = 0;
            master = ( ( SimpleSlaveDof * ) idof )->giveMasterDofManagerNum();
            for ( j = 1; j <= numberOfMDM; j++ ) {
                if ( masterDofManagers.at(j) == master ) {
                    alreadyFound = 1;
                    break;
                }
            }

            if ( alreadyFound == 0 ) {
                // check master for same coordinate system
                // first mark master as checked
                numberOfMDM++;
                masterDofManagers.at(numberOfMDM) = master;
                // compare coordinate systems
                masterNode = dynamic_cast< Node * >( domain->giveDofManager(master) );
                if ( !masterNode ) {
                    _warning2("checkConsistency: master dofManager is not compatible", 1);
                    result = 0;
                } else if ( !this->hasSameLCS(masterNode) ) {
                    _warning2("checkConsistency: different lcs for master/slave nodes", 1);
                    result = 0;
                }
            }
        }
    }

    return result;
}


bool
Node :: hasSameLCS(Node *remote)
{
    FloatMatrix *thisLcs, *masterLcs;
    thisLcs = this->giveLocalCoordinateTriplet();
    masterLcs = remote->giveLocalCoordinateTriplet();
    int k, l;

    if ( ( this->hasLocalCS() ) && ( remote->hasLocalCS() ) ) {
        for ( k = 1; k <= 3; k++ ) {
            for ( l = 1; l <= 3; l++ ) {
                if ( fabs( thisLcs->at(k, l) - masterLcs->at(k, l) ) > 1.e-4 ) {
                    return false;
                }
            }
        }
    } else if ( this->hasLocalCS() ) {
        for ( k = 1; k <= 3; k++ ) {
            for ( l = 1; l <= 3; l++ ) {
                if ( fabs( thisLcs->at(k, l) - ( k == l ) ) > 1.e-4 ) {
                    return false;
                }
            }
        }
    } else if ( remote->hasLocalCS() ) {
        for ( k = 1; k <= 3; k++ ) {
            for ( l = 1; l <= 3; l++ ) {
                if ( fabs( masterLcs->at(k, l) - ( k == l ) ) > 1.e-4 ) {
                    return false;
                }
            }
        }
    }

    return true;
}


bool
Node :: computeL2GTransformation(FloatMatrix &answer, const IntArray &dofIDArry)
{
    //
    // computes transformation of receiver from global cs to nodal (user-defined) cs.
    // Note: implementation rely on D_u, D_v and D_w (R_u, R_v, R_w) order in cltypes.h
    // file. Do not change their order and do not insert any values between these values.
    //

    int i, j;
    DofIDItem id, id2;

    if ( localCoordinateSystem == NULL ) {
        answer.beEmptyMtrx();
        return false;
    } else {
        if ( dofIDArry.isEmpty() ) {
            // response for all local dofs is computed

            answer.resize(numberOfDofs, numberOfDofs);
            answer.zero();

            for ( i = 1; i <= numberOfDofs; i++ ) {
                // test for vector quantities
                switch ( id = giveDof(i)->giveDofID() ) {
                case D_u:
                case D_v:
                case D_w:
                    for ( j = 1; j <= numberOfDofs; j++ ) {
                        id2 = giveDof(j)->giveDofID();
                        if ( ( id2 == D_u ) || ( id2 == D_v ) || ( id2 == D_w ) ) {
                            answer.at(j, i) = localCoordinateSystem->at( ( int ) ( id ) - ( int ) ( D_u ) + 1,
                                                                        ( int ) ( id2 ) - ( int ) ( D_u ) + 1 );
                        }
                    }

                    break;

                case V_u:
                case V_v:
                case V_w:
                    for ( j = 1; j <= numberOfDofs; j++ ) {
                        id2 = giveDof(j)->giveDofID();
                        if ( ( id2 == V_u ) || ( id2 == V_v ) || ( id2 == V_w ) ) {
                            answer.at(j, i) = localCoordinateSystem->at( ( int ) ( id ) - ( int ) ( V_u ) + 1,
                                                                        ( int ) ( id2 ) - ( int ) ( V_u ) + 1 );
                        }
                    }

                    break;

                case R_u:
                case R_v:
                case R_w:
                    for ( j = 1; j <= numberOfDofs; j++ ) {
                        id2 = giveDof(j)->giveDofID();
                        if ( ( id2 == R_u ) || ( id2 == R_v ) || ( id2 == R_w ) ) {
                            answer.at(j, i) = localCoordinateSystem->at( ( int ) ( id ) - ( int ) ( R_u ) + 1,
                                                                        ( int ) ( id2 ) - ( int ) ( R_u ) + 1 );
                        }
                    }

                    break;

                case T_f:
                case P_f:
                    // scalar quantities
                    answer.at(i, i) = 1.0;
                    break;

                default:
                    _error2( "computeGNTransformation: unknown dofID (%s)", __DofIDItemToString(id) );
                }
            }
        } else { // end if (dofIDArry.isEmpty())
            // map is provided -> assemble for requested dofs
            int size = dofIDArry.giveSize();
            answer.resize(size, size);
            answer.zero();

            for ( i = 1; i <= size; i++ ) {
                // test for vector quantities
                switch ( id = ( DofIDItem ) dofIDArry.at(i) ) {
                case D_u:
                case D_v:
                case D_w:
                    for ( j = 1; j <= size; j++ ) {
                        id2 = ( DofIDItem ) dofIDArry.at(j);
                        if ( ( id2 == D_u ) || ( id2 == D_v ) || ( id2 == D_w ) ) {
                            answer.at(j, i) = localCoordinateSystem->at( ( int ) ( id ) - ( int ) ( D_u ) + 1, ( int ) ( id2 ) - ( int ) ( D_u ) + 1 );
                        }
                    }

                    break;

                case V_u:
                case V_v:
                case V_w:
                    for ( j = 1; j <= size; j++ ) {
                        id2 = ( DofIDItem ) dofIDArry.at(j);
                        if ( ( id2 == V_u ) || ( id2 == V_v ) || ( id2 == V_w ) ) {
                            answer.at(j, i) = localCoordinateSystem->at( ( int ) ( id ) - ( int ) ( V_u ) + 1, ( int ) ( id2 ) - ( int ) ( V_u ) + 1 );
                        }
                    }

                    break;

                case R_u:
                case R_v:
                case R_w:
                    for ( j = 1; j <= size; j++ ) {
                        id2 = ( DofIDItem ) dofIDArry.at(j);
                        if ( ( id2 == R_u ) || ( id2 == R_v ) || ( id2 == R_w ) ) {
                            answer.at(j, i) = localCoordinateSystem->at( ( int ) ( id ) - ( int ) ( R_u ) + 1, ( int ) ( id2 ) - ( int ) ( R_u ) + 1 );
                        }
                    }

                    break;

                case T_f:
                case P_f:
                    // scalar quantities
                    answer.at(i, i) = 1.0;
                    break;

                default:
                    _error2( "computeGNTransformation: unknown dofID (%s)", __DofIDItemToString(id) );
                }
            }
        } // end map is provided -> assemble for requested dofs

    } // end localCoordinateSystem defined
    return true;
}


contextIOResultType
Node :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = DofManager :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        int _haslcs = hasLocalCS();
        if ( ( iores = coordinates.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream->write(& _haslcs, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( _haslcs ) {
            if ( ( iores = localCoordinateSystem->storeYourself(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        }
    }

    return CIO_OK;
}


contextIOResultType
Node :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("restoreContex : can't write into NULL stream");
    }

    if ( ( iores = DofManager :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        int _haslcs;
        if ( ( iores = coordinates.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream->read(& _haslcs, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( _haslcs ) {
            if ( localCoordinateSystem == NULL ) {
                localCoordinateSystem = new FloatMatrix();
            }

            if ( ( iores = localCoordinateSystem->restoreYourself(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        } else {
            localCoordinateSystem = NULL;
        }
    }

    return CIO_OK;
}


#ifdef __OOFEG
void
Node :: drawYourself(oofegGraphicContext &gc)
//
// draws graphics representation of receiver
//
{
    GraphicObj *go;
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( ( mode == OGC_nodeGeometry ) || ( mode == OGC_nodeAnnotation ) ) {
        WCRec p [ 1 ]; /* point */
        p [ 0 ].x = ( FPNum ) this->giveCoordinate(1);
        p [ 0 ].y = ( FPNum ) this->giveCoordinate(2);
        p [ 0 ].z = ( FPNum ) this->giveCoordinate(3);

        EASValsSetLayer(OOFEG_NODE_ANNOTATION_LAYER);
        EASValsSetMType(FILLED_CIRCLE_MARKER);
 #if 1
        if ( this->giveDomain()->giveEngngModel()->hasXfemManager(1) ) {
            XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
            int i;
            for ( i = 1; i <= xf->giveNumberOfEnrichmentItems(); i++ ) {
                if ( xf->giveEnrichmentItem(i)->isDofManEnriched(this->number) ) {
                    EASValsSetMType(SQUARE_MARKER);
                }
            }
        }

 #endif

 #ifdef __PARALLEL_MODE
        if ( this->giveParallelMode() == DofManager_local ) {
            EASValsSetColor( gc.getNodeColor() );
        } else if ( this->giveParallelMode() == DofManager_shared ) {
            EASValsSetColor( gc.getDeformedElementColor() );
        } else {
            EASValsSetColor( gc.getCrackPatternColor() );
        }

 #else
        EASValsSetColor( gc.getNodeColor() );
 #endif

        bool ordinary = true;

        int idof;
        for ( idof = 1; idof <= this->giveNumberOfDofs(); idof++ ) {
            if ( this->giveDof(1)->giveClassID() == SimpleSlaveDofClass ) {
                ordinary = false;
                break;
            }
        }

        if ( !ordinary ) {
            EASValsSetColor( gc.getBcIcColor() );
        }

        EASValsSetMSize(8);
        go = CreateMarker3D(p);
        EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | MTYPE_MASK | MSIZE_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }

    if ( mode == OGC_nodeAnnotation ) {
        char num [ 30 ];
        WCRec p [ 1 ]; /* point */
        EASValsSetColor( gc.getNodeColor() );
        EASValsSetLayer(OOFEG_NODE_ANNOTATION_LAYER);
        p [ 0 ].x = ( FPNum ) this->giveCoordinate(1);
        p [ 0 ].y = ( FPNum ) this->giveCoordinate(2);
        p [ 0 ].z = ( FPNum ) this->giveCoordinate(3);
 #ifdef __PARALLEL_MODE
        sprintf( num, "%d(%d)", this->giveNumber(), this->giveGlobalNumber() );
 #else
        sprintf( num, "%d", this->giveLabel() );
 #endif
        go = CreateAnnText3D(p, num);
        EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    } else if ( mode == OGC_essentialBC ) {
        int i, hasDisplSupport [ 3 ], hasRotSupport [ 3 ], hasAny = 0;
        TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
        if ( !tStep ) {
            TimeStep __temp( domain->giveEngngModel() );
            tStep = & __temp;
        }

        WCRec pp [ 2 ];

        for ( i = 0; i < 3; i++ ) {
            hasDisplSupport [ i ] = 0;
            hasRotSupport [ i ] = 0;
        }

        for ( i = 1; i <= numberOfDofs; i++ ) {
            if ( this->giveDof(i)->hasBc(tStep) ) {
                hasAny = 1;
                switch ( giveDof(i)->giveDofID() ) {
                case D_u: hasDisplSupport [ 0 ] = 1;
                    break;
                case D_v: hasDisplSupport [ 1 ] = 1;
                    break;
                case D_w: hasDisplSupport [ 2 ] = 1;
                    break;
                case R_u: hasRotSupport [ 0 ] = 1;
                    break;
                case R_v: hasRotSupport [ 1 ] = 1;
                    break;
                case R_w: hasRotSupport [ 2 ] = 1;
                    break;
                default: break;
                }
            }
        }

        if ( hasAny != 0 ) {
            EASValsSetColor( gc.getBcIcColor() );
            EASValsSetLayer(OOFEG_BCIC_ANNOTATION_LAYER);
            pp [ 0 ].x = ( FPNum ) this->giveCoordinate(1);
            pp [ 0 ].y = ( FPNum ) this->giveCoordinate(2);
            pp [ 0 ].z = ( FPNum ) this->giveCoordinate(3);

            /* primary bc */
            for ( i = 0; i < 3; i++ ) {
                if ( hasDisplSupport [ i ] || hasRotSupport [ i ] ) {
                    pp [ 1 ].x = 0.;
                    pp [ 1 ].y = 0.;
                    pp [ 1 ].z = 0.;

                    if ( !this->hasLocalCS() ) {
                        if ( i == 0 ) {
                            pp [ 1 ].x = 1.0;
                        }

                        if ( i == 1 ) {
                            pp [ 1 ].y = 1.0;
                        }

                        if ( i == 2 ) {
                            pp [ 1 ].z = 1.0;
                        }
                    } else {
                        FloatMatrix *T = this->giveLocalCoordinateTriplet();
                        ;
                        pp [ 1 ].x = T->at(i + 1, 1);
                        pp [ 1 ].y = T->at(i + 1, 2);
                        pp [ 1 ].z = T->at(i + 1, 3);
                    }


                    if ( hasDisplSupport [ i ] && hasRotSupport [ i ] ) {
                        EASValsSetVecMType(TRIPLE_ARROW_VECMARKER);
                    }

                    if ( hasDisplSupport [ i ] ) {
                        EASValsSetVecMType(ARROW_VECMARKER);
                    }

                    if ( hasRotSupport [ i ] ) {
                        EASValsSetVecMType(DOUBLE_ARROW_VECMARKER);
                    }

                    go = CreateVecMarker3D(pp);
                    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | VECMTYPE_MASK, go);
                    EMAddGraphicsToModel(ESIModel(), go);
                }
            }
        }
    } else if ( mode == OGC_naturalBC ) {
        int i;
        TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
        if ( !tStep ) {
            TimeStep __temp( domain->giveEngngModel() );
            tStep = & __temp;
        }

        WCRec pp [ 2 ];
        /* load */
        if ( !this->giveLoadArray()->isEmpty() ) {
            double defScale = gc.getDefScale();
            pp [ 0 ].x = ( FPNum ) this->giveCoordinate(1);
            pp [ 0 ].y = ( FPNum ) this->giveCoordinate(2);
            pp [ 0 ].z = ( FPNum ) this->giveCoordinate(3);
            pp [ 1 ].x = pp [ 1 ].x = pp [ 1 ].x = 0.0;

            FloatArray load;
            FloatMatrix t;
            IntArray dofIDArry(0);
            computeLoadVectorAt(load, tStep, VM_Total);
            if (computeL2GTransformation(t, dofIDArry)) {
                load.rotatedWith(t,'n');
            }

            FloatArray force(3), momentum(3);
            for ( i = 1; i <= numberOfDofs; i++ ) {
                switch ( giveDof(i)->giveDofID() ) {
                case D_u: force.at(1) = defScale * load.at(i);
                    break;
                case D_v: force.at(2) = defScale * load.at(i);
                    break;
                case D_w: force.at(3) = defScale * load.at(i);
                    break;
                case R_u: momentum.at(1) = defScale * load.at(i);
                    break;
                case R_v: momentum.at(2) = defScale * load.at(i);
                    break;
                case R_w: momentum.at(3) = defScale * load.at(i);
                    break;
                default: break;
                }
            }

            EASValsSetColor( gc.getBcForceColor() );
            EASValsSetLayer(OOFEG_NATURALBC_LAYER);

            // draw force
            EASValsSetVecMType(ARROW_VECMARKER);
            pp [ 1 ].x = force.at(1);
            pp [ 1 ].y = force.at(2);
            pp [ 1 ].z = force.at(3);
            go = CreateVector3D(pp);
            EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | VECMTYPE_MASK, go);
            EMAddGraphicsToModel(ESIModel(), go);
            // draw moment
            EASValsSetVecMType(DOUBLE_ARROW_VECMARKER);
            pp [ 1 ].x = momentum.at(1);
            pp [ 1 ].y = momentum.at(2);
            pp [ 1 ].z = momentum.at(3);
            go = CreateVector3D(pp);
            EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | VECMTYPE_MASK, go);
            EMAddGraphicsToModel(ESIModel(), go);
        }
    } else if ( mode == OGC_nodeVectorPlot ) {
        GraphicObj *go;
        TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
        int i;

        if ( gc.giveIntVarType() == IST_Velocity ) {
            WCRec p [ 2 ]; /* point */
            double defScale = gc.getDefScale();

            p [ 0 ].x = p [ 1 ].x = ( FPNum ) this->giveCoordinate(1);
            p [ 0 ].y = p [ 1 ].y = ( FPNum ) this->giveCoordinate(2);
            p [ 0 ].z = p [ 1 ].z = ( FPNum ) this->giveCoordinate(3);

            //p[1].x = p[1].y = p[1].z = 0.0;
            for ( i = 1; i <= numberOfDofs; i++ ) {
                if ( this->giveDof(i)->giveDofID() == V_u ) {
                    p [ 1 ].x = defScale * this->giveDof(i)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                } else if ( this->giveDof(i)->giveDofID() == V_v ) {
                    p [ 1 ].y = defScale * this->giveDof(i)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                } else if ( this->giveDof(i)->giveDofID() == V_w ) {
                    p [ 1 ].z = defScale * this->giveDof(i)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                }
            }

            EASValsSetColor( gc.getDeformedElementColor() );
            EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
            go = CreateVector3D(p);
            EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK, go);
            EMAddGraphicsToModel(ESIModel(), go);
        }
    }
}

#endif
} // end namespace oofem
