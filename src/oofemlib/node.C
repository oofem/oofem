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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#include "node.h"
#include "dof.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "nodalload.h"
#include "timestep.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "verbose.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "domain.h"
#include "engngm.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "xfem/enrichmentitem.h"
 #include "xfem/xfemmanager.h"
#endif

namespace oofem {
REGISTER_DofManager(Node);

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
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int size;

#  ifdef VERBOSE
    // VERBOSE_PRINT1("Instanciating node ",number)
#  endif

    result = DofManager :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    IR_GIVE_FIELD(ir, coordinates, _IFT_Node_coords);

    //
    // scaling of coordinates if necessary
    //
    if ( domain->giveEngngModel()->giveEquationScalingFlag() ) {
        double lscale = domain->giveEngngModel()->giveVariableScale(VST_Length);
        this->coordinates.times(1. / lscale);
    }


    // Read if available local coordinate system in this node
    if ( ir->hasField(_IFT_Node_lcs) ) {
        FloatArray triplets;
        IR_GIVE_FIELD(ir, triplets, _IFT_Node_lcs);
        size = triplets.giveSize();
        if ( size != 6 ) {
            OOFEM_WARNING("lcs in node %d is not properly defined, will be ignored", this->giveNumber() );
        }

        double n1 = 0.0, n2 = 0.0;
        localCoordinateSystem = new FloatMatrix(3, 3);

        for ( int j = 1; j <= 3; j++ ) {
            localCoordinateSystem->at(1, j) = triplets.at(j);
            n1 += triplets.at(j) * triplets.at(j);
            localCoordinateSystem->at(2, j) = triplets.at(j + 3);
            n2 += triplets.at(j + 3) * triplets.at(j + 3);
        }

        n1 = sqrt(n1);
        n2 = sqrt(n2);
        if ( ( n1 <= 1.e-6 ) || ( n2 <= 1.e-6 ) ) {
            OOFEM_ERROR("lcs input error");
        }

        for ( int j = 1; j <= 3; j++ ) { // normalize e1' e2'
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

void Node :: giveInputRecord(DynamicInputRecord &input)
{
    DofManager :: giveInputRecord(input);

    input.setField(coordinates, _IFT_Node_coords);

    if ( localCoordinateSystem != NULL ) {
        input.setField(* localCoordinateSystem, _IFT_Node_lcs);
    }
}


void
Node :: computeLoadVector(FloatArray &answer, Load *load, CharType type, TimeStep *tStep, ValueModeType mode)
{
    answer.clear();
    if ( type != ExternalForcesVector ) {
        return;
    }

    NodalLoad *loadN = dynamic_cast< NodalLoad * >(load);
    if ( !loadN ) {
        OOFEM_ERROR("incompatible load");
    }

    if ( loadN->giveBCGeoType() != NodalLoadBGT ) {
        OOFEM_ERROR("incompatible load type applied");
    }

    load->computeComponentArrayAt(answer, tStep, mode);

    // Transform from Global to Local c.s.
    if ( loadN->giveCoordSystMode() == NodalLoad :: CST_Global ) {
        FloatMatrix L2G;
        if ( this->computeL2GTransformation(L2G, loadN->giveDofIDs()) ) {
            answer.rotatedWith(L2G, 't');
        }
    }
}


void
Node :: printYourself()
// Prints the receiver on screen.
{
    printf("Node %d    coord : x %f  y %fz< %f\n", number, this->giveCoordinate(1), this->giveCoordinate(2), this->giveCoordinate(3));
    for ( Dof *dof: *this ) {
        dof->printYourself();
    }

    printf("load array : ");
    loadArray.printYourself();
    printf("\n");
}


void
Node :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    DofManager :: updateYourself(tStep);

    fMode mode = domain->giveEngngModel()->giveFormulation();

    if ( mode == AL ) { // updated Lagrange
        for ( Dof *d: *this ) {
            DofIDItem id = d->giveDofID();
            if ( id == D_u || id == D_v || id == D_w ) {
                int ic = id - D_u + 1;
                coordinates.at(ic) += d->giveUnknown(VM_Incremental, tStep);
            } else if ( id == V_u || id == V_v || id == V_w ) {
                int ic = id - V_u + 1;
                coordinates.at(ic) += d->giveUnknown(VM_Total, tStep) * tStep->giveTimeIncrement();
            }
        }
    }
}


double
Node :: giveUpdatedCoordinate(int ic, TimeStep *tStep, double scale)
//
// returns coordinate + scale * displacement
//
{
#ifdef DEBUG
    if ( ( ic < 1 ) || ( ic > 3 ) ) {
        OOFEM_ERROR("Can't return non-existing coordinate (index not in range 1..3)");
        return 0.;
    }
#endif

    if ( tStep->isTheCurrentTimeStep() ) {
        double coordinate = this->giveCoordinate(ic);
        if ( !this->hasLocalCS() ) {
            // this has no local cs.
            for ( Dof *d: *this ) {
                DofIDItem id = d->giveDofID();
                if ( id == ic ) {
                    coordinate += scale * d->giveUnknown(VM_Total, tStep);
                    break;
                } else if ( id - V_u + 1 == ic ) {
                    coordinate += scale * d->giveUnknown(VM_Total, tStep) * tStep->giveTimeIncrement();
                    break;
                }
            }
        } else {
            //
            // this has local cs.
            // We must perform transformation of displacements DOFs
            // in to global c.s and then to add them to global coordinates.
            //
            FloatMatrix *T = this->giveLocalCoordinateTriplet();
            FloatArray displacements( T->giveNumberOfRows() );
            displacements.zero();

            for ( Dof *d: *this ) {
                DofIDItem id = d->giveDofID();
                if ( id == D_u || id == D_v || id == D_w ) {
                    int ic2 = id - D_u + 1;
                    displacements.at(ic2) = scale * d->giveUnknown(VM_Total, tStep);
                } else if ( id == V_u || id == V_v || id == V_w ) {
                    int ic2 = id - V_u + 1;
                    displacements.at(ic2) = scale * d->giveUnknown(VM_Total, tStep) * tStep->giveTimeIncrement();
                }
            }

            // perform transformation for desired displacement
            for ( int i = 1; i <= 3; i++ ) {
                coordinate += displacements.at(i) * T->at(i, ic);
            }
        }

        return coordinate;
    } else {
        OOFEM_ERROR("Can't return updatedCoordinate for non-current timestep");
    }

    return 0.;
}


void
Node :: giveUpdatedCoordinates(FloatArray &coord, TimeStep *tStep, double scale)
//
// returns coordinate + scale * displacement
//
{
    if ( tStep->isTheCurrentTimeStep() ) {
        FloatArray vec;
        coord = this->coordinates;
        this->giveUnknownVectorOfType(vec, DisplacementVector, VM_Total, tStep);
        for ( int i = 1; i <= coord.giveSize(); i++ ) {
            coord.at(i) += scale * vec.at(i);
        }
    }
}


int
Node :: checkConsistency()
{
    /*
     * Checks internal data consistency in node.
     * Current implementation checks (when receiver has slave dofs) if receiver has the same
     * coordinate system as master dofManager of slave dof.
     */
    int result;
    int nslaves = 0;

    result = DofManager :: checkConsistency();

    for ( Dof *dof: *this ) {
        if ( dynamic_cast< SimpleSlaveDof * >( dof ) ) {
            nslaves++;
        }
    }

    if ( nslaves == 0 ) {
        return result;             // return o.k. if no slaves exists
    }

    IntArray masterDofManagers(nslaves);
    int numberOfMDM = 0; // counter of different master dofManagers
    int master, alreadyFound = 0;
    Node *masterNode;

    for ( Dof *dof: *this ) {
        SimpleSlaveDof *sdof = dynamic_cast< SimpleSlaveDof * >( dof );
        if ( sdof ) {
            alreadyFound  = 0;
            master = sdof->giveMasterDofManagerNum();
            for ( int j = 1; j <= numberOfMDM; j++ ) {
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
                    OOFEM_WARNING("master dofManager is not compatible", 1);
                    result = 0;
                } else if ( !this->hasSameLCS(masterNode) ) {
                    OOFEM_WARNING("different lcs for master/slave nodes", 1);
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

    if ( ( this->hasLocalCS() ) && ( remote->hasLocalCS() ) ) {
        for ( int k = 1; k <= 3; k++ ) {
            for ( int l = 1; l <= 3; l++ ) {
                if ( fabs( thisLcs->at(k, l) - masterLcs->at(k, l) ) > 1.e-4 ) {
                    return false;
                }
            }
        }
    } else if ( this->hasLocalCS() ) {
        for ( int k = 1; k <= 3; k++ ) {
            for ( int l = 1; l <= 3; l++ ) {
                if ( fabs( thisLcs->at(k, l) - ( k == l ) ) > 1.e-4 ) {
                    return false;
                }
            }
        }
    } else if ( remote->hasLocalCS() ) {
        for ( int k = 1; k <= 3; k++ ) {
            for ( int l = 1; l <= 3; l++ ) {
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
    DofIDItem id;

    if ( localCoordinateSystem == NULL ) {
        answer.clear();
        return false;
    } else {
        ///@todo This relies on the order of the dofs, not good.. / Mikael
        if ( dofIDArry.isEmpty() ) {
            // response for all local dofs is computed

            int numberOfDofs = this->giveNumberOfDofs();
            answer.resize(numberOfDofs, numberOfDofs);
            answer.zero();

            int i = 0;
            for ( Dof *dof: *this ) {
                // test for vector quantities
                i++;
                int j = 0;
                switch ( id = dof->giveDofID() ) {
                case D_u:
                case D_v:
                case D_w:
                    for ( Dof *dof2: *this ) {
                        DofIDItem id2 = dof2->giveDofID();
                        j++;
                        if ( ( id2 == D_u ) || ( id2 == D_v ) || ( id2 == D_w ) ) {
                            answer.at(j, i) = localCoordinateSystem->at( ( int ) ( id ) - ( int ) ( D_u ) + 1,
                                                                        ( int ) ( id2 ) - ( int ) ( D_u ) + 1 );
                        }
                    }

                    break;

                case V_u:
                case V_v:
                case V_w:
                    for ( Dof *dof2: *this ) {
                        DofIDItem id2 = dof2->giveDofID();
                        j++;
                        if ( ( id2 == V_u ) || ( id2 == V_v ) || ( id2 == V_w ) ) {
                            answer.at(j, i) = localCoordinateSystem->at( ( int ) ( id ) - ( int ) ( V_u ) + 1,
                                                                        ( int ) ( id2 ) - ( int ) ( V_u ) + 1 );
                        }
                    }

                    break;

                case R_u:
                case R_v:
                case R_w:
                    for ( Dof *dof2: *this ) {
                        DofIDItem id2 = dof2->giveDofID();
                        j++;
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
                    OOFEM_ERROR("unknown dofID (%s)", __DofIDItemToString(id).c_str());
                }
            }
        } else { // end if (dofIDArry.isEmpty())
            // map is provided -> assemble for requested dofs
            int size = dofIDArry.giveSize();
            answer.resize(size, size);
            answer.zero();

            for ( int i = 1; i <= size; i++ ) {
                // test for vector quantities
                switch ( id = ( DofIDItem ) dofIDArry.at(i) ) {
                case D_u:
                case D_v:
                case D_w:
                    for ( int j = 1; j <= size; j++ ) {
                        DofIDItem id2 = ( DofIDItem ) dofIDArry.at(j);
                        if ( ( id2 == D_u ) || ( id2 == D_v ) || ( id2 == D_w ) ) {
                            answer.at(j, i) = localCoordinateSystem->at( ( int ) ( id ) - ( int ) ( D_u ) + 1, ( int ) ( id2 ) - ( int ) ( D_u ) + 1 );
                        }
                    }

                    break;

                case V_u:
                case V_v:
                case V_w:
                    for ( int j = 1; j <= size; j++ ) {
                        DofIDItem id2 = ( DofIDItem ) dofIDArry.at(j);
                        if ( ( id2 == V_u ) || ( id2 == V_v ) || ( id2 == V_w ) ) {
                            answer.at(j, i) = localCoordinateSystem->at( ( int ) ( id ) - ( int ) ( V_u ) + 1, ( int ) ( id2 ) - ( int ) ( V_u ) + 1 );
                        }
                    }

                    break;

                case R_u:
                case R_v:
                case R_w:
                    for ( int j = 1; j <= size; j++ ) {
                        DofIDItem id2 = ( DofIDItem ) dofIDArry.at(j);
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
                    OOFEM_ERROR("unknown dofID (%s)", __DofIDItemToString(id).c_str());
                }
            }
        } // end map is provided -> assemble for requested dofs
    } // end localCoordinateSystem defined
    return true;
}


contextIOResultType
Node :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = DofManager :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        int _haslcs = hasLocalCS();
        if ( ( iores = coordinates.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream.write(_haslcs) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( _haslcs ) {
            if ( ( iores = localCoordinateSystem->storeYourself(stream) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        }
    }

    return CIO_OK;
}


contextIOResultType
Node :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = DofManager :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        int _haslcs;
        if ( ( iores = coordinates.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream.read(_haslcs) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( _haslcs ) {
            if ( localCoordinateSystem == NULL ) {
                localCoordinateSystem = new FloatMatrix();
            }

            if ( ( iores = localCoordinateSystem->restoreYourself(stream) ) != CIO_OK ) {
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
Node :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
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
        if ( this->giveDomain()->hasXfemManager() ) {
            XfemManager *xf = this->giveDomain()->giveXfemManager();
            for ( int i = 1; i <= xf->giveNumberOfEnrichmentItems(); i++ ) {
                if ( xf->giveEnrichmentItem(i)->isDofManEnriched(* this) ) {
                    EASValsSetMType(SQUARE_MARKER);
                }
            }
        }

 #endif

        if ( this->giveParallelMode() == DofManager_local ) {
            EASValsSetColor( gc.getNodeColor() );
        } else if ( this->giveParallelMode() == DofManager_shared ) {
            EASValsSetColor( gc.getDeformedElementColor() );
        } else {
            EASValsSetColor( gc.getCrackPatternColor() );
        }

        bool ordinary = true;

        for ( Dof *dof: *this ) {
            if ( dof->isPrimaryDof() ) {
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
 
        sprintf( num, "%d(%d)", this->giveNumber(), this->giveGlobalNumber() );
        go = CreateAnnText3D(p, num);
        EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    } else if ( mode == OGC_essentialBC ) {
        int i, hasDisplSupport [ 3 ], hasRotSupport [ 3 ], hasAny = 0;
        if ( !tStep ) {
            TimeStep __temp( domain->giveEngngModel() );
            tStep = & __temp;
        }

        WCRec pp [ 2 ];

        for ( i = 0; i < 3; i++ ) {
            hasDisplSupport [ i ] = 0;
            hasRotSupport [ i ] = 0;
        }

        for ( Dof *dof: *this ) {
            if ( dof->hasBc(tStep) ) {
                hasAny = 1;
                switch ( dof->giveDofID() ) {
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
            pp [ 1 ].x = pp [ 1 ].y = pp [ 1 ].z = 0.0;

            FloatArray load, f;
            FloatMatrix t;
            IntArray dofIDArry(0);
            
            load.clear();
            for ( int iload : *this->giveLoadArray() ) {   // to more than one load
                Load *loadN = domain->giveLoad(iload);
                this->computeLoadVector(f, loadN, ExternalForcesVector, tStep, VM_Total);
                load.add(f);
            }
            if ( computeL2GTransformation(t, dofIDArry) ) {
                load.rotatedWith(t, 'n');
            }

            FloatArray force(3), momentum(3);
            int i = 0;
            for ( Dof *dof: *this ) {
                i++;
                switch ( dof->giveDofID() ) {
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

        if ( gc.giveIntVarType() == IST_Velocity ) {
            WCRec p [ 2 ]; /* point */
            double defScale = gc.getDefScale();

            p [ 0 ].x = p [ 1 ].x = ( FPNum ) this->giveCoordinate(1);
            p [ 0 ].y = p [ 1 ].y = ( FPNum ) this->giveCoordinate(2);
            p [ 0 ].z = p [ 1 ].z = ( FPNum ) this->giveCoordinate(3);

            //p[1].x = p[1].y = p[1].z = 0.0;
            for ( Dof *dof: *this ) {
                if ( dof->giveDofID() == V_u ) {
                    p [ 1 ].x = defScale * dof->giveUnknown(VM_Total, tStep);
                } else if ( dof->giveDofID() == V_v ) {
                    p [ 1 ].y = defScale * dof->giveUnknown(VM_Total, tStep);
                } else if ( dof->giveDofID() == V_w ) {
                    p [ 1 ].z = defScale * dof->giveUnknown(VM_Total, tStep);
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
