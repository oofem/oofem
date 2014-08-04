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

#include "refinedelement.h"
#include "element.h"
#include "node.h"
#include "dof.h"
#include "mathfem.h"

#include <cstdarg>
#include <cstdlib> // For abort
///@todo Replace abort() with OOFEM_ERROR

namespace oofem {
RefinedElement :: RefinedElement(Domain *d, int elem, int level) : fineNodeList(0), boundaryFlag()
    // Constructor
{
    Element *element;
    int nodes, sides, dim, len;

    this->elementId = elem;

    element = d->giveElement(elem);
    nodes = element->giveNumberOfDofManagers();
    sides = element->giveNumberOfBoundarySides();
    dim = element->giveSpatialDimension();

    len = 1;
    while ( dim-- ) {
        len *= ( level + 2 );
    }

    fineNodeList.resize(nodes);
    for ( int inode = 0; inode < nodes; inode++ ) {
        fineNodeList[inode].resize(len);
    }

    this->boundaryFlag.resize(sides);
    // this -> boundaryFlag.zero();
}



RefinedElement :: ~RefinedElement()
// Destructor.
{ }


IntArray *
RefinedElement :: giveFineNodeArray(int node)
{
    return & this->fineNodeList[node-1];
}


void
RefinedElement :: giveBoundaryFlagArray(int inode, Element *element, IntArray &answer)
{
    /* note: number of connected edges must correspond to edge ordering in face_ed_nd and quad_ed_nd (see refinedmesh.C);
     * ordering of edges at a particular node is given by fine node ordering {m = 0, n = 0};
     * 1-based indexing is used contrary to 0-based indexing in refinedmesh.C */

    static int face_con_ed [ 3 ] [ 2 ] = { { 3, 1 }, { 1, 2 }, { 2, 3 } };
    static int quad_con_ed [ 4 ] [ 2 ] = { { 4, 1 }, { 1, 2 }, { 2, 3 }, { 3, 4 } };

    /* note: number of connected faces must correspond to face ordering in tetra_fc_nd and hexa_fc_nd (see refinedmesh.C);
     * ordering of faces at a particular node is given by fine node ordering {m = 0, n = 0, k = 0};
     * 1-based indexing is used contrary to 0-based indexing in refinedmesh.C */

    static int tetra_con_fc [ 4 ] [ 3 ] = { { 4, 2, 1 }, { 2, 3, 1 }, { 3, 4, 1 }, { 3, 2, 4 } };
    static int hexa_con_fc [ 8 ] [ 3 ] = { { 5, 2, 1 }, { 2, 3, 1 }, { 3, 4, 1 }, { 4, 5, 1 }, { 2, 5, 6 }, { 3, 2, 6 }, { 4, 3, 6 }, { 5, 4, 6 } };

    int i, dim, *con = NULL;

#ifdef DEBUG
    if ( inode > element->giveNumberOfNodes() ) {
        abort();
    }

#endif

    answer.resize( dim = element->giveSpatialDimension() );

    switch ( element->giveGeometryType() ) {
    case EGT_line_1:
    case EGT_line_2:
        answer.at(1) = boundaryFlag.at(inode);
        return;

        break;
    case EGT_triangle_1:
    case EGT_triangle_2:
        con = face_con_ed [ inode - 1 ];
        break;
    case EGT_quad_1:
        // case EGT_quad_2:
        con = quad_con_ed [ inode - 1 ];
        break;
    case EGT_tetra_1:
        //  case EGT_tetra_2:
        con = tetra_con_fc [ inode - 1 ];
        break;
    case EGT_hexa_1:
        // case EGT_hexa_2:
        con = hexa_con_fc [ inode - 1 ];
        break;
    default:
        OOFEM_ERROR("Unsupported geometry type");
    }

    for ( i = 0; i < dim; i++ ) {
        answer.at(i + 1) = boundaryFlag.at(con [ i ]);
    }
}



bool
RefinedElement :: giveBcDofArray1D(int inode, Element *element, IntArray &sideBcDofId, int &sideNumBc, TimeStep *tStep)
{
    static int edge_con_nd [ 2 ] = {
        2, 1
    };

    int nodeNumBc;
    Node *node;
    IntArray nodeBcDofId;

    node = element->giveNode(inode);

    nodeNumBc = 0;
    nodeBcDofId.resize( node->giveNumberOfDofs() );
    for ( Dof *dof: *node ) {
        if ( dof->hasBc(tStep) != 0 ) {
            nodeBcDofId.at(++nodeNumBc) = dof->giveDofID();
        }
    }

    if ( nodeNumBc == 0 ) {
        return false;
    }

    sideNumBc = this->giveCompatibleBcDofArray(element->giveNode(edge_con_nd [ inode - 1 ]), node, nodeBcDofId, nodeNumBc,
                                               sideBcDofId, VM_Total, tStep);
    return true;
}


bool
RefinedElement :: giveBcDofArray2D(int inode, Element *element, std::vector< IntArray > &sideBcDofIdList, IntArray &sideNumBc, TimeStep *tStep)
{
    /* note: ordering of connected nodes is given by fine node ordering {m = 0, n = 0};
     * 1-based indexing is used contrary to 0-based indexing in refinedmesh.C */

    static int face_con_nd [ 3 ] [ 2 ] = { { 3, 2 }, { 1, 3 }, { 2, 1 } };
    static int quad_con_nd [ 4 ] [ 2 ] = { { 4, 2 }, { 1, 3 }, { 2, 4 }, { 3, 1 } };

    int *con = NULL, iside, nodeNumBc;
    Node *node;
    IntArray nodeBcDofId;

    node = element->giveNode(inode);

    nodeNumBc = 0;
    nodeBcDofId.resize( node->giveNumberOfDofs() );
    for ( Dof *dof: *node ) {
        if ( dof->hasBc(tStep) != 0 ) {
            nodeBcDofId.at(++nodeNumBc) = dof->giveDofID();
        }
    }

    if ( nodeNumBc == 0 ) {
        return false;
    }

    switch ( element->giveGeometryType() ) {
    case EGT_triangle_1:
    case EGT_triangle_2:
        con = face_con_nd [ inode - 1 ];
        break;
    case EGT_quad_1:
        // case EGT_quad_2:
        con = quad_con_nd [ inode - 1 ];
        break;
    default:
        OOFEM_ERROR("Unsupported geometry type");
    }

    for ( iside = 0; iside < 2; iside++ ) {
        sideNumBc.at(iside + 1) = this->giveCompatibleBcDofArray(element->giveNode(con [ iside ]),
                                                                 node, nodeBcDofId, nodeNumBc,
                                                                 sideBcDofIdList[iside],
                                                                 VM_Total, tStep);
    }

    return true;
}


bool
RefinedElement :: giveBcDofArray3D(int inode, Element *element, std::vector< IntArray > &sideBcDofIdList, IntArray &sideNumBc,
                                   std::vector< IntArray > &faceBcDofIdList, IntArray &faceNumBc, TimeStep *tStep)
{
    /* note: ordering of connected nodes is given by fine node ordering {n = k = 0, m = k = 0, m = n = 0};
     * 1-based indexing is used contrary to 0-based indexing in refinedmesh.C */

    static int tetra_con_nd [ 4 ] [ 3 ] = { { 2, 3, 4 }, { 3, 1, 4 }, { 1, 2, 4 }, { 1, 3, 2 } };
    static int hexa_con_nd [ 8 ] [ 3 ] = { { 2, 4, 5 }, { 3, 1, 6 }, { 4, 2, 7 }, { 1, 3, 8 }, { 8, 6, 1 }, { 5, 7, 2 }, { 6, 8, 3 }, { 7, 5, 4 } };

    /* note: number of connected faces must correspond to face ordering in tetra_fc_nd and hexa_fc_nd;
     * ordering of faces at a particular node is given by fine node ordering {m = 0, n = 0, k = 0};
     * 1-based indexing is used contrary to 0-based indexing in refinedmesh.C */

    static int tetra_con_fc [ 4 ] [ 3 ] = { { 4, 2, 1 }, { 2, 3, 1 }, { 3, 4, 1 }, { 3, 2, 4 } };
    static int hexa_con_fc [ 8 ] [ 3 ] = { { 5, 2, 1 }, { 2, 3, 1 }, { 3, 4, 1 }, { 4, 5, 1 }, { 2, 5, 6 }, { 3, 2, 6 }, { 4, 3, 6 }, { 5, 4, 6 } };

    /* note: ordering of nodes on faces is irrelevant;
     * 1-based indexing is used contrary to 0-based indexing in refinedmesh.C */

    static int tetra_fc_nd [ 4 ] [ 3 ] = { { 1, 2, 3 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 1, 4 } };
    static int hexa_fc_nd [ 6 ] [ 4 ] = { { 1, 2, 3, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 4, 8, 7 }, { 4, 1, 5, 8 }, { 8, 7, 6, 5 } };

    int *con = NULL, iside, iface, jnode, nodeNumBc, fcNumBc;
    Node *node;
    IntArray nodeBcDofId, faceBcDofId;
    bool hasBc = false;

    node = element->giveNode(inode);

    nodeNumBc = 0;
    nodeBcDofId.resize( node->giveNumberOfDofs() );
    for ( Dof *dof: *node ) {
        if ( dof->hasBc(tStep) != 0 ) {
            nodeBcDofId.at(++nodeNumBc) = dof->giveDofID();
        }
    }

    if ( nodeNumBc == 0 ) {
        return false;
    }

    switch ( element->giveGeometryType() ) {
    case EGT_tetra_1:
        // case EGT_tetra_2:
        con = tetra_con_nd [ inode - 1 ];
        break;
    case EGT_hexa_1:
        // case EGT_hexa_2:
        con = hexa_con_nd [ inode - 1 ];
        break;
    default:
        OOFEM_ERROR("Unsupported geometry type");
    }

    for ( iside = 0; iside < 3; iside++ ) {
        sideNumBc.at(iside + 1) = this->giveCompatibleBcDofArray(element->giveNode(con [ iside ]),
                                                                 node, nodeBcDofId, nodeNumBc,
                                                                 sideBcDofIdList[iside],
                                                                 VM_Total, tStep);
        if ( sideNumBc.at(iside + 1) != 0 ) {
            hasBc = true;
        }
    }

    if ( hasBc == true ) {
        faceBcDofId.resize(nodeNumBc);

        switch ( element->giveGeometryType() ) {
        case EGT_tetra_1:
            //  case EGT_tetra_2:

            for ( int i = 0; i < 3; i++ ) {        // there are always 3 faces at a node
                iface = tetra_con_fc [ inode - 1 ] [ i ];
                fcNumBc = nodeNumBc;
                for ( int idof = 1; idof <= fcNumBc; idof++ ) {
                    faceBcDofId.at(idof) = nodeBcDofId.at(idof);
                }

                for ( int j = 0; j < 3; j++ ) {    // there 3 nodes per face
                    jnode = tetra_fc_nd [ iface - 1 ] [ j ];
                    if ( jnode == inode ) {
                        continue;
                    }

                    fcNumBc = this->giveCompatibleBcDofArray(element->giveNode(jnode), node, faceBcDofId, fcNumBc,
                                                             faceBcDofIdList[i],
                                                             VM_Total, tStep);
                    for ( int idof = 1; idof <= fcNumBc; idof++ ) {
                        faceBcDofId.at(idof) = faceBcDofIdList[i].at(idof);
                    }
                }

                faceNumBc.at(i + 1) = fcNumBc;
            }

            break;
        case EGT_hexa_1:
            //  case EGT_hexa_2:

            for ( int i = 0; i < 3; i++ ) {        // there are always 3 faces at a node
                iface = hexa_con_fc [ inode - 1 ] [ i ];
                fcNumBc = nodeNumBc;
                for ( int idof = 1; idof <= fcNumBc; idof++ ) {
                    faceBcDofId.at(idof) = nodeBcDofId.at(idof);
                }

                for ( int j = 0; j < 4; j++ ) {    // there 4 nodes per face
                    jnode = hexa_fc_nd [ iface - 1 ] [ j ];
                    if ( jnode == inode ) {
                        continue;
                    }

                    fcNumBc = this->giveCompatibleBcDofArray(element->giveNode(jnode), node, faceBcDofId, fcNumBc,
                                                             faceBcDofIdList[i],
                                                             VM_Total, tStep);
                    for ( int idof = 1; idof <= fcNumBc; idof++ ) {
                        faceBcDofId.at(idof) = faceBcDofIdList[i].at(idof);
                    }
                }

                faceNumBc.at(i + 1) = fcNumBc;
            }

            break;
        default:
            OOFEM_ERROR("Unsupported geometry type");
        }
    }

    return true;
}



bool
RefinedElement :: giveBoundaryLoadArray1D(int inode, Element *element, IntArray &boundaryLoadArray)
{
    int iload, loads, bloads;
    IntArray *loadArray;

    if ( ( loads = ( loadArray = element->giveBoundaryLoadArray() )->giveSize() ) == 0 ) {
        return false;
    }

#ifdef DEBUG
    if ( loads / 2 * 2 != loads ) {
        abort();
    }

#endif

    boundaryLoadArray.resize(loads);
    bloads = 0;

    for ( iload = 1; iload <= loads; iload += 2 ) {
        if ( loadArray->at(iload + 1) != 1 ) {
            continue;
        }

        bloads += 2;
        boundaryLoadArray.at(bloads - 1) = loadArray->at(iload);
        boundaryLoadArray.at(bloads) = 1;
    }

    boundaryLoadArray.resize(bloads);

    return true;
}



bool
RefinedElement :: giveBoundaryLoadArray2D(int inode, Element *element, std::vector< IntArray > &boundaryLoadList)
{
    /* note: number of connected edges must correspond to OOFEM element side numbering;
     * ordering of edges at a particular node is given by fine node ordering {m = 0, n = 0};
     * 1-based indexing is used */

    static int face_con_ed [ 3 ] [ 2 ] = { { 3, 1 }, { 1, 2 }, { 2, 3 } };
    static int quad_con_ed [ 4 ] [ 2 ] = { { 4, 1 }, { 1, 2 }, { 2, 3 }, { 3, 4 } };

    static int fine_quad_side [ 2 ] = {
        4, 1
    };                                                           // {m = 0, n = 0}

    int iside, iload, loads, bloads, side, *con = NULL;
    IntArray *loadArray;

    if ( ( loads = ( loadArray = element->giveBoundaryLoadArray() )->giveSize() ) == 0 ) {
        return false;
    }

#ifdef DEBUG
    if ( loads / 2 * 2 != loads ) {
        abort();
    }

#endif

    switch ( element->giveGeometryType() ) {
    case EGT_triangle_1:
    case EGT_triangle_2:
        con = face_con_ed [ inode - 1 ];
        break;
    case EGT_quad_1:
        // case EGT_quad_2:
        con = quad_con_ed [ inode - 1 ];
        break;
    default:
        OOFEM_ERROR("Unsupported geometry type");
    }

    for ( iside = 0; iside < 2; iside++ ) {
        IntArray &boundaryLoadArray = boundaryLoadList[iside];
        boundaryLoadArray.resize(loads);
        bloads = 0;

        side = con [ iside ];
        for ( iload = 1; iload <= loads; iload += 2 ) {
            if ( loadArray->at(iload + 1) != side ) {
                continue;
            }

            bloads += 2;
            boundaryLoadArray.at(bloads - 1) = loadArray->at(iload);
            boundaryLoadArray.at(bloads) = fine_quad_side [ iside ];
        }

        boundaryLoadArray.resize(bloads);
    }

    return true;
}



bool
RefinedElement :: giveBoundaryLoadArray3D(int inode, Element *element, std::vector< IntArray > &boundaryLoadList)
{
    /* note: number of connected faces must correspond to OOFEM element side numbering;
     * ordering of faces at a particular node is given by fine node ordering {m = 0, n = 0, k = 0};
     * 1-based indexing is used */

    static int tetra_con_fc [ 4 ] [ 3 ] = { { 4, 2, 1 }, { 2, 3, 1 }, { 3, 4, 1 }, { 3, 2, 4 } };
    static int hexa_con_fc [ 8 ] [ 3 ] = { { 6, 3, 1 }, { 3, 4, 1 }, { 4, 5, 1 }, { 5, 6, 1 }, { 3, 6, 2 }, { 4, 3, 2 }, { 5, 4, 2 }, { 6, 5, 2 } };

    static int fine_hexa_side [ 3 ] = {
        6, 3, 1
    };                                                           // {m = 0, n = 0, k = 0}

    int iside, iload, loads, bloads, side, *con = NULL;
    IntArray *loadArray;

    if ( ( loads = ( loadArray = element->giveBoundaryLoadArray() )->giveSize() ) == 0 ) {
        return false;
    }

#ifdef DEBUG
    if ( loads / 2 * 2 != loads ) {
        abort();
    }

#endif

    switch ( element->giveGeometryType() ) {
    case EGT_tetra_1:
        //  case EGT_tetra_2:
        con = tetra_con_fc [ inode - 1 ];
        break;
    case EGT_hexa_1:
        // case EGT_hexa_2:
        con = hexa_con_fc [ inode - 1 ];
        break;
    default:
        OOFEM_ERROR("Unsupported geometry type");
    }

    for ( iside = 0; iside < 3; iside++ ) {
        IntArray &boundaryLoadArray = boundaryLoadList[iside];
        boundaryLoadArray.resize(loads);
        bloads = 0;

        side = con [ iside ];
        for ( iload = 1; iload <= loads; iload += 2 ) {
            if ( loadArray->at(iload + 1) != side ) {
                continue;
            }

            bloads += 2;
            boundaryLoadArray.at(bloads - 1) = loadArray->at(iload);
            boundaryLoadArray.at(bloads) = fine_hexa_side [ iside ];
        }

        boundaryLoadArray.resize(bloads);
    }

    return true;
}

int
RefinedElement :: giveCompatibleBcDofArray(Node *master_node, Node *slave_node, IntArray &dofIDArray, int dofs,
                                           IntArray &answer, ValueModeType mode, TimeStep *tStep)
{
    FloatMatrix *Lcs, *nodeLcs, trFromNodeLcsToLcs;
    bool compatibleCS, newLcs, newNodeLcs;
    double epsilon = 1.0e-9;
    int compDofs = 0;
    DofIDItem dofId;
    double bcValue;
    int bcId;

    newLcs = newNodeLcs = false;
    compatibleCS = true;

    Lcs = master_node->giveLocalCoordinateTriplet();
    nodeLcs = slave_node->giveLocalCoordinateTriplet();
    if ( Lcs != NULL || nodeLcs != NULL ) {
        // check the compatibility of CSs

        if ( Lcs == NULL ) {
            Lcs = new FloatMatrix(3, 3);
            Lcs->zero();
            Lcs->at(1, 1) = Lcs->at(2, 2) = Lcs->at(3, 3) = 1.0;
            newLcs = true;
        }

        if ( nodeLcs == NULL ) {
            nodeLcs = new FloatMatrix(3, 3);
            nodeLcs->zero();
            nodeLcs->at(1, 1) = nodeLcs->at(2, 2) = nodeLcs->at(3, 3) = 1.0;
            newNodeLcs = true;
        }

        trFromNodeLcsToLcs.beProductTOf(* Lcs, * nodeLcs);

        if ( newLcs == true ) {
            delete Lcs;
        }

        if ( newNodeLcs == true ) {
            delete nodeLcs;
        }

        for ( int i = 1; i <= 3; i++ ) {
            if ( fabs(trFromNodeLcsToLcs.at(i, i) - 1.0) > epsilon ) {
                compatibleCS = false;
                break;
            }
        }
    }

    answer.resize(dofs);

    if ( compatibleCS == true ) {
        for ( int i = 1; i <= dofs; i++ ) {
            Dof *nodeDof = slave_node->giveDofWithID( dofIDArray.at(i) );

#ifdef DEBUG
            if ( nodeDof->hasBc(tStep) == false ) {
                OOFEM_ERROR("dof has no BC");
            }

#endif

            dofId = nodeDof->giveDofID();
            bcId = nodeDof->giveBcId();
            bcValue = nodeDof->giveBcValue(mode, tStep);

            for ( Dof *dof: *master_node ) {
                if ( dof->hasBc(tStep) == false ) {
                    continue;
                }

                if ( dof->giveDofID() != dofId ) {
                    continue;
                }

                if ( dof->giveBcId() == bcId ) {
                    answer.at(++compDofs) = dofIDArray.at(i);
                    break;
                }

                if ( dof->giveBcValue(mode, tStep) == bcValue ) {
                    answer.at(++compDofs) = dofIDArray.at(i);
                    break;
                }
            }
        }
    }

    /*
     * else{
     * if(orthogonal == true){
     * }
     * else{
     * }
     * }
     */

    answer.resize(compDofs);

    return ( compDofs );
}


std :: string RefinedElement :: errorInfo(const char *func) const
{
    return std :: string("RefinedElement::") + func + ", number: " + std::to_string(this->elementId);
}
} // end namespace oofem
