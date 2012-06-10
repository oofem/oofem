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

#include "refinedelement.h"
#include "element.h"
#include "node.h"
#include "dof.h"
#include "mathfem.h"
#include "oofem_limits.h"

#ifndef __MAKEDEPEND
 #include <stdarg.h>
#endif

namespace oofem {
RefinedElement :: RefinedElement(Domain *d, int elem, int level) : fineNodeList(0), boundaryFlag()
    // Constructor
{
    Element *element;
    int inode, nodes, sides, dim, len;
    IntArray *connectivity;

    this->elementId = elem;

    element = d->giveElement(elem);
    nodes = element->giveNumberOfDofManagers();
    sides = element->giveNumberOfBoundarySides();
    dim = element->giveSpatialDimension();

    len = 1;
    while ( dim-- ) {
        len *= ( level + 2 );
    }

    fineNodeList.growTo(nodes);
    for ( inode = 1; inode <= nodes; inode++ ) {
        connectivity = new IntArray(len);
        fineNodeList.put(inode, connectivity);
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
    if ( this->fineNodeList.includes(node) ) {
        return this->fineNodeList.at(node);
    }

    /*
     * else {
     *   _errori ("giveNodeAssocFineNodeList: No such node list defined: ", node);
     * }
     */
    return NULL;
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
        _error("giveBoundaryFlagArray: Unsupported geometry type");
    }

    for ( i = 0; i < dim; i++ ) {
        answer.at(i + 1) = boundaryFlag.at(con [ i ]);
    }
}



bool
RefinedElement :: giveBcDofArray1D(int inode, Element *element, IntArray *sideBcDofId, int &sideNumBc, TimeStep *tStep)
{
    static int edge_con_nd [ 2 ] = {
        2, 1
    };

    int idof, nodeNumBc;
    Node *node;
    IntArray nodeBcDofId;

    node = element->giveNode(inode);

    nodeNumBc = 0;
    nodeBcDofId.resize( node->giveNumberOfDofs() );
    for ( idof = 1; idof <= node->giveNumberOfDofs(); idof++ ) {
        if ( node->giveDof(idof)->hasBc(tStep) != 0 ) {
            nodeBcDofId.at(++nodeNumBc) = idof;
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
RefinedElement :: giveBcDofArray2D(int inode, Element *element, AList< IntArray > &sideBcDofIdList, IntArray &sideNumBc, TimeStep *tStep)
{
    /* note: ordering of connected nodes is given by fine node ordering {m = 0, n = 0};
     * 1-based indexing is used contrary to 0-based indexing in refinedmesh.C */

    static int face_con_nd [ 3 ] [ 2 ] = { { 3, 2 }, { 1, 3 }, { 2, 1 } };
    static int quad_con_nd [ 4 ] [ 2 ] = { { 4, 2 }, { 1, 3 }, { 2, 4 }, { 3, 1 } };

    int idof, *con = NULL, iside, nodeNumBc;
    Node *node;
    IntArray nodeBcDofId;

    node = element->giveNode(inode);

    nodeNumBc = 0;
    nodeBcDofId.resize( node->giveNumberOfDofs() );
    for ( idof = 1; idof <= node->giveNumberOfDofs(); idof++ ) {
        if ( node->giveDof(idof)->hasBc(tStep) != 0 ) {
            nodeBcDofId.at(++nodeNumBc) = idof;
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
        _error("giveBcDofArray2D: Unsupported geometry type");
    }

    for ( iside = 0; iside < 2; iside++ ) {
        sideNumBc.at(iside + 1) = this->giveCompatibleBcDofArray(element->giveNode(con [ iside ]),
                                                                 node, nodeBcDofId, nodeNumBc,
                                                                 sideBcDofIdList.at(iside + 1),
                                                                 VM_Total, tStep);
    }

    return true;
}


bool
RefinedElement :: giveBcDofArray3D(int inode, Element *element, AList< IntArray > &sideBcDofIdList, IntArray &sideNumBc,
                                   AList< IntArray > &faceBcDofIdList, IntArray &faceNumBc, TimeStep *tStep)
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

    int i, j, idof, *con = NULL, iside, iface, jnode, nodeNumBc, fcNumBc;
    Node *node;
    IntArray nodeBcDofId, faceBcDofId;
    bool hasBc = false;

    node = element->giveNode(inode);

    nodeNumBc = 0;
    nodeBcDofId.resize( node->giveNumberOfDofs() );
    for ( idof = 1; idof <= node->giveNumberOfDofs(); idof++ ) {
        if ( node->giveDof(idof)->hasBc(tStep) != 0 ) {
            nodeBcDofId.at(++nodeNumBc) = idof;
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
        _error("giveBcDofArray3D: Unsupported geometry type");
    }

    for ( iside = 0; iside < 3; iside++ ) {
        sideNumBc.at(iside + 1) = this->giveCompatibleBcDofArray(element->giveNode(con [ iside ]),
                                                                 node, nodeBcDofId, nodeNumBc,
                                                                 sideBcDofIdList.at(iside + 1),
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

            for ( i = 0; i < 3; i++ ) {        // there are always 3 faces at a node
                iface = tetra_con_fc [ inode - 1 ] [ i ];
                fcNumBc = nodeNumBc;
                for ( idof = 1; idof <= fcNumBc; idof++ ) {
                    faceBcDofId.at(idof) = nodeBcDofId.at(idof);
                }

                for ( j = 0; j < 3; j++ ) {    // there 3 nodes per face
                    jnode = tetra_fc_nd [ iface - 1 ] [ j ];
                    if ( jnode == inode ) {
                        continue;
                    }

                    fcNumBc = this->giveCompatibleBcDofArray(element->giveNode(jnode), node, faceBcDofId, fcNumBc,
                                                             faceBcDofIdList.at(i + 1),
                                                             VM_Total, tStep);
                    for ( idof = 1; idof <= fcNumBc; idof++ ) {
                        faceBcDofId.at(idof) = faceBcDofIdList.at(i + 1)->at(idof);
                    }
                }

                faceNumBc.at(i + 1) = fcNumBc;
            }

            break;
        case EGT_hexa_1:
            //  case EGT_hexa_2:

            for ( i = 0; i < 3; i++ ) {        // there are always 3 faces at a node
                iface = hexa_con_fc [ inode - 1 ] [ i ];
                fcNumBc = nodeNumBc;
                for ( idof = 1; idof <= fcNumBc; idof++ ) {
                    faceBcDofId.at(idof) = nodeBcDofId.at(idof);
                }

                for ( j = 0; j < 4; j++ ) {    // there 4 nodes per face
                    jnode = hexa_fc_nd [ iface - 1 ] [ j ];
                    if ( jnode == inode ) {
                        continue;
                    }

                    fcNumBc = this->giveCompatibleBcDofArray(element->giveNode(jnode), node, faceBcDofId, fcNumBc,
                                                             faceBcDofIdList.at(i + 1),
                                                             VM_Total, tStep);
                    for ( idof = 1; idof <= fcNumBc; idof++ ) {
                        faceBcDofId.at(idof) = faceBcDofIdList.at(i + 1)->at(idof);
                    }
                }

                faceNumBc.at(i + 1) = fcNumBc;
            }

            break;
        default:
            _error("giveBcDofArray3D: Unsupported geometry type");
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
RefinedElement :: giveBoundaryLoadArray2D(int inode, Element *element, AList< IntArray > &boundaryLoadList)
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
    IntArray *loadArray, *boundaryLoadArray;

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
        _error("giveBoundaryLoadArray2D: Unsupported geometry type");
    }

    for ( iside = 0; iside < 2; iside++ ) {
        boundaryLoadArray = boundaryLoadList.at(iside + 1);
        boundaryLoadArray->resize(loads);
        bloads = 0;

        side = con [ iside ];
        for ( iload = 1; iload <= loads; iload += 2 ) {
            if ( loadArray->at(iload + 1) != side ) {
                continue;
            }

            bloads += 2;
            boundaryLoadArray->at(bloads - 1) = loadArray->at(iload);
            boundaryLoadArray->at(bloads) = fine_quad_side [ iside ];
        }

        boundaryLoadArray->resize(bloads);
    }

    return true;
}



bool
RefinedElement :: giveBoundaryLoadArray3D(int inode, Element *element, AList< IntArray > &boundaryLoadList)
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
    IntArray *loadArray, *boundaryLoadArray;

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
        _error("giveBoundaryLoadArray3D: Unsupported geometry type");
    }

    for ( iside = 0; iside < 3; iside++ ) {
        boundaryLoadArray = boundaryLoadList.at(iside + 1);
        boundaryLoadArray->resize(loads);
        bloads = 0;

        side = con [ iside ];
        for ( iload = 1; iload <= loads; iload += 2 ) {
            if ( loadArray->at(iload + 1) != side ) {
                continue;
            }

            bloads += 2;
            boundaryLoadArray->at(bloads - 1) = loadArray->at(iload);
            boundaryLoadArray->at(bloads) = fine_hexa_side [ iside ];
        }

        boundaryLoadArray->resize(bloads);
    }

    return true;
}

int
RefinedElement :: giveCompatibleBcDofArray(Node *master_node, Node *slave_node, IntArray &dofArray, int dofs,
                                           IntArray *answer, ValueModeType mode, TimeStep *tStep)
{
    Dof *dof, *nodeDof;
    FloatMatrix *Lcs, *nodeLcs, trFromNodeLcsToLcs;
    bool compatibleCS = false, orthogonalCS = false, newLcs = false, newNodeLcs = false;
    double epsilon = 1.0e-9;
    int i, j, compDofs = 0;
    DofIDItem dofId;
    double bcValue;
    int bcId;

    compatibleCS = orthogonalCS = true;

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

        for ( i = 1; i <= 3; i++ ) {
            if ( fabs(trFromNodeLcsToLcs.at(i, i) - 1.0) > epsilon ) {
                compatibleCS = false;
                /*
                 *  for(i = 1; i <= 3; i++){
                 *   for(j = 1; j <= 3; j++){
                 *    entry = fabs(trFromNodeLcsToLcs.at(i, j));
                 *    if(entry > epsilon && fabs(1.0 - entry) > epsilon){
                 *     orthogonalCS = false;
                 *     break;
                 *    }
                 *   }
                 *   if(orthogonalCS == false)break;
                 *  }
                 */
                break;
            }
        }
    }

    answer->resize(dofs);

    if ( compatibleCS == true ) {
        for ( i = 1; i <= dofs; i++ ) {
            nodeDof = slave_node->giveDof( dofArray.at(i) );

#ifdef DEBUG
            if ( nodeDof->hasBc(tStep) == false ) {
                _error("extractCompatibleBcDof: dof has no BC");
            }

#endif

            dofId = nodeDof->giveDofID();
            bcId = nodeDof->giveBcId();
            bcValue = nodeDof->giveBcValue(mode, tStep);

            for ( j = 1; j <= master_node->giveNumberOfDofs(); j++ ) {
                dof = master_node->giveDof(j);
                if ( dof->hasBc(tStep) == false ) {
                    continue;
                }

                if ( dof->giveDofID() != dofId ) {
                    continue;
                }

                if ( dof->giveBcId() == bcId ) {
                    answer->at(++compDofs) = dofArray.at(i);
                    break;
                }

                if ( dof->giveBcValue(mode, tStep) == bcValue ) {
                    answer->at(++compDofs) = dofArray.at(i);
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

    answer->resize(compDofs);

    return ( compDofs );
}


/*
 *
 * int element -> giveNode(1) -> giveNumberOfDofs()
 * Dof* element -> giveNode(1) -> giveDof(1)
 * int element -> giveNode(1) -> giveDof(1) -> hasBc(tStep)
 * int element -> giveNode(1) -> giveDof(1) -> hasIc(tStep)
 * int element -> giveNode(1) -> giveDof(1) -> giveBcIdValue()
 * double element -> giveNode(1) -> giveDof(1) -> giveBcValue(UnknownMode_Total, tStep)
 * DofIDItem element -> giveNode(1) -> giveDof(1) -> giveDofID()
 *
 * int element -> giveNode(1) -> hasLocalCS()
 * FloatMatrix* element -> giveNode(1) -> giveLocalCoordinateTriplet ()
 * prvni dva radky je normalizovany vstup za lcs
 *
 * BoundaryCondition* = element -> giveNode(1) -> giveDof(1) -> giveBc()
 * InitialCondition* = element -> giveNode(1) -> giveDof(1) -> giveIc()
 *
 * int element -> giveNode(1) -> giveDof(1) -> giveBc() -> isImposed(tStep)
 */

void RefinedElement :: error(const char *file, int line, const char *format, ...) const
{
    char buffer [ MAX_ERROR_MSG_LENGTH ];
    va_list args;

    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    __OOFEM_ERROR3(file, line, "Class: RefinedElement, number: %d\n%s", this->elementId, buffer);
}
} // end namespace oofem
