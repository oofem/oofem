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

#include "refinedmesh.h"
#include "domain.h"
#include "element.h"
#include "node.h"
#include "refinedelement.h"

namespace oofem {
#define EDGE_ELEM          1
#define FACE_ELEM          2
#define QUAD_ELEM          3
#define TETRA_ELEM         4
#define HEXA_ELEM          5


#define CUMUL_EDGES          ( fe_edges )
#define CUMUL_FACES          ( fe_edges + fe_faces )
#define CUMUL_QUADS          ( fe_edges + fe_faces + fe_quads )
#define CUMUL_TETRAS         ( fe_edges + fe_faces + fe_quads + fe_tetras )
#define CUMUL_HEXAS          ( fe_edges + fe_faces + fe_quads + fe_tetras + fe_hexas )


/*  < : C style */
#define elem_type(ID)    ( ( ( ID ) <= CUMUL_EDGES ) ? EDGE_ELEM :    \
                          ( ( ID ) <= CUMUL_FACES ) ? FACE_ELEM :    \
                          ( ( ID ) <= CUMUL_QUADS ) ? QUAD_ELEM :    \
                          ( ( ID ) <= CUMUL_TETRAS ) ? TETRA_ELEM :  \
                          ( ( ID ) <= CUMUL_HEXAS ) ? HEXA_ELEM : 0 )

#define is_edge(ID)       ( ( elem_type(ID) == EDGE_ELEM ) ? 1 : 0 )
#define is_face(ID)       ( ( elem_type(ID) == FACE_ELEM ) ? 1 : 0 )
#define is_quad(ID)       ( ( elem_type(ID) == QUAD_ELEM ) ? 1 : 0 )
#define is_tetra(ID)      ( ( elem_type(ID) == TETRA_ELEM ) ? 1 : 0 )
#define is_hexa(ID)       ( ( elem_type(ID) == HEXA_ELEM ) ? 1 : 0 )

/*  <, >= : C style */
/*
 * #define is_edge(ID)    (((ID) <= CUMUL_EDGES && (ID) > 0) ? 1 : 0)
 * #define is_face(ID)    (((ID) <= CUMUL_FACES && (ID) > CUMUL_EDGES) ? 1 : 0)
 * #define is_quad(ID)    (((ID) <= CUMUL_QUADS && (ID) > CUMUL_FACES) ? 1 : 0)
 * #define is_tetra(ID)   (((ID) <= CUMUL_TETRAS && (ID) > CUMUL_QUADS) ? 1 : 0)
 * #define is_hexa(ID)    (((ID) <= CUMUL_HEXAS && (ID) > CUMUL_TETRAS) ? 1 : 0)
 */

#define global_edge_id(ID)    ( ( ID ) )
#define global_face_id(ID)    ( ( ID ) + CUMUL_EDGES )
#define global_quad_id(ID)    ( ( ID ) + CUMUL_FACES )
#define global_tetra_id(ID)   ( ( ID ) + CUMUL_QUADS )
#define global_hexa_id(ID)    ( ( ID ) + CUMUL_TETRAS )


#define local_edge_id(ID)     ( ( ID ) )
#define local_face_id(ID)     ( ( ID ) -CUMUL_EDGES )
#define local_quad_id(ID)     ( ( ID ) -CUMUL_FACES )
#define local_tetra_id(ID)    ( ( ID ) -CUMUL_QUADS )
#define local_hexa_id(ID)     ( ( ID ) -CUMUL_TETRAS )


#define matrix_2d(ARRAY, U, V)      ( ARRAY ) [ ( level + 2 ) * ( V ) + ( U ) ]
#define matrix_3d(ARRAY, U, V, W)    ( ARRAY ) [ ( level + 2 ) * ( level + 2 ) * ( W ) + ( level + 2 ) * ( V ) + ( U ) ]


#define error_message(MSG)  { OOFEM_LOG_RELEVANT( "refineMeshGlobally: %s\n", ( MSG ) ); return ( -1 ); }


/* I do my own local and global numbering;
 * the reason is to mimic t3d output file numbering;
 *
 * local numbering: each type of elem (edge, face, quad, tetra, hexa) is numbered separately
 *                  but respecting the relative position in global element list
 * global numbering: local numbering with element types ordered as edge face quad tetra hexa
 *
 * example:
 * no:    1   2   3   4   5   6   7   8   9
 * tp:    t   q   h   f   t   f   h   q   f
 *
 * local numbering:  faces  1(4) 2(6) 3(9)
 *                   quads  1(2) 2(8)
 *                   tetras 1(1) 2(5)
 *                   hexas  1(3) 2(7)
 * global numbering: 1(4) 2(6) 3(9) 4(2) 5(8) 6(1) 7(5) 8(3) 9(7) */



int
RefinedMesh :: refineMeshGlobally(Domain *d, int level, std :: vector< RefinedElement > &refinedElementList)
{
    int jj, kk = 0, j1, j2, j3, j4, type;
    long node, elem, nd, pos = 0, p, number, size;
    long i, j, k, m, n, node1, node2, node3, node4, elem1, elem2;
    long refine_nodes, fine_node_id, refine_node_id;
    long mesh_edge_id, fine_edge_id, fine_quad_id, fine_hexa_id;
    long *tmp_array_start [ 4 ], *tmp_array_end [ 4 ];
    long *tmp_array_cen [ 6 ], *tmp_array [ 12 ];
    int swap [ 6 ], flag;

    Element *element = NULL;
    IntArray *connectivity = NULL, *boundary = NULL;

    long fe_nodes = 0, fe_elems = 0;
    long fe_edges = 0, fe_faces = 0, fe_quads = 0, fe_tetras = 0, fe_hexas = 0, fe_pyrams = 0, fe_wedges = 0;
    long edge_id = 0, face_id = 0, quad_id = 0, tetra_id = 0, hexa_id = 0; // pyram_id = 0, wedge_id = 0;
    long mesh_edges = 0, mesh_faces = 0, mesh_quads = 0;
    long fine_nodes = 0, fine_edges = 0, fine_quads = 0, fine_hexas = 0;

    fe_node_rec *fe_node = NULL, *fe_node_array = NULL;
    fe_edge_rec *fe_edge = NULL, *fe_edge_array = NULL;
    fe_face_rec *fe_face = NULL, *fe_face_array = NULL;
    fe_quad_rec *fe_quad = NULL, *fe_quad_array = NULL;
    fe_tetra_rec *fe_tetra = NULL, *fe_tetra_array = NULL;
    fe_hexa_rec *fe_hexa = NULL, *fe_hexa_array = NULL;

    fe_node_rec *fine_node = NULL, *fine_node_array = NULL;
    /*
     * fe_node_rec *refine_node = NULL, *refine_node_array = NULL;
     */
    fe_node_rec *nd1 = NULL, *nd2 = NULL, *nd3 = NULL, *nd4 = NULL;

    mesh_edge_rec *mesh_edge = NULL, *mesh_edge_array = NULL;
    mesh_face_rec *mesh_face = NULL, *mesh_face_array = NULL, *mesh_fc [ 4 ];
    mesh_quad_rec *mesh_quad = NULL, *mesh_quad_array = NULL, *mesh_qd [ 6 ];

#ifdef COMPLETE_DATA_STRUCTURE
    tmp_face_rec *tmp_face = NULL, *tmp_face_array = NULL, *tmp_fc = NULL;
    tmp_quad_rec *tmp_quad = NULL, *tmp_quad_array = NULL, *tmp_qd = NULL;
#endif

    tmp_tetra_rec *tmp_tetra = NULL, *tmp_tetra_array = NULL, *tmp_tet = NULL;
    tmp_hexa_rec *tmp_hexa = NULL, *tmp_hexa_array = NULL, *tmp_hx = NULL;

    fine_edge_rec *fine_edge = NULL, *fine_edge_array = NULL;
    fine_quad_rec *fine_quad = NULL, *fine_quad_array = NULL;
    fine_hexa_rec *fine_hexa = NULL, *fine_hexa_array = NULL;

    long *node_num_elems = NULL, *node_con_elems = NULL;
    long *node_num_nodes = NULL, *node_con_nodes = NULL;

    /* caution: side and face ordering on face and tetra is not consistent with T3D */

    short face_ed_nd [ 3 ] [ 2 ] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
    short quad_ed_nd [ 4 ] [ 2 ] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };

    /* note: ordering of nodes on tetra faces is important (normal: inner, outer, outer, outer) */

    short tetra_fc_nd [ 4 ] [ 3 ] = { { 0, 1, 2 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 0, 3 } };

    /* note: ordering of nodes on hexa faces is important (normal: inner, outer, outer, outer, outer, inner) */

    short hexa_fc_nd [ 6 ] [ 4 ] = { { 0, 1, 2, 3 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 }, { 7, 6, 5, 4 } };

    short quad_con_nd [ 4 ] [ 2 ] = { { 1, 3 }, { 2, 0 }, { 3, 1 }, { 0, 2 } };
    short hexa_con_nd [ 8 ] [ 3 ] = { { 1, 3, 4 }, { 2, 0, 5 }, { 3, 1, 6 }, { 0, 2, 7 }, { 7, 5, 0 }, { 4, 6, 1 }, { 5, 7, 2 }, { 6, 4, 3 } };

    fe_nodes = d->giveNumberOfDofManagers();
    fe_elems = d->giveNumberOfElements();

    for ( i = 0; i < fe_elems; i++ ) {
        element = d->giveElement(i + 1);
        switch ( element->giveGeometryType() ) {
        case EGT_line_1:
            fe_edges++;
            break;
        case EGT_triangle_1:
            fe_faces++;
            break;
        case EGT_quad_1:
            fe_quads++;
            break;
        case EGT_tetra_1:
            fe_tetras++;
            break;
        case EGT_hexa_1:
            fe_hexas++;
            break;
        default:
            error_message("Not supported element type");
            break;
        }
    }

    if ( fe_nodes == 0 ) {
        error_message("Not enough nodes found");
    }

    if ( fe_nodes != 0 ) {
        if ( ( fe_node_array = ( fe_node_rec * ) calloc( fe_nodes, sizeof( fe_node_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( fe_edges != 0 ) {
        if ( ( fe_edge_array = ( fe_edge_rec * ) calloc( fe_edges, sizeof( fe_edge_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( fe_faces != 0 ) {
        if ( ( fe_face_array = ( fe_face_rec * ) calloc( fe_faces, sizeof( fe_face_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( fe_quads != 0 ) {
        if ( ( fe_quad_array = ( fe_quad_rec * ) calloc( fe_quads, sizeof( fe_quad_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( fe_tetras != 0 ) {
        if ( ( fe_tetra_array = ( fe_tetra_rec * ) calloc( fe_tetras, sizeof( fe_tetra_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( fe_hexas != 0 ) {
        if ( ( fe_hexa_array = ( fe_hexa_rec * ) calloc( fe_hexas, sizeof( fe_hexa_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    fe_node = fe_node_array;
    for ( i = 0; i < fe_nodes; i++, fe_node++ ) {
        fe_node->id = i + 1;
    }

    edge_id = face_id = quad_id = tetra_id = hexa_id = 0;
    for ( i = 0; i < fe_elems; i++ ) {
        element = d->giveElement(i + 1);
        switch ( element->giveGeometryType() ) {
        case EGT_line_1:
            fe_edge = & ( fe_edge_array [ edge_id++ ] );
            for ( j = 0; j < 2; j++ ) {
                nd = element->giveNode(j + 1)->giveNumber();
                if ( nd <= 0 || nd > fe_nodes ) {
                    error_message("Invalid edge nodes");
                }

                fe_edge->node [ j ] = & ( fe_node_array [ nd - 1 ] );
            }

            break;
        case EGT_triangle_1:
            fe_face = & ( fe_face_array [ face_id++ ] );
            for ( j = 0; j < 3; j++ ) {
                nd = element->giveNode(j + 1)->giveNumber();
                if ( nd <= 0 || nd > fe_nodes ) {
                    error_message("Invalid face nodes");
                }

                fe_face->node [ j ] = & ( fe_node_array [ nd - 1 ] );
            }

            break;
        case EGT_quad_1:
            fe_quad = & ( fe_quad_array [ quad_id++ ] );
            for ( j = 0; j < 4; j++ ) {
                nd = element->giveNode(j + 1)->giveNumber();
                if ( nd <= 0 || nd > fe_nodes ) {
                    error_message("Invalid quad nodes");
                }

                fe_quad->node [ j ] = & ( fe_node_array [ nd - 1 ] );
            }

            break;
        case EGT_tetra_1:
            fe_tetra = & ( fe_tetra_array [ tetra_id++ ] );
            for ( j = 0; j < 4; j++ ) {
                nd = element->giveNode(j + 1)->giveNumber();
                if ( nd <= 0 || nd > fe_nodes ) {
                    error_message("Invalid tetra nodes");
                }

                fe_tetra->node [ j ] = & ( fe_node_array [ nd - 1 ] );
            }

            break;
        case EGT_hexa_1:
            fe_hexa = & ( fe_hexa_array [ hexa_id++ ] );
            for ( j = 0; j < 8; j++ ) {
                nd = element->giveNode(j + 1)->giveNumber();
                if ( nd <= 0 || nd > fe_nodes ) {
                    error_message("Invalid hexa nodes");
                }

                fe_hexa->node [ j ] = & ( fe_node_array [ nd - 1 ] );
            }

            break;
        default:
            error_message("Not supported element type");
            break;
        }
    }

    /* count number of elements incident to nodes */

    if ( ( node_num_elems = ( long * ) calloc( fe_nodes + 1, sizeof( long ) ) ) == NULL ) {
        error_message("Memory allocation error");
    }

    fe_edge = fe_edge_array;
    for ( i = 0; i < fe_edges; i++, fe_edge++ ) {
        for ( j = 0; j < 2; j++ ) {
            nd = fe_edge->node [ j ]->id;
            node_num_elems [ nd - 1 ]++;
        }
    }

    fe_face = fe_face_array;
    for ( i = 0; i < fe_faces; i++, fe_face++ ) {
        for ( j = 0; j < 3; j++ ) {
            nd = fe_face->node [ j ]->id;
            node_num_elems [ nd - 1 ]++;
        }
    }

    fe_quad = fe_quad_array;
    for ( i = 0; i < fe_quads; i++, fe_quad++ ) {
        for ( j = 0; j < 4; j++ ) {
            nd = fe_quad->node [ j ]->id;
            node_num_elems [ nd - 1 ]++;
        }
    }

    fe_tetra = fe_tetra_array;
    for ( i = 0; i < fe_tetras; i++, fe_tetra++ ) {
        for ( j = 0; j < 4; j++ ) {
            nd = fe_tetra->node [ j ]->id;
            node_num_elems [ nd - 1 ]++;
        }
    }

    fe_hexa = fe_hexa_array;
    for ( i = 0; i < fe_hexas; i++, fe_hexa++ ) {
        for ( j = 0; j < 8; j++ ) {
            nd = fe_hexa->node [ j ]->id;
            node_num_elems [ nd - 1 ]++;
        }
    }

    /* recalculate the number of elements to current addresses */

    pos = 0;
    for ( i = 0; i < fe_nodes; i++ ) {
        number = node_num_elems [ i ];
        node_num_elems [ i ] = pos;
        pos += number;
    }

    node_num_elems [ fe_nodes ] = size = pos;
    if ( ( node_con_elems = ( long * ) calloc( size, sizeof( long ) ) ) == NULL ) {
        error_message("Memory allocation error");
    }

    /* store element numbers incident to nodes */

    fe_edge = fe_edge_array;
    for ( i = 0; i < fe_edges; i++, fe_edge++ ) {
        for ( j = 0; j < 2; j++ ) {
            nd = fe_edge->node [ j ]->id;
            node_con_elems [ node_num_elems [ nd - 1 ]++ ] = i + 1;
        }
    }

    fe_face = fe_face_array;
    for ( i = 0; i < fe_faces; i++, fe_face++ ) {
        for ( j = 0; j < 3; j++ ) {
            nd = fe_face->node [ j ]->id;
            node_con_elems [ node_num_elems [ nd - 1 ]++ ] = i + 1 + fe_edges;
        }
    }

    fe_quad = fe_quad_array;
    for ( i = 0; i < fe_quads; i++, fe_quad++ ) {
        for ( j = 0; j < 4; j++ ) {
            nd = fe_quad->node [ j ]->id;
            node_con_elems [ node_num_elems [ nd - 1 ]++ ] = i + 1 + fe_edges + fe_faces;
        }
    }

    fe_tetra = fe_tetra_array;
    for ( i = 0; i < fe_tetras; i++, fe_tetra++ ) {
        for ( j = 0; j < 4; j++ ) {
            nd = fe_tetra->node [ j ]->id;
            node_con_elems [ node_num_elems [ nd - 1 ]++ ] = i + 1 + fe_edges + fe_faces + fe_quads;
        }
    }

    fe_hexa = fe_hexa_array;
    for ( i = 0; i < fe_hexas; i++, fe_hexa++ ) {
        for ( j = 0; j < 8; j++ ) {
            nd = fe_hexa->node [ j ]->id;
            node_con_elems [ node_num_elems [ nd - 1 ]++ ] = i + 1 + fe_edges + fe_faces + fe_quads + fe_tetras;
        }
    }

    /* recalculate the addresses to address of the first element */

    pos = 0;
    for ( i = 0; i < fe_nodes; i++ ) {
        number = node_num_elems [ i ] - pos;
        node_num_elems [ i ] = pos;
        pos += number;
    }

#ifdef DEBUG
    for ( i = 0; i < fe_nodes; i++ ) {
        OOFEM_LOG_DEBUG("node %ld: %ld:", i + 1, node_num_elems [ i + 1 ] - node_num_elems [ i ]);
        for ( j = node_num_elems [ i ]; j < node_num_elems [ i + 1 ]; j++ ) {
            switch ( elem_type(node_con_elems [ j ]) ) {
            case EDGE_ELEM:
                OOFEM_LOG_DEBUG( " %ld (e)", local_edge_id(node_con_elems [ j ]) );
                break;
            case FACE_ELEM:
                OOFEM_LOG_DEBUG( " %ld (f)", local_face_id(node_con_elems [ j ]) );
                break;
            case QUAD_ELEM:
                OOFEM_LOG_DEBUG( " %ld (q)", local_quad_id(node_con_elems [ j ]) );
                break;
            case TETRA_ELEM:
                OOFEM_LOG_DEBUG( " %ld (t)", local_tetra_id(node_con_elems [ j ]) );
                break;
            case HEXA_ELEM:
                OOFEM_LOG_DEBUG( " %ld (h)", local_hexa_id(node_con_elems [ j ]) );
                break;
            default:
                OOFEM_LOG_DEBUG(" %ld (?)", node_con_elems [ j ]);
            }
        }

        OOFEM_LOG_DEBUG("\n");
    }

    OOFEM_LOG_DEBUG("\n");
#endif

    /* count number of nodes incident to nodes */

    if ( ( node_num_nodes = ( long * ) calloc( fe_nodes + 1, sizeof( long ) ) ) == NULL ) {
        error_message("Memory allocation error");
    }

    fe_node = fe_node_array;
    for ( i = 0; i < fe_nodes; i++, fe_node++ ) {
        node = i + 1;
        for ( j = node_num_elems [ i ]; j < node_num_elems [ node ]; j++ ) {
            elem = node_con_elems [ j ];
            type = elem_type(elem);                     /* macro */
            switch ( type ) {
            case EDGE_ELEM:
                elem = local_edge_id(elem);            /* macro */
                fe_edge = & ( fe_edge_array [ elem - 1 ] );
                for ( k = 0; k < 2; k++ ) {
                    nd = fe_edge->node [ k ]->id;
                    if ( nd != node && nd > 0 ) {
                        fe_edge->node [ k ]->id = -nd;
                        node_num_nodes [ i ]++;
                    }
                }

                break;
            case FACE_ELEM:
                elem = local_face_id(elem);            /* macro */
                fe_face = & ( fe_face_array [ elem - 1 ] );
                for ( k = 0; k < 3; k++ ) {
                    nd = fe_face->node [ k ]->id;
                    if ( nd != node && nd > 0 ) {
                        fe_face->node [ k ]->id = -nd;
                        node_num_nodes [ i ]++;
                    }
                }

                break;
            case QUAD_ELEM:
                elem = local_quad_id(elem);            /* macro */
                fe_quad = & ( fe_quad_array [ elem - 1 ] );
                for ( k = 0; k < 4; k++ ) {
                    if ( fe_node == fe_quad->node [ k ] ) {
                        for ( m = 0; m < 2; m++ ) {
                            nd = fe_quad->node [ quad_con_nd [ k ] [ m ] ]->id;
                            if ( nd > 0 ) {
                                fe_quad->node [ quad_con_nd [ k ] [ m ] ]->id = -nd;
                                node_num_nodes [ i ]++;
                            }
                        }

                        break;
                    }
                }

                break;
            case TETRA_ELEM:
                elem = local_tetra_id(elem);           /* macro */
                fe_tetra = & ( fe_tetra_array [ elem - 1 ] );
                for ( k = 0; k < 4; k++ ) {
                    nd = fe_tetra->node [ k ]->id;
                    if ( nd != node && nd > 0 ) {
                        fe_tetra->node [ k ]->id = -nd;
                        node_num_nodes [ i ]++;
                    }
                }

                break;
            case HEXA_ELEM:
                elem = local_hexa_id(elem);            /* macro */
                fe_hexa = & ( fe_hexa_array [ elem - 1 ] );
                for ( k = 0; k < 8; k++ ) {
                    if ( fe_node == fe_hexa->node [ k ] ) {
                        for ( m = 0; m < 3; m++ ) {
                            nd = fe_hexa->node [ hexa_con_nd [ k ] [ m ] ]->id;
                            if ( nd > 0 ) {
                                fe_hexa->node [ hexa_con_nd [ k ] [ m ] ]->id = -nd;
                                node_num_nodes [ i ]++;
                            }
                        }

                        break;
                    }
                }

                break;
            default:
                error_message("Invalid element number");
                break;
            }
        }

        /* restore connected node id to positive value */

        for ( j = node_num_elems [ i ]; j < node_num_elems [ node ]; j++ ) {
            elem = node_con_elems [ j ];
            type = elem_type(elem);                     /* macro */
            switch ( type ) {
            case EDGE_ELEM:
                elem = local_edge_id(elem);            /* macro */
                fe_edge = & ( fe_edge_array [ elem - 1 ] );
                for ( k = 0; k < 2; k++ ) {
                    fe_edge->node [ k ]->id = abs(fe_edge->node [ k ]->id);
                }

                break;
            case FACE_ELEM:
                elem = local_face_id(elem);            /* macro */
                fe_face = & ( fe_face_array [ elem - 1 ] );
                for ( k = 0; k < 3; k++ ) {
                    fe_face->node [ k ]->id = abs(fe_face->node [ k ]->id);
                }

                break;
            case QUAD_ELEM:
                elem = local_quad_id(elem);            /* macro */
                fe_quad = & ( fe_quad_array [ elem - 1 ] );
                for ( k = 0; k < 4; k++ ) {
                    fe_quad->node [ k ]->id = abs(fe_quad->node [ k ]->id);
                }

                break;
            case TETRA_ELEM:
                elem = local_tetra_id(elem);           /* macro */
                fe_tetra = & ( fe_tetra_array [ elem - 1 ] );
                for ( k = 0; k < 4; k++ ) {
                    fe_tetra->node [ k ]->id = abs(fe_tetra->node [ k ]->id);
                }

                break;
            case HEXA_ELEM:
                elem = local_hexa_id(elem);            /* macro */
                fe_hexa = & ( fe_hexa_array [ elem - 1 ] );
                for ( k = 0; k < 8; k++ ) {
                    fe_hexa->node [ k ]->id = abs(fe_hexa->node [ k ]->id);
                }

                break;
            default:
                error_message("Invalid element number");
                break;
            }
        }
    }

    /* recalculate the number of nodes to current addresses */

    pos = 0;
    for ( i = 0; i < fe_nodes; i++ ) {
        number = node_num_nodes [ i ];
        node_num_nodes [ i ] = pos;
        pos += number;
    }

    node_num_nodes [ fe_nodes ] = size = pos;
    if ( ( node_con_nodes = ( long * ) calloc( size, sizeof( long ) ) ) == NULL ) {
        error_message("Memory allocation error");
    }

    mesh_edges = size / 2;

    /* store nodes incident to nodes */

    fe_node = fe_node_array;
    for ( i = 0; i < fe_nodes; i++, fe_node++ ) {
        node = i + 1;
        for ( j = node_num_elems [ i ]; j < node_num_elems [ node ]; j++ ) {
            elem = node_con_elems [ j ];
            type = elem_type(elem);                     /* macro */
            switch ( type ) {
            case EDGE_ELEM:
                elem = local_edge_id(elem);            /* macro */
                fe_edge = & ( fe_edge_array [ elem - 1 ] );
                for ( k = 0; k < 2; k++ ) {
                    nd = fe_edge->node [ k ]->id;
                    if ( nd != node && nd > 0 ) {
                        fe_edge->node [ k ]->id = -nd;
                        node_con_nodes [ node_num_nodes [ i ]++ ] = nd;
                    }
                }

                break;
            case FACE_ELEM:
                elem = local_face_id(elem);            /* macro */
                fe_face = & ( fe_face_array [ elem - 1 ] );
                for ( k = 0; k < 3; k++ ) {
                    nd = fe_face->node [ k ]->id;
                    if ( nd != node && nd > 0 ) {
                        fe_face->node [ k ]->id = -nd;
                        node_con_nodes [ node_num_nodes [ i ]++ ] = nd;
                    }
                }

                break;
            case QUAD_ELEM:
                elem = local_quad_id(elem);            /* macro */
                fe_quad = & ( fe_quad_array [ elem - 1 ] );
                for ( k = 0; k < 4; k++ ) {
                    if ( fe_node == fe_quad->node [ k ] ) {
                        for ( m = 0; m < 2; m++ ) {
                            nd = fe_quad->node [ quad_con_nd [ k ] [ m ] ]->id;
                            if ( nd > 0 ) {
                                fe_quad->node [ quad_con_nd [ k ] [ m ] ]->id = -nd;
                                node_con_nodes [ node_num_nodes [ i ]++ ] = nd;
                            }
                        }

                        break;
                    }
                }

                break;
            case TETRA_ELEM:
                elem = local_tetra_id(elem);           /* macro */
                fe_tetra = & ( fe_tetra_array [ elem - 1 ] );
                for ( k = 0; k < 4; k++ ) {
                    nd = fe_tetra->node [ k ]->id;
                    if ( nd != node && nd > 0 ) {
                        fe_tetra->node [ k ]->id = -nd;
                        node_con_nodes [ node_num_nodes [ i ]++ ] = nd;
                    }
                }

                break;
            case HEXA_ELEM:
                elem = local_hexa_id(elem);            /* macro */
                fe_hexa = & ( fe_hexa_array [ elem - 1 ] );
                for ( k = 0; k < 8; k++ ) {
                    if ( fe_node == fe_hexa->node [ k ] ) {
                        for ( m = 0; m < 3; m++ ) {
                            nd = fe_hexa->node [ hexa_con_nd [ k ] [ m ] ]->id;
                            if ( nd > 0 ) {
                                fe_hexa->node [ hexa_con_nd [ k ] [ m ] ]->id = -nd;
                                node_con_nodes [ node_num_nodes [ i ]++ ] = nd;
                            }
                        }

                        break;
                    }
                }

                break;
            default:
                error_message("Invalid element number");
                break;
            }
        }

        /* restore connected node id to positive value */

        for ( j = node_num_elems [ i ]; j < node_num_elems [ node ]; j++ ) {
            elem = node_con_elems [ j ];
            type = elem_type(elem);                     /* macro */
            switch ( type ) {
            case EDGE_ELEM:
                elem = local_edge_id(elem);            /* macro */
                fe_edge = & ( fe_edge_array [ elem - 1 ] );
                for ( k = 0; k < 2; k++ ) {
                    fe_edge->node [ k ]->id = abs(fe_edge->node [ k ]->id);
                }

                break;
            case FACE_ELEM:
                elem = local_face_id(elem);            /* macro */
                fe_face = & ( fe_face_array [ elem - 1 ] );
                for ( k = 0; k < 3; k++ ) {
                    fe_face->node [ k ]->id = abs(fe_face->node [ k ]->id);
                }

                break;
            case QUAD_ELEM:
                elem = local_quad_id(elem);            /* macro */
                fe_quad = & ( fe_quad_array [ elem - 1 ] );
                for ( k = 0; k < 4; k++ ) {
                    fe_quad->node [ k ]->id = abs(fe_quad->node [ k ]->id);
                }

                break;
            case TETRA_ELEM:
                elem = local_tetra_id(elem);           /* macro */
                fe_tetra = & ( fe_tetra_array [ elem - 1 ] );
                for ( k = 0; k < 4; k++ ) {
                    fe_tetra->node [ k ]->id = abs(fe_tetra->node [ k ]->id);
                }

                break;
            case HEXA_ELEM:
                elem = local_hexa_id(elem);            /* macro */
                fe_hexa = & ( fe_hexa_array [ elem - 1 ] );
                for ( k = 0; k < 8; k++ ) {
                    fe_hexa->node [ k ]->id = abs(fe_hexa->node [ k ]->id);
                }

                break;
            default:
                error_message("Invalid element number");
                break;
            }
        }
    }

    /* recalculate the addresses to address of the first node */

    pos = 0;
    for ( i = 0; i < fe_nodes; i++ ) {
        number = node_num_nodes [ i ] - pos;
        node_num_nodes [ i ] = pos;
        pos += number;
    }

#ifdef DEBUG
    for ( i = 0; i < fe_nodes; i++ ) {
        OOFEM_LOG_DEBUG("node %ld: %ld:", i + 1, node_num_nodes [ i + 1 ] - node_num_nodes [ i ]);
        for ( j = node_num_nodes [ i ]; j < node_num_nodes [ i + 1 ]; j++ ) {
            OOFEM_LOG_DEBUG(" %ld", node_con_nodes [ j ]);
        }

        OOFEM_LOG_DEBUG("\n");
    }

    OOFEM_LOG_DEBUG("\n");
#endif

#ifdef COMPLETE_DATA_STRUCTURE
    if ( fe_faces != 0 ) {
        if ( ( tmp_face_array = ( tmp_face_rec * ) calloc( fe_faces, sizeof( tmp_face_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( fe_quads != 0 ) {
        if ( ( tmp_quad_array = ( tmp_quad_rec * ) calloc( fe_quads, sizeof( tmp_quad_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

#endif

    if ( fe_tetras != 0 ) {
        if ( ( tmp_tetra_array = ( tmp_tetra_rec * ) calloc( fe_tetras, sizeof( tmp_tetra_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( fe_hexas != 0 ) {
        if ( ( tmp_hexa_array = ( tmp_hexa_rec * ) calloc( fe_hexas, sizeof( tmp_hexa_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

#ifdef COMPLETE_DATA_STRUCTURE

    /* setup 2D ngb elems */

    fe_face = fe_face_array;
    for ( i = 0; i < fe_faces; i++, fe_face++ ) {
        elem1 = global_face_id(i) + 1;                  /* macro */
        for ( j = 0; j < 3; j++ ) {
            j1 = face_ed_nd [ j ] [ 0 ];
            j2 = face_ed_nd [ j ] [ 1 ];
            node1 = fe_face->node [ j1 ]->id;
            node2 = fe_face->node [ j2 ]->id;
            for ( k = node_num_elems [ node1 - 1 ]; k < node_num_elems [ node1 ]; k++ ) {
                elem2 = node_con_elems [ k ];
                if ( elem2 == elem1 ) {
                    continue;
                }

                if ( is_face(elem2) == 1 || is_quad(elem2) == 1 ) {           /* macro */
                    for ( m = node_num_elems [ node2 - 1 ]; m < node_num_elems [ node2 ]; m++ ) {
                        if ( elem2 == node_con_elems [ m ] ) {
                            /* this error message is generated to prevent multiple creation of the non-manifold edge;
                             * therefore do not break the outer for cycle to detect the non-manifold edge */

                            if ( ( & ( tmp_face_array [ i ] ) )->ngb_elem_id [ j ] != 0 ) {
                                error_message("Domain is not 2-manifold");
                            }

                            ( & ( tmp_face_array [ i ] ) )->ngb_elem_id [ j ] = elem2;
                            break;
                        }
                    }
                }
            }
        }
    }

    fe_quad = fe_quad_array;
    for ( i = 0; i < fe_quads; i++, fe_quad++ ) {
        elem1 = global_quad_id(i) + 1;                  /* macro */
        for ( j = 0; j < 4; j++ ) {
            j1 = quad_ed_nd [ j ] [ 0 ];
            j2 = quad_ed_nd [ j ] [ 1 ];
            node1 = fe_quad->node [ j1 ]->id;
            node2 = fe_quad->node [ j2 ]->id;
            for ( k = node_num_elems [ node1 - 1 ]; k < node_num_elems [ node1 ]; k++ ) {
                elem2 = node_con_elems [ k ];
                if ( elem2 == elem1 ) {
                    continue;
                }

                if ( is_face(elem2) == 1 || is_quad(elem2) == 1 ) {           /* macro */
                    for ( m = node_num_elems [ node2 - 1 ]; m < node_num_elems [ node2 ]; m++ ) {
                        if ( elem2 == node_con_elems [ m ] ) {
                            /* this error message is generated to prevent multiple creation of the non-manifold edge;
                             * therefore do not break the outer for cycle to detect the non-manifold edge */

                            if ( ( & ( tmp_quad_array [ i ] ) )->ngb_elem_id [ j ] != 0 ) {
                                error_message("Domain is not 2-manifold");
                            }

                            ( & ( tmp_quad_array [ i ] ) )->ngb_elem_id [ j ] = elem2;
                            break;
                        }
                    }
                }
            }
        }
    }

 #ifdef DEBUG
    tmp_face = tmp_face_array;
    for ( i = 0; i < fe_faces; i++, tmp_face++ ) {
        OOFEM_LOG_DEBUG("face %ld: ", i + 1);
        for ( j = 0; j < 3; j++ ) {
            switch ( elem_type(tmp_face->ngb_elem_id [ j ]) ) {
            case FACE_ELEM:
                OOFEM_LOG_DEBUG( " %ld (f)", local_face_id(tmp_face->ngb_elem_id [ j ]) );
                break;
            case QUAD_ELEM:
                OOFEM_LOG_DEBUG( " %ld (q)", local_quad_id(tmp_face->ngb_elem_id [ j ]) );
                break;
            }
        }

        OOFEM_LOG_DEBUG("\n");
    }

    tmp_quad = tmp_quad_array;
    for ( i = 0; i < fe_quads; i++, tmp_quad++ ) {
        OOFEM_LOG_DEBUG("quad %ld: ", i + 1);
        for ( j = 0; j < 4; j++ ) {
            switch ( elem_type(tmp_quad->ngb_elem_id [ j ]) ) {
            case FACE_ELEM:
                OOFEM_LOG_DEBUG( " %ld (f)", local_face_id(tmp_quad->ngb_elem_id [ j ]) );
                break;
            case QUAD_ELEM:
                OOFEM_LOG_DEBUG( " %ld (q)", local_quad_id(tmp_quad->ngb_elem_id [ j ]) );
                break;
            }
        }

        OOFEM_LOG_DEBUG("\n");
    }

    OOFEM_LOG_DEBUG("\n");
 #endif

#endif

    /* setup 3D ngb elems */

    fe_tetra = fe_tetra_array;
    for ( i = 0; i < fe_tetras; i++, fe_tetra++ ) {
        elem1 = global_tetra_id(i) + 1;                          /* macro */
        mesh_faces += 4;
        for ( j = 0; j < 4; j++ ) {
            j1 = tetra_fc_nd [ j ] [ 0 ];
            j2 = tetra_fc_nd [ j ] [ 1 ];
            j3 = tetra_fc_nd [ j ] [ 2 ];
            node1 = fe_tetra->node [ j1 ]->id;
            node2 = fe_tetra->node [ j2 ]->id;
            node3 = fe_tetra->node [ j3 ]->id;
            for ( k = node_num_elems [ node1 - 1 ]; k < node_num_elems [ node1 ]; k++ ) {
                elem2 = node_con_elems [ k ];
                if ( elem2 == elem1 ) {
                    continue;
                }

                if ( is_tetra(elem2) == 0 ) {
                    continue;                                                   /* macro */
                }

                for ( m = node_num_elems [ node2 - 1 ]; m < node_num_elems [ node2 ]; m++ ) {
                    if ( elem2 != node_con_elems [ m ] ) {
                        continue;
                    }

                    for ( n = node_num_elems [ node3 - 1 ]; n < node_num_elems [ node3 ]; n++ ) {
                        if ( elem2 == node_con_elems [ n ] ) {
                            ( & ( tmp_tetra_array [ i ] ) )->ngb_elem_id [ j ] = elem2;
                            if ( elem1 > elem2 ) {
                                mesh_faces--;
                            }

                            k = node_num_elems [ node1 ];
                            m = node_num_elems [ node2 ];
                            break;
                        }
                    }
                }
            }
        }
    }

    fe_hexa = fe_hexa_array;
    for ( i = 0; i < fe_hexas; i++, fe_hexa++ ) {
        elem1 = global_hexa_id(i) + 1;                           /* macro */
        mesh_quads += 6;
        for ( j = 0; j < 6; j++ ) {
            j1 = hexa_fc_nd [ j ] [ 0 ];
            j2 = hexa_fc_nd [ j ] [ 1 ];
            j3 = hexa_fc_nd [ j ] [ 2 ];
            j4 = hexa_fc_nd [ j ] [ 3 ];
            node1 = fe_hexa->node [ j1 ]->id;
            node2 = fe_hexa->node [ j2 ]->id;
            node3 = fe_hexa->node [ j3 ]->id;
            node4 = fe_hexa->node [ j4 ]->id;
            for ( k = node_num_elems [ node1 - 1 ]; k < node_num_elems [ node1 ]; k++ ) {
                elem2 = node_con_elems [ k ];
                if ( elem2 == elem1 ) {
                    continue;
                }

                if ( is_hexa(elem2) == 0 ) {
                    continue;                                                 /* macro */
                }

                for ( m = node_num_elems [ node2 - 1 ]; m < node_num_elems [ node2 ]; m++ ) {
                    if ( elem2 != node_con_elems [ m ] ) {
                        continue;
                    }

                    for ( n = node_num_elems [ node3 - 1 ]; n < node_num_elems [ node3 ]; n++ ) {
                        if ( elem2 != node_con_elems [ n ] ) {
                            continue;
                        }

                        for ( p = node_num_elems [ node4 - 1 ]; p < node_num_elems [ node4 ]; p++ ) {
                            if ( elem2 == node_con_elems [ p ] ) {
                                ( & ( tmp_hexa_array [ i ] ) )->ngb_elem_id [ j ] = elem2;
                                if ( elem1 > elem2 ) {
                                    mesh_quads--;
                                }

                                k = node_num_elems [ node1 ];
                                m = node_num_elems [ node2 ];
                                n = node_num_elems [ node3 ];
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

#ifdef DEBUG
    tmp_tetra = tmp_tetra_array;
    for ( i = 0; i < fe_tetras; i++, tmp_tetra++ ) {
        OOFEM_LOG_DEBUG("tetra %ld: ", i + 1);
        for ( j = 0; j < 4; j++ ) {
            if ( tmp_tetra->ngb_elem_id [ j ] == 0 ) {
                continue;
            }

            OOFEM_LOG_DEBUG( " %ld (t)", local_tetra_id(tmp_tetra->ngb_elem_id [ j ]) );
        }

        OOFEM_LOG_DEBUG("\n");
    }

    tmp_hexa = tmp_hexa_array;
    for ( i = 0; i < fe_hexas; i++, tmp_hexa++ ) {
        OOFEM_LOG_DEBUG("hexa %ld: ", i + 1);
        for ( j = 0; j < 6; j++ ) {
            if ( tmp_hexa->ngb_elem_id [ j ] == 0 ) {
                continue;
            }

            OOFEM_LOG_DEBUG( " %ld (h)", local_hexa_id(tmp_hexa->ngb_elem_id [ j ]) );
        }

        OOFEM_LOG_DEBUG("\n");
    }

    OOFEM_LOG_DEBUG("\n");
#endif

    if ( mesh_edges != 0 ) {
        if ( ( mesh_edge_array = ( mesh_edge_rec * ) calloc( mesh_edges, sizeof( mesh_edge_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( mesh_faces != 0 ) {
        if ( ( mesh_face_array = ( mesh_face_rec * ) calloc( mesh_faces, sizeof( mesh_face_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( mesh_quads != 0 ) {
        if ( ( mesh_quad_array = ( mesh_quad_rec * ) calloc( mesh_quads, sizeof( mesh_quad_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    fine_nodes = mesh_edges + mesh_faces + mesh_quads + fe_faces + fe_quads + fe_tetras + fe_pyrams + fe_wedges + fe_hexas;

    fine_edges = fe_edges * 2;
    fine_quads = fe_faces * 3 + fe_quads * 4;
    fine_hexas = fe_tetras * 4 + fe_pyrams * 8 + fe_wedges * 8 + fe_hexas * 8;

    if ( fine_nodes != 0 ) {
        if ( ( fine_node_array = ( fe_node_rec * ) calloc( fine_nodes, sizeof( fe_node_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( fine_edges != 0 ) {
        if ( ( fine_edge_array = ( fine_edge_rec * ) calloc( fine_edges, sizeof( fine_edge_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( fine_quads != 0 ) {
        if ( ( fine_quad_array = ( fine_quad_rec * ) calloc( fine_quads, sizeof( fine_quad_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    if ( fine_hexas != 0 ) {
        if ( ( fine_hexa_array = ( fine_hexa_rec * ) calloc( fine_hexas, sizeof( fine_hexa_rec ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    fine_node_id = fe_nodes;
    fine_node = fine_node_array;

    mesh_edge_id = 0;
    mesh_edge = mesh_edge_array;

    if ( fe_edges != 0 ) {
        /* ensure that mesh_edges for fe_edges are created first */

        fe_edge = fe_edge_array;
        for ( i = 0; i < fe_edges; i++, fe_edge++ ) {
            node1 = ( nd1 = fe_edge->node [ 0 ] )->id;
            node2 = ( nd2 = fe_edge->node [ 1 ] )->id;

            fine_node->id = ++fine_node_id;

            /* note: do not swap nodes, keep tham in the same order as on fe_edge !!! */

            mesh_edge->node [ 0 ] = nd1;
            mesh_edge->node [ 1 ] = nd2;

            mesh_edge->mid_node = fine_node++;
            mesh_edge_id++;

            for ( k = node_num_nodes [ node1 - 1 ]; k < node_num_nodes [ node1 ]; k++ ) {
                if ( node_con_nodes [ k ] == node2 ) {
                    node_con_nodes [ k ] = -mesh_edge_id;
                    break;
                }
            }

            for ( k = node_num_nodes [ node2 - 1 ]; k < node_num_nodes [ node2 ]; k++ ) {
                if ( node_con_nodes [ k ] == node1 ) {
                    node_con_nodes [ k ] = -mesh_edge_id;
                    break;
                }
            }

            mesh_edge++;
        }
    }

#ifdef COMPLETE_DATA_STRUCTURE

    /* create mesh edges */

    fe_face = fe_face_array;
    tmp_face = tmp_face_array;
    for ( i = 0; i < fe_faces; i++, fe_face++, tmp_face++ ) {
        elem1 = global_face_id(i) + 1;                          /* macro */
        for ( j = 0; j < 3; j++ ) {
            elem2 = tmp_face->ngb_elem_id [ j ];

            if ( elem2 == 0 || elem1 < elem2 ) {
                j1 = face_ed_nd [ j ] [ 0 ];
                j2 = face_ed_nd [ j ] [ 1 ];
                node1 = ( nd1 = fe_face->node [ j1 ] )->id;
                node2 = ( nd2 = fe_face->node [ j2 ] )->id;

                fine_node->id = ++fine_node_id;

                /* note: ensure that starting node has smaller id */

                if ( node1 < node2 ) {
                    mesh_edge->node [ 0 ] = nd1;
                    mesh_edge->node [ 1 ] = nd2;
                } else {
                    mesh_edge->node [ 0 ] = nd2;
                    mesh_edge->node [ 1 ] = nd1;
                }

                mesh_edge->mid_node = fine_node++;
                mesh_edge_id++;

                for ( k = node_num_nodes [ node1 - 1 ]; k < node_num_nodes [ node1 ]; k++ ) {
                    node = node_con_nodes [ k ];
                    if ( node == node2 ) {
                        node_con_nodes [ k ] = -mesh_edge_id;
                        break;
                    }
                }

                for ( k = node_num_nodes [ node2 - 1 ]; k < node_num_nodes [ node2 ]; k++ ) {
                    node = node_con_nodes [ k ];
                    if ( node == node1 ) {
                        node_con_nodes [ k ] = -mesh_edge_id;
                        break;
                    }
                }

                tmp_face->edge [ j ] = mesh_edge;

                if ( elem2 != 0 ) {
                    if ( is_face(elem2) == 1 ) {
                        elem = local_face_id(elem2);      /* macro */
                        tmp_fc = & ( tmp_face_array [ elem - 1 ] );
                        for ( k = 0; k < 3; k++ ) {
                            if ( tmp_fc->ngb_elem_id [ k ] == elem1 ) {
                                tmp_fc->edge [ k ] = mesh_edge;
                                break;
                            }
                        }
                    }

                    if ( is_quad(elem2) == 1 ) {
                        elem = local_quad_id(elem2);      /* macro */
                        tmp_qd = & ( tmp_quad_array [ elem - 1 ] );
                        for ( k = 0; k < 4; k++ ) {
                            if ( tmp_qd->ngb_elem_id [ k ] == elem1 ) {
                                tmp_qd->edge [ k ] = mesh_edge;
                                break;
                            }
                        }
                    }
                }

                mesh_edge++;
            }
        }
    }

    fe_quad = fe_quad_array;
    tmp_quad = tmp_quad_array;
    for ( i = 0; i < fe_quads; i++, fe_quad++, tmp_quad++ ) {
        elem1 = global_quad_id(i) + 1;                          /* macro */
        for ( j = 0; j < 4; j++ ) {
            elem2 = tmp_quad->ngb_elem_id [ j ];

            if ( elem2 == 0 || elem1 < elem2 ) {
                j1 = quad_ed_nd [ j ] [ 0 ];
                j2 = quad_ed_nd [ j ] [ 1 ];
                node1 = ( nd1 = fe_quad->node [ j1 ] )->id;
                node2 = ( nd2 = fe_quad->node [ j2 ] )->id;

                fine_node->id = ++fine_node_id;

                /* note: ensure that starting node has smaller id */

                if ( node1 < node2 ) {
                    mesh_edge->node [ 0 ] = nd1;
                    mesh_edge->node [ 1 ] = nd2;
                } else {
                    mesh_edge->node [ 0 ] = nd2;
                    mesh_edge->node [ 1 ] = nd1;
                }

                mesh_edge->mid_node = fine_node++;
                mesh_edge_id++;

                for ( k = node_num_nodes [ node1 - 1 ]; k < node_num_nodes [ node1 ]; k++ ) {
                    node = node_con_nodes [ k ];
                    if ( node == node2 ) {
                        node_con_nodes [ k ] = -mesh_edge_id;
                        break;
                    }
                }

                for ( k = node_num_nodes [ node2 - 1 ]; k < node_num_nodes [ node2 ]; k++ ) {
                    node = node_con_nodes [ k ];
                    if ( node == node1 ) {
                        node_con_nodes [ k ] = -mesh_edge_id;
                        break;
                    }
                }

                tmp_quad->edge [ j ] = mesh_edge;

                if ( elem2 != 0 ) {
                    if ( is_face(elem2) == 1 ) {
                        elem = local_face_id(elem2);      /* macro */
                        tmp_fc = & ( tmp_face_array [ elem - 1 ] );
                        for ( k = 0; k < 3; k++ ) {
                            if ( tmp_fc->ngb_elem_id [ k ] == elem1 ) {
                                tmp_fc->edge [ k ] = mesh_edge;
                                break;
                            }
                        }
                    }

                    if ( is_quad(elem2) == 1 ) {
                        elem = local_quad_id(elem2);      /* macro */
                        tmp_qd = & ( tmp_quad_array [ elem - 1 ] );
                        for ( k = 0; k < 4; k++ ) {
                            if ( tmp_qd->ngb_elem_id [ k ] == elem1 ) {
                                tmp_qd->edge [ k ] = mesh_edge;
                                break;
                            }
                        }
                    }
                }

                mesh_edge++;
            }
        }
    }

#endif

    /* create mesh edges (only on remaining 3D mesh if COMPLETE_DATA_STRUCTURE is defined) */

    for ( i = 0; i < fe_nodes; i++ ) {
        node = i + 1;
        for ( j = node_num_nodes [ i ]; j < node_num_nodes [ node ]; j++ ) {
            /* note: negative nd means that there is already mesh_edge number */

            nd = node_con_nodes [ j ];
            if ( nd < 0 ) {
                node_con_nodes [ j ] = -nd;
                continue;
            }

#ifdef DEBUG
            if ( node > nd ) {
                error_message("Unexpected situation");
            }

#endif

            nd1 = & ( fe_node_array [ node - 1 ] );
            nd2 = & ( fe_node_array [ nd - 1 ] );

            fine_node->id = ++fine_node_id;

            /* note: nd1 should be smaller than nd2 */

            mesh_edge->node [ 0 ] = nd1;
            mesh_edge->node [ 1 ] = nd2;
            mesh_edge->mid_node = fine_node++;
            mesh_edge++;

            node_con_nodes [ j ] = ++mesh_edge_id;

            for ( k = node_num_nodes [ nd - 1 ]; k < node_num_nodes [ nd ]; k++ ) {
                if ( node_con_nodes [ k ] == node ) {
                    node_con_nodes [ k ] = -mesh_edge_id;
                    break;
                }
            }
        }
    }

    /* create mesh faces */

    mesh_face = mesh_face_array;

    fe_tetra = fe_tetra_array;
    tmp_tetra = tmp_tetra_array;
    for ( i = 0; i < fe_tetras; i++, fe_tetra++, tmp_tetra++ ) {
        elem1 = global_tetra_id(i) + 1;                          /* macro */
        for ( j = 0; j < 4; j++ ) {
            elem2 = tmp_tetra->ngb_elem_id [ j ];

            if ( elem2 == 0 || elem1 < elem2 ) {
                j1 = tetra_fc_nd [ j ] [ 0 ];
                j2 = tetra_fc_nd [ j ] [ 1 ];
                j3 = tetra_fc_nd [ j ] [ 2 ];
                nd1 = fe_tetra->node [ j1 ];
                nd2 = fe_tetra->node [ j2 ];
                nd3 = fe_tetra->node [ j3 ];

                fine_node->id = ++fine_node_id;

                mesh_face->node [ 0 ] = nd1;
                mesh_face->node [ 1 ] = nd2;
                mesh_face->node [ 2 ] = nd3;
                mesh_face->mid_node = fine_node++;

                tmp_tetra->face [ j ] = mesh_face;

                if ( elem2 != 0 ) {
                    elem = local_tetra_id(elem2);         /* macro */
                    tmp_tet = & ( tmp_tetra_array [ elem - 1 ] );
                    for ( k = 0; k < 4; k++ ) {
                        if ( tmp_tet->ngb_elem_id [ k ] == elem1 ) {
                            tmp_tet->face [ k ] = mesh_face;
                            break;
                        }
                    }
                }

                mesh_face++;
            }
        }
    }

    /* create mesh quads */

    mesh_quad = mesh_quad_array;

    fe_hexa = fe_hexa_array;
    tmp_hexa = tmp_hexa_array;
    for ( i = 0; i < fe_hexas; i++, fe_hexa++, tmp_hexa++ ) {
        elem1 = global_hexa_id(i) + 1;                          /* macro */
        for ( j = 0; j < 6; j++ ) {
            elem2 = tmp_hexa->ngb_elem_id [ j ];

            if ( elem2 == 0 || elem1 < elem2 ) {
                j1 = hexa_fc_nd [ j ] [ 0 ];
                j2 = hexa_fc_nd [ j ] [ 1 ];
                j3 = hexa_fc_nd [ j ] [ 2 ];
                j4 = hexa_fc_nd [ j ] [ 3 ];
                nd1 = fe_hexa->node [ j1 ];
                nd2 = fe_hexa->node [ j2 ];
                nd3 = fe_hexa->node [ j3 ];
                nd4 = fe_hexa->node [ j4 ];

                fine_node->id = ++fine_node_id;

                mesh_quad->node [ 0 ] = nd1;
                mesh_quad->node [ 1 ] = nd2;
                mesh_quad->node [ 2 ] = nd3;
                mesh_quad->node [ 3 ] = nd4;
                mesh_quad->mid_node = fine_node++;

                tmp_hexa->quad [ j ] = mesh_quad;

                if ( elem2 != 0 ) {
                    elem = local_hexa_id(elem2);         /* macro */
                    tmp_hx = & ( tmp_hexa_array [ elem - 1 ] );
                    for ( k = 0; k < 6; k++ ) {
                        if ( tmp_hx->ngb_elem_id [ k ] == elem1 ) {
                            tmp_hx->quad [ k ] = mesh_quad;
                            break;
                        }
                    }
                }

                mesh_quad++;
            }
        }
    }

    refine_nodes = 0;
    refine_nodes += mesh_edges * level * 2;
    refine_nodes += mesh_faces * ( level * 3 + level * level * 3 );
    refine_nodes += mesh_quads * ( level * 4 + level * level * 4 );
    refine_nodes += fe_faces * ( level * 3 + level * level * 3 );
    refine_nodes += fe_quads * ( level * 4 + level * level * 4 );
    refine_nodes += fe_tetras * ( level * 4 + level * level * 6 + level * level * level * 4 );
    refine_nodes += fe_hexas * ( level * 6 + level * level * 12 + level * level * level * 8 );

    refine_node_id = fe_nodes + fine_nodes;

    /* refine mesh edges */

    mesh_edge = mesh_edge_array;
    for ( i = 0; i < mesh_edges; i++, mesh_edge++ ) {
        for ( j = 0; j < 2; j++ ) {
            if ( ( mesh_edge->fine_id [ j ] = ( int * ) calloc( level + 2, sizeof( int ) ) ) == NULL ) {
                error_message("Memory allocation error");
            }

            mesh_edge->fine_id [ j ] [ 0 ] = mesh_edge->node [ j ]->id;
            for ( k = 1; k < level + 1; k++ ) {
                mesh_edge->fine_id [ j ] [ k ] = ++refine_node_id;
            }

            mesh_edge->fine_id [ j ] [ level + 1 ] = mesh_edge->mid_node->id;
        }
    }

    /*
     *    o-------o-------o                     o
     |       |       |                    / \
     |       |       |                   / \
     |       |       |                  / \
     *    o-------o-------o                 o       o
     |       |       |                /  \   / \
     |      cen      |               /     o \
     |       |       |              /      | \
     *    o-start-o--end--o             o-------o-------o
     */

    for ( j = 0; j < 4; j++ ) {
        if ( ( tmp_array_start [ j ] = ( long * ) calloc( level + 2, sizeof( long ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }

        if ( ( tmp_array_end [ j ] = ( long * ) calloc( level + 2, sizeof( long ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }

        if ( ( tmp_array_cen [ j ] = ( long * ) calloc( level + 2, sizeof( long ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    /* refine mesh faces */

    mesh_face = mesh_face_array;
    for ( i = 0; i < mesh_faces; i++, mesh_face++ ) {
        for ( j = 0; j < 3; j++ ) {
            j1 = face_ed_nd [ j ] [ 0 ];
            j2 = face_ed_nd [ j ] [ 1 ];

            node1 = mesh_face->node [ j1 ]->id;
            node2 = mesh_face->node [ j2 ]->id;

            /* I do rely on the fact that mesh_edge starts with node with smaller id */

            if ( node1 < node2 ) {
                for ( k = node_num_nodes [ node1 - 1 ]; k < node_num_nodes [ node1 ]; k++ ) {
                    mesh_edge = & ( mesh_edge_array [ node_con_nodes [ k ] - 1 ] );
                    if ( mesh_edge->node [ 1 ]->id == node2 ) {
                        for ( m = 0; m < level + 2; m++ ) {
                            tmp_array_start [ j ] [ m ] = mesh_edge->fine_id [ 0 ] [ m ];
                            tmp_array_end [ j ] [ m ] = mesh_edge->fine_id [ 1 ] [ m ];
                        }

                        break;
                    }
                }
            } else {
                for ( k = node_num_nodes [ node2 - 1 ]; k < node_num_nodes [ node2 ]; k++ ) {
                    mesh_edge = & ( mesh_edge_array [ node_con_nodes [ k ] - 1 ] );
                    if ( mesh_edge->node [ 1 ]->id == node1 ) {
                        for ( m = 0; m < level + 2; m++ ) {
                            tmp_array_start [ j ] [ m ] = mesh_edge->fine_id [ 1 ] [ m ];
                            tmp_array_end [ j ] [ m ] = mesh_edge->fine_id [ 0 ] [ m ];
                        }

                        break;
                    }
                }
            }

            tmp_array_cen [ j ] [ 0 ] = mesh_edge->mid_node->id;
            for ( m = 1; m < level + 1; m++ ) {
                tmp_array_cen [ j ] [ m ] = ++refine_node_id;
            }

            tmp_array_cen [ j ] [ level + 1 ] = mesh_face->mid_node->id;
        }

        for ( j = 0; j < 3; j++ ) {
            if ( ( mesh_face->fine_id [ j ] = ( int * ) calloc( ( level + 2 ) * ( level + 2 ), sizeof( int ) ) ) == NULL ) {
                error_message("Memory allocation error");
            }

            if ( ( jj = j - 1 ) < 0 ) {
                jj = 2;
            }

            j1 = face_ed_nd [ jj ] [ 0 ];
            j2 = face_ed_nd [ jj ] [ 1 ];

            pos = 0;
            for ( k = 0; k < level + 2; k++ ) {
                mesh_face->fine_id [ j ] [ pos++ ] = tmp_array_start [ j2 ] [ k ];
            }

            for ( k = 1; k < level + 1; k++ ) {
                mesh_face->fine_id [ j ] [ pos++ ] = tmp_array_end [ j1 ] [ k ];
                for ( m = 1; m < level + 1; m++ ) {
                    mesh_face->fine_id [ j ] [ pos++ ] = ++refine_node_id;
                }

                mesh_face->fine_id [ j ] [ pos++ ] = tmp_array_cen [ j2 ] [ k ];
            }

            for ( k = 0; k < level + 2; k++ ) {
                mesh_face->fine_id [ j ] [ pos++ ] = tmp_array_cen [ j1 ] [ k ];
            }
        }
    }

    /* refine mesh quads */

    mesh_quad = mesh_quad_array;
    for ( i = 0; i < mesh_quads; i++, mesh_quad++ ) {
        for ( j = 0; j < 4; j++ ) {
            j1 = quad_ed_nd [ j ] [ 0 ];
            j2 = quad_ed_nd [ j ] [ 1 ];

            node1 = mesh_quad->node [ j1 ]->id;
            node2 = mesh_quad->node [ j2 ]->id;

            /* I do rely on the fact that mesh_edge starts with node with smaller id */

            if ( node1 < node2 ) {
                for ( k = node_num_nodes [ node1 - 1 ]; k < node_num_nodes [ node1 ]; k++ ) {
                    mesh_edge = & ( mesh_edge_array [ node_con_nodes [ k ] - 1 ] );
                    if ( mesh_edge->node [ 1 ]->id == node2 ) {
                        for ( m = 0; m < level + 2; m++ ) {
                            tmp_array_start [ j ] [ m ] = mesh_edge->fine_id [ 0 ] [ m ];
                            tmp_array_end [ j ] [ m ] = mesh_edge->fine_id [ 1 ] [ m ];
                        }

                        break;
                    }
                }
            } else {
                for ( k = node_num_nodes [ node2 - 1 ]; k < node_num_nodes [ node2 ]; k++ ) {
                    mesh_edge = & ( mesh_edge_array [ node_con_nodes [ k ] - 1 ] );
                    if ( mesh_edge->node [ 1 ]->id == node1 ) {
                        for ( m = 0; m < level + 2; m++ ) {
                            tmp_array_start [ j ] [ m ] = mesh_edge->fine_id [ 1 ] [ m ];
                            tmp_array_end [ j ] [ m ] = mesh_edge->fine_id [ 0 ] [ m ];
                        }

                        break;
                    }
                }
            }

            tmp_array_cen [ j ] [ 0 ] = mesh_edge->mid_node->id;
            for ( m = 1; m < level + 1; m++ ) {
                tmp_array_cen [ j ] [ m ] = ++refine_node_id;
            }

            tmp_array_cen [ j ] [ level + 1 ] = mesh_quad->mid_node->id;
        }

        for ( j = 0; j < 4; j++ ) {
            if ( ( mesh_quad->fine_id [ j ] = ( int * ) calloc( ( level + 2 ) * ( level + 2 ), sizeof( int ) ) ) == NULL ) {
                error_message("Memory allocation error");
            }

            if ( ( jj = j - 1 ) < 0 ) {
                jj = 3;
            }

            j1 = quad_ed_nd [ jj ] [ 0 ];
            j2 = quad_ed_nd [ jj ] [ 1 ];

            pos = 0;
            for ( k = 0; k < level + 2; k++ ) {
                mesh_quad->fine_id [ j ] [ pos++ ] = tmp_array_start [ j2 ] [ k ];
            }

            for ( k = 1; k < level + 1; k++ ) {
                mesh_quad->fine_id [ j ] [ pos++ ] = tmp_array_end [ j1 ] [ k ];
                for ( m = 1; m < level + 1; m++ ) {
                    mesh_quad->fine_id [ j ] [ pos++ ] = ++refine_node_id;
                }

                mesh_quad->fine_id [ j ] [ pos++ ] = tmp_array_cen [ j2 ] [ k ];
            }

            for ( k = 0; k < level + 2; k++ ) {
                mesh_quad->fine_id [ j ] [ pos++ ] = tmp_array_cen [ j1 ] [ k ];
            }
        }
    }

#ifdef COMPLETE_DATA_STRUCTURE
    if ( node_num_nodes != NULL ) {
        free(node_num_nodes);
    }

    if ( node_con_nodes != NULL ) {
        free(node_con_nodes);
    }

#endif

    /* create fine edges */

    /* this is really not necessary because everything is stored in mesh_edge;
     * however, mesh_edge_array is going to be freed */

    fine_edge_id = 0;

    mesh_edge = mesh_edge_array;
    for ( i = 0; i < fe_edges; i++, mesh_edge++ ) {
        for ( j = 0; j < 2; j++ ) {
            fine_edge = & ( fine_edge_array [ fine_edge_id++ ] );
            if ( ( fine_edge->fine_id = ( int * ) calloc( ( level + 2 ), sizeof( int ) ) ) == NULL ) {
                error_message("Memory allocation error");
            }

            for ( k = 0; k < level + 2; k++ ) {
                fine_edge->fine_id [ k ] = mesh_edge->fine_id [ j ] [ k ];
            }
        }
    }

    /* create fine quads */

    fine_quad_id = 0;

    fe_face = fe_face_array;
    for ( i = 0; i < fe_faces; i++, fe_face++, fine_node++ ) {
        fine_node->id = ++fine_node_id;

#ifdef COMPLETE_DATA_STRUCTURE
        tmp_face = & tmp_face_array [ i ];
#endif

        for ( j = 0; j < 3; j++ ) {
            j1 = face_ed_nd [ j ] [ 0 ];
            j2 = face_ed_nd [ j ] [ 1 ];

            node1 = fe_face->node [ j1 ]->id;
            node2 = fe_face->node [ j2 ]->id;

#ifdef COMPLETE_DATA_STRUCTURE
            mesh_edge = tmp_face->edge [ j ];
#endif

            /* I do rely on the fact that mesh_edge starts with node with smaller id */

            if ( node1 < node2 ) {
#ifndef COMPLETE_DATA_STRUCTURE
                for ( k = node_num_nodes [ node1 - 1 ]; k < node_num_nodes [ node1 ]; k++ ) {
                    mesh_edge = & ( mesh_edge_array [ node_con_nodes [ k ] - 1 ] );
                    if ( mesh_edge->node [ 1 ]->id == node2 ) {
                        break;
                    }
                }

#endif

                for ( k = 0; k < level + 2; k++ ) {
                    tmp_array_start [ j ] [ k ] = mesh_edge->fine_id [ 0 ] [ k ];
                    tmp_array_end [ j ] [ k ] = mesh_edge->fine_id [ 1 ] [ k ];
                }
            } else {
#ifndef COMPLETE_DATA_STRUCTURE
                for ( k = node_num_nodes [ node2 - 1 ]; k < node_num_nodes [ node2 ]; k++ ) {
                    mesh_edge = & ( mesh_edge_array [ node_con_nodes [ k ] - 1 ] );
                    if ( mesh_edge->node [ 1 ]->id == node1 ) {
                        break;
                    }
                }

#endif

                for ( k = 0; k < level + 2; k++ ) {
                    tmp_array_start [ j ] [ k ] = mesh_edge->fine_id [ 1 ] [ k ];
                    tmp_array_end [ j ] [ k ] = mesh_edge->fine_id [ 0 ] [ k ];
                }
            }

            tmp_array_cen [ j ] [ 0 ] = mesh_edge->mid_node->id;
            for ( k = 1; k < level + 1; k++ ) {
                tmp_array_cen [ j ] [ k ] = ++refine_node_id;
            }

            tmp_array_cen [ j ] [ level + 1 ] = fine_node->id;
        }

        for ( j = 0; j < 3; j++ ) {
            fine_quad = & ( fine_quad_array [ fine_quad_id++ ] );
            if ( ( fine_quad->fine_id = ( int * ) calloc( ( level + 2 ) * ( level + 2 ), sizeof( int ) ) ) == NULL ) {
                error_message("Memory allocation error");
            }

            if ( ( jj = j - 1 ) < 0 ) {
                jj = 2;
            }

            j1 = face_ed_nd [ jj ] [ 0 ];
            j2 = face_ed_nd [ jj ] [ 1 ];

            pos = 0;
            for ( k = 0; k < level + 2; k++ ) {
                fine_quad->fine_id [ pos++ ] = tmp_array_start [ j2 ] [ k ];
            }

            for ( k = 1; k < level + 1; k++ ) {
                fine_quad->fine_id [ pos++ ] = tmp_array_end [ j1 ] [ k ];
                for ( m = 1; m < level + 1; m++ ) {
                    fine_quad->fine_id [ pos++ ] = ++refine_node_id;
                }

                fine_quad->fine_id [ pos++ ] = tmp_array_cen [ j2 ] [ k ];
            }

            for ( k = 0; k < level + 2; k++ ) {
                fine_quad->fine_id [ pos++ ] = tmp_array_cen [ j1 ] [ k ];
            }
        }
    }

    fe_quad = fe_quad_array;
    for ( i = 0; i < fe_quads; i++, fe_quad++, fine_node++ ) {
        fine_node->id = ++fine_node_id;

#ifdef COMPLETE_DATA_STRUCTURE
        tmp_quad = & ( tmp_quad_array [ i ] );
#endif

        for ( j = 0; j < 4; j++ ) {
            j1 = quad_ed_nd [ j ] [ 0 ];
            j2 = quad_ed_nd [ j ] [ 1 ];

            node1 = fe_quad->node [ j1 ]->id;
            node2 = fe_quad->node [ j2 ]->id;

#ifdef COMPLETE_DATA_STRUCTURE
            mesh_edge = tmp_quad->edge [ j ];
#endif

            /* I do rely on the fact that mesh_edge starts with node with smaller id */

            if ( node1 < node2 ) {
#ifndef COMPLETE_DATA_STRUCTURE
                for ( k = node_num_nodes [ node1 - 1 ]; k < node_num_nodes [ node1 ]; k++ ) {
                    mesh_edge = & ( mesh_edge_array [ node_con_nodes [ k ] - 1 ] );
                    if ( mesh_edge->node [ 1 ]->id == node2 ) {
                        break;
                    }
                }

#endif

                for ( k = 0; k < level + 2; k++ ) {
                    tmp_array_start [ j ] [ k ] = mesh_edge->fine_id [ 0 ] [ k ];
                    tmp_array_end [ j ] [ k ] = mesh_edge->fine_id [ 1 ] [ k ];
                }
            } else {
#ifndef COMPLETE_DATA_STRUCTURE
                for ( k = node_num_nodes [ node2 - 1 ]; k < node_num_nodes [ node2 ]; k++ ) {
                    mesh_edge = & ( mesh_edge_array [ node_con_nodes [ k ] - 1 ] );
                    if ( mesh_edge->node [ 1 ]->id == node1 ) {
                        break;
                    }
                }

#endif

                for ( k = 0; k < level + 2; k++ ) {
                    tmp_array_start [ j ] [ k ] = mesh_edge->fine_id [ 1 ] [ k ];
                    tmp_array_end [ j ] [ k ] = mesh_edge->fine_id [ 0 ] [ k ];
                }
            }

            tmp_array_cen [ j ] [ 0 ] = mesh_edge->mid_node->id;
            for ( k = 1; k < level + 1; k++ ) {
                tmp_array_cen [ j ] [ k ] = ++refine_node_id;
            }

            tmp_array_cen [ j ] [ level + 1 ] = fine_node->id;
        }

        for ( j = 0; j < 4; j++ ) {
            fine_quad = & ( fine_quad_array [ fine_quad_id++ ] );
            if ( ( fine_quad->fine_id = ( int * ) calloc( ( level + 2 ) * ( level + 2 ), sizeof( int ) ) ) == NULL ) {
                error_message("Memory allocation error");
            }

            if ( ( jj = j - 1 ) < 0 ) {
                jj = 3;
            }

            j1 = quad_ed_nd [ jj ] [ 0 ];
            j2 = quad_ed_nd [ jj ] [ 1 ];

            pos = 0;
            for ( k = 0; k < level + 2; k++ ) {
                fine_quad->fine_id [ pos++ ] = tmp_array_start [ j2 ] [ k ];
            }

            for ( k = 1; k < level + 1; k++ ) {
                fine_quad->fine_id [ pos++ ] = tmp_array_end [ j1 ] [ k ];
                for ( m = 1; m < level + 1; m++ ) {
                    fine_quad->fine_id [ pos++ ] = ++refine_node_id;
                }

                fine_quad->fine_id [ pos++ ] = tmp_array_cen [ j2 ] [ k ];
            }

            for ( k = 0; k < level + 2; k++ ) {
                fine_quad->fine_id [ pos++ ] = tmp_array_cen [ j1 ] [ k ];
            }
        }
    }

    for ( j = 0; j < 4; j++ ) {
        free(tmp_array_start [ j ]);
        free(tmp_array_end [ j ]);
        free(tmp_array_cen [ j ]);
    }

#ifdef COMPLETE_DATA_STRUCTURE
    if ( tmp_face_array != NULL ) {
        free(tmp_face_array);
    }

    if ( tmp_quad_array != NULL ) {
        free(tmp_quad_array);
    }

#else
    if ( node_num_nodes != NULL ) {
        free(node_num_nodes);
    }

    if ( node_con_nodes != NULL ) {
        free(node_con_nodes);
    }

#endif

    mesh_edge = mesh_edge_array;
    for ( i = 0; i < mesh_edges; i++, mesh_edge++ ) {
        for ( j = 0; j < 2; j++ ) {
            if ( mesh_edge->fine_id [ j ] != NULL ) {
                free(mesh_edge->fine_id [ j ]);
            }
        }
    }

    if ( mesh_edge_array != NULL ) {
        free(mesh_edge_array);
    }

    /* create fine hexas */

    for ( j = 0; j < 6; j++ ) {
        if ( ( tmp_array_cen [ j ] = ( long * ) calloc( level + 2, sizeof( long ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    for ( j = 0; j < 12; j++ ) {
        if ( ( tmp_array [ j ] = ( long * ) calloc( ( level + 2 ) * ( level + 2 ), sizeof( long ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }
    }

    fine_hexa_id = 0;

    fe_tetra = fe_tetra_array;
    tmp_tetra = tmp_tetra_array;
    for ( i = 0; i < fe_tetras; i++, fe_tetra++, tmp_tetra++, fine_node++ ) {
        fine_node->id = ++fine_node_id;

        for ( j = 0; j < 4; j++ ) {
            tmp_array_cen [ j ] [ 0 ] = tmp_tetra->face [ j ]->mid_node->id;
            for ( k = 1; k < level + 1; k++ ) {
                tmp_array_cen [ j ] [ k ] = ++refine_node_id;
            }

            tmp_array_cen [ j ] [ level + 1 ] = fine_node->id;
        }

        /* mark swapped faces (see ordering of nodes on tetra faces above);
         * bottom face is swapped if normal is outer;
         * side faces are swapped if normal is inner */

        for ( j = 0; j < 4; j++ ) {
            mesh_face = mesh_fc [ j ] = tmp_tetra->face [ j ];
            for ( k = 0; k < 3; k++ ) {
                if ( mesh_face->node [ k ] == fe_tetra->node [ tetra_fc_nd [ j ] [ 0 ] ] ) {
                    if ( ( kk = k + 1 ) == 3 ) {
                        kk = 0;
                    }

                    if ( mesh_face->node [ kk ] == fe_tetra->node [ tetra_fc_nd [ j ] [ 1 ] ] ) {
                        swap [ j ] = 0;
                    } else {
                        swap [ j ] = 1;
                    }

                    break;
                }
            }
        }

        /* form inner tmp fine quads "perpendicular" to bottom face;
         * each one is oriented u: from center of bottom edge to center of bottom face
         *                      v: from center of bottom edge to center of side face */

        for ( j = 0; j < 3; j++ ) {
            pos = 0;

            p = ( level + 2 ) * ( level + 1 );
            mesh_face = mesh_fc [ 0 ];
            for ( k = 0; k < 3; k++ ) {
                if ( mesh_face->node [ k ] == fe_tetra->node [ j ] ) {
                    if ( swap [ 0 ] == 1 ) {
                        kk = k;
                    } else {
                        if ( ( kk = k + 1 ) == 3 ) {
                            kk = 0;
                        }
                    }

                    break;
                }
            }

            for ( m = 0; m < level + 2; m++ ) {
                tmp_array [ j ] [ pos++ ] = mesh_face->fine_id [ kk ] [ p++ ];
            }

            if ( level != 0 ) {
                p = ( level + 2 ) * ( level + 1 );
                mesh_face = mesh_fc [ j + 1 ];
                for ( k = 0; k < 3; k++ ) {
                    if ( mesh_face->node [ k ] == fe_tetra->node [ j ] ) {
                        if ( swap [ j + 1 ] == 1 ) {
                            kk = k;
                        } else {
                            if ( ( kk = k + 1 ) == 3 ) {
                                kk = 0;
                            }
                        }

                        break;
                    }
                }

                for ( n = 1; n < level + 1; n++ ) {
                    tmp_array [ j ] [ pos++ ] = mesh_face->fine_id [ kk ] [ ++p ];
                    for ( m = 0; m < level; m++ ) {
                        tmp_array [ j ] [ pos++ ] = ++refine_node_id;
                    }

                    tmp_array [ j ] [ pos++ ] = tmp_array_cen [ 0 ] [ n ];
                }
            }

            for ( m = 0; m < level + 2; m++ ) {
                tmp_array [ j ] [ pos++ ] = tmp_array_cen [ j + 1 ] [ m ];
            }
        }

        /* form inner tmp fine quads "parallel" to bottom face;
         * each one is oriented as fine quads on bottom face if not swapped (this is 0-1-2) */

        for ( j = 0; j < 3; j++ ) {
            j1 = face_ed_nd [ j ] [ 0 ];
            j2 = face_ed_nd [ j ] [ 1 ];

            pos = 0;

            p = ( level + 2 ) * ( level + 1 );
            mesh_face = mesh_fc [ j + 1 ];
            for ( k = 0; k < 3; k++ ) {
                if ( mesh_face->node [ k ] == fe_tetra->node [ 3 ] ) {
                    if ( swap [ j + 1 ] == 1 ) {
                        kk = k;
                    } else {
                        if ( ( kk = k + 1 ) == 3 ) {
                            kk = 0;
                        }
                    }

                    break;
                }
            }

            for ( m = 0; m < level + 2; m++ ) {
                tmp_array [ j + 3 ] [ pos++ ] = mesh_face->fine_id [ kk ] [ p++ ];
            }

            if ( ( jj = j ) == 0 ) {
                jj = 3;
            }

            if ( level != 0 ) {
                p = ( level + 2 ) * ( level + 1 );
                mesh_face = mesh_fc [ jj ];
                for ( k = 0; k < 3; k++ ) {
                    if ( mesh_face->node [ k ] == fe_tetra->node [ j ] ) {
                        if ( swap [ jj ] == 1 ) {
                            kk = k;
                        } else {
                            if ( ( kk = k + 1 ) == 3 ) {
                                kk = 0;
                            }
                        }

                        break;
                    }
                }

                for ( n = 1; n < level + 1; n++ ) {
                    tmp_array [ j + 3 ] [ pos++ ] = mesh_face->fine_id [ kk ] [ ++p ];
                    for ( m = 0; m < level; m++ ) {
                        tmp_array [ j + 3 ] [ pos++ ] = ++refine_node_id;
                    }

                    tmp_array [ j + 3 ] [ pos++ ] = tmp_array_cen [ j + 1 ] [ n ];
                }
            }

            for ( m = 0; m < level + 2; m++ ) {
                tmp_array [ j + 3 ] [ pos++ ] = tmp_array_cen [ jj ] [ m ];
            }
        }

        /* form bottom fine hexas */

        /* filling 3d matrix in a sequential way would be too difficult
         * with respect to switching among boundary and inner sources;
         * therefore it is much easier to traverse boundary sources
         * and localize them into 3d matrix;
         * this is done using macros matrix_3d and possible matrix_2d;
         * when traversing whole boudary source it is better to use direct
         * access instead of macro matrix_2d */

        /* note: I ignore that some entries of 3d matrix are rewritten several times;
         * in order to enable direct access */

        for ( j = 0; j < 3; j++ ) {
            fine_hexa = & ( fine_hexa_array [ fine_hexa_id++ ] );
            if ( ( fine_hexa->fine_id = ( int * ) calloc( ( level + 2 ) * ( level + 2 ) * ( level + 2 ), sizeof( int ) ) ) == NULL ) {
                error_message("Memory allocation error");
            }

            mesh_face = mesh_fc [ 0 ];
            for ( k = 0; k < 3; k++ ) {
                if ( mesh_face->node [ k ] == fe_tetra->node [ j ] ) {
                    break;
                }
            }

            if ( swap [ 0 ] == 1 ) {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, n, m, 0) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            } else {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, m, n, 0) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            }

            mesh_face = mesh_fc [ j + 1 ];
            for ( k = 0; k < 3; k++ ) {
                if ( mesh_face->node [ k ] == fe_tetra->node [ j ] ) {
                    break;
                }
            }

            if ( swap [ j + 1 ] == 1 ) {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, n, 0, m) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            } else {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, m, 0, n) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            }

            if ( ( jj = j ) == 0 ) {
                jj = 3;
            }

            mesh_face = mesh_fc [ jj ];
            for ( k = 0; k < 3; k++ ) {
                if ( mesh_face->node [ k ] == fe_tetra->node [ j ] ) {
                    break;
                }
            }

            if ( swap [ jj ] == 1 ) {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, 0, m, n) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            } else {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, 0, n, m) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            }

            for ( pos = 0, n = 0; n < level + 2; n++ ) {
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    matrix_3d(fine_hexa->fine_id, level + 1, m, n) = tmp_array [ j ] [ pos ]; /* macro */
                    matrix_3d(fine_hexa->fine_id, m, level + 1, n) = tmp_array [ jj - 1 ] [ pos ]; /* macro */
                    matrix_3d(fine_hexa->fine_id, m, n, level + 1) = tmp_array [ j + 3 ] [ pos ]; /* macro */
                }
            }

            if ( level != 0 ) {
                for ( k = 1; k < level + 1; k++ ) {
                    for ( n = 1; n < level + 1; n++ ) {
                        for ( m = 1; m < level + 1; m++ ) {
                            matrix_3d(fine_hexa->fine_id, m, n, k) = ++refine_node_id;  /* macro */
                        }
                    }
                }
            }
        }

        /* form top hexa */

        fine_hexa = & ( fine_hexa_array [ fine_hexa_id++ ] );
        if ( ( fine_hexa->fine_id = ( int * ) calloc( ( level + 2 ) * ( level + 2 ) * ( level + 2 ), sizeof( int ) ) ) == NULL ) {
            error_message("Memory allocation error");
        }

        mesh_face = mesh_fc [ 1 ];
        for ( k = 0; k < 3; k++ ) {
            if ( mesh_face->node [ k ] == fe_tetra->node [ 3 ] ) {
                break;
            }
        }

        if ( swap [ 1 ] == 1 ) {
            for ( pos = 0, n = 0; n < level + 2; n++ ) {
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    matrix_3d(fine_hexa->fine_id, n, 0, m) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                }
            }
        } else {
            for ( pos = 0, n = 0; n < level + 2; n++ ) {
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    matrix_3d(fine_hexa->fine_id, m, 0, n) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                }
            }
        }

        mesh_face = mesh_fc [ 2 ];
        for ( k = 0; k < 3; k++ ) {
            if ( mesh_face->node [ k ] == fe_tetra->node [ 3 ] ) {
                break;
            }
        }

        if ( swap [ 2 ] == 1 ) {
            for ( pos = 0, n = 0; n < level + 2; n++ ) {
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    matrix_3d(fine_hexa->fine_id, 0, m, n) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                }
            }
        } else {
            for ( pos = 0, n = 0; n < level + 2; n++ ) {
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    matrix_3d(fine_hexa->fine_id, 0, n, m) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                }
            }
        }

        mesh_face = mesh_fc [ 3 ];
        for ( k = 0; k < 3; k++ ) {
            if ( mesh_face->node [ k ] == fe_tetra->node [ 3 ] ) {
                break;
            }
        }

        if ( swap [ 3 ] == 1 ) {
            for ( pos = 0, n = 0; n < level + 2; n++ ) {
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    matrix_3d(fine_hexa->fine_id, m, n, 0) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                }
            }
        } else {
            for ( pos = 0, n = 0; n < level + 2; n++ ) {
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    matrix_3d(fine_hexa->fine_id, n, m, 0) = mesh_face->fine_id [ k ] [ pos ]; /* macro */
                }
            }
        }

        for ( pos = 0, n = 0; n < level + 2; n++ ) {
            for ( m = 0; m < level + 2; m++, pos++ ) {
                matrix_3d(fine_hexa->fine_id, level + 1, n, m) = tmp_array [ 0 + 3 ] [ pos ]; /* macro */
                matrix_3d(fine_hexa->fine_id, m, level + 1, n) = tmp_array [ 2 + 3 ] [ pos ]; /* macro */
                matrix_3d(fine_hexa->fine_id, n, m, level + 1) = tmp_array [ 1 + 3 ] [ pos ]; /* macro */
            }
        }

        if ( level != 0 ) {
            for ( k = 1; k < level + 1; k++ ) {
                for ( n = 1; n < level + 1; n++ ) {
                    for ( m = 1; m < level + 1; m++ ) {
                        matrix_3d(fine_hexa->fine_id, m, n, k) = ++refine_node_id;    /* macro */
                    }
                }
            }
        }
    }

    fe_hexa = fe_hexa_array;
    tmp_hexa = tmp_hexa_array;
    for ( i = 0; i < fe_hexas; i++, fe_hexa++, tmp_hexa++, fine_node++ ) {
        fine_node->id = ++fine_node_id;

        for ( j = 0; j < 6; j++ ) {
            tmp_array_cen [ j ] [ 0 ] = tmp_hexa->quad [ j ]->mid_node->id;
            for ( k = 1; k < level + 1; k++ ) {
                tmp_array_cen [ j ] [ k ] = ++refine_node_id;
            }

            tmp_array_cen [ j ] [ level + 1 ] = fine_node->id;
        }

        /* mark swapped faces (see ordering of nodes on hexa faces above);
         * bottom and top faces are swapped if normal is outer;
         * side faces are swapped if normal is inner */

        for ( j = 0; j < 6; j++ ) {
            mesh_quad = mesh_qd [ j ] = tmp_hexa->quad [ j ];
            for ( k = 0; k < 4; k++ ) {
                if ( mesh_quad->node [ k ] == fe_hexa->node [ hexa_fc_nd [ j ] [ 0 ] ] ) {
                    if ( ( kk = k + 1 ) == 4 ) {
                        kk = 0;
                    }

                    if ( mesh_quad->node [ kk ] == fe_hexa->node [ hexa_fc_nd [ j ] [ 1 ] ] ) {
                        swap [ j ] = 0;
                    } else {
                        swap [ j ] = 1;
                    }

                    break;
                }
            }
        }

        /* form inner bottom tmp fine quads "perpendicular" to bottom face;
         * each one is oriented u: from center of bottom edge to center of bottom face
         *                      v: from center of bottom edge to center of side face */

        for ( j = 0; j < 4; j++ ) {
            pos = 0;

            p = ( level + 2 ) * ( level + 1 );
            mesh_quad = mesh_qd [ 0 ];
            for ( k = 0; k < 4; k++ ) {
                if ( mesh_quad->node [ k ] == fe_hexa->node [ j ] ) {
                    if ( swap [ 0 ] == 1 ) {
                        kk = k;
                    } else {
                        if ( ( kk = k + 1 ) == 4 ) {
                            kk = 0;
                        }
                    }

                    break;
                }
            }

            for ( m = 0; m < level + 2; m++ ) {
                tmp_array [ j ] [ pos++ ] = mesh_quad->fine_id [ kk ] [ p++ ];
            }

            if ( level != 0 ) {
                p = ( level + 2 ) * ( level + 1 );
                mesh_quad = mesh_qd [ j + 1 ];
                for ( k = 0; k < 4; k++ ) {
                    if ( mesh_quad->node [ k ] == fe_hexa->node [ j ] ) {
                        if ( swap [ j + 1 ] == 1 ) {
                            kk = k;
                        } else {
                            if ( ( kk = k + 1 ) == 4 ) {
                                kk = 0;
                            }
                        }

                        break;
                    }
                }

                for ( n = 1; n < level + 1; n++ ) {
                    tmp_array [ j ] [ pos++ ] = mesh_quad->fine_id [ kk ] [ ++p ];
                    for ( m = 0; m < level; m++ ) {
                        tmp_array [ j ] [ pos++ ] = ++refine_node_id;
                    }

                    tmp_array [ j ] [ pos++ ] = tmp_array_cen [ 0 ] [ n ];
                }
            }

            for ( m = 0; m < level + 2; m++ ) {
                tmp_array [ j ] [ pos++ ] = tmp_array_cen [ j + 1 ] [ m ];
            }
        }

        /* form top tmp fine quads "perpendicular" to top face;
         * each one is oriented u: from center of top edge to center of top face
         *                      v: from center of top edge to center of side face */

        for ( j = 0; j < 4; j++ ) {
            pos = 0;

            p = ( level + 2 ) * ( level + 1 );
            mesh_quad = mesh_qd [ 5 ];
            for ( k = 0; k < 4; k++ ) {
                if ( mesh_quad->node [ k ] == fe_hexa->node [ j + 4 ] ) {
                    if ( swap [ 5 ] == 0 ) {
                        kk = k;
                    } else {
                        if ( ( kk = k + 1 ) == 4 ) {
                            kk = 0;
                        }
                    }

                    break;
                }
            }

            for ( m = 0; m < level + 2; m++ ) {
                tmp_array [ j + 4 ] [ pos++ ] = mesh_quad->fine_id [ kk ] [ p++ ];
            }

            if ( level != 0 ) {
                p = ( level + 2 ) * ( level + 1 );
                mesh_quad = mesh_qd [ j + 1 ];
                for ( k = 0; k < 4; k++ ) {
                    if ( mesh_quad->node [ k ] == fe_hexa->node [ j + 4 ] ) {
                        if ( swap [ j + 1 ] == 0 ) {
                            kk = k;
                        } else {
                            if ( ( kk = k + 1 ) == 4 ) {
                                kk = 0;
                            }
                        }

                        break;
                    }
                }

                for ( n = 1; n < level + 1; n++ ) {
                    tmp_array [ j + 4 ] [ pos++ ] = mesh_quad->fine_id [ kk ] [ ++p ];
                    for ( m = 0; m < level; m++ ) {
                        tmp_array [ j + 4 ] [ pos++ ] = ++refine_node_id;
                    }

                    tmp_array [ j + 4 ] [ pos++ ] = tmp_array_cen [ 5 ] [ n ];
                }
            }

            for ( m = 0; m < level + 2; m++ ) {
                tmp_array [ j + 4 ] [ pos++ ] = tmp_array_cen [ j + 1 ] [ m ];
            }
        }

        /* form inner tmp fine quads "parallel" to bottom and top face;
         * each one is oriented as fine quads on bottom face if not swapped (this is 0-1-2-3) */

        for ( j = 0; j < 4; j++ ) {
            pos = 0;

            p = ( level + 2 ) * ( level + 1 );
            mesh_quad = mesh_qd [ j + 1 ];
            for ( k = 0; k < 4; k++ ) {
                if ( mesh_quad->node [ k ] == fe_hexa->node [ j + 4 ] ) {
                    if ( swap [ j + 1 ] == 1 ) {
                        kk = k;
                    } else {
                        if ( ( kk = k + 1 ) == 4 ) {
                            kk = 0;
                        }
                    }

                    break;
                }
            }

            for ( m = 0; m < level + 2; m++ ) {
                tmp_array [ j + 8 ] [ pos++ ] = mesh_quad->fine_id [ kk ] [ p++ ];
            }

            if ( ( jj = j ) == 0 ) {
                jj = 4;
            }

            if ( level != 0 ) {
                p = ( level + 2 ) * ( level + 1 );
                mesh_quad = mesh_qd [ jj ];
                for ( k = 0; k < 4; k++ ) {
                    if ( mesh_quad->node [ k ] == fe_hexa->node [ j ] ) {
                        if ( swap [ jj ] == 1 ) {
                            kk = k;
                        } else {
                            if ( ( kk = k + 1 ) == 4 ) {
                                kk = 0;
                            }
                        }

                        break;
                    }
                }

                for ( n = 1; n < level + 1; n++ ) {
                    tmp_array [ j + 8 ] [ pos++ ] = mesh_quad->fine_id [ kk ] [ ++p ];
                    for ( m = 0; m < level; m++ ) {
                        tmp_array [ j + 8 ] [ pos++ ] = ++refine_node_id;
                    }

                    tmp_array [ j + 8 ] [ pos++ ] = tmp_array_cen [ j + 1 ] [ n ];
                }
            }

            for ( m = 0; m < level + 2; m++ ) {
                tmp_array [ j + 8 ] [ pos++ ] = tmp_array_cen [ jj ] [ m ];
            }
        }

        /* form bottom fine hexas */

        /* filling 3d matrix in a sequential way would be too difficult
         * with respect to switching among boundary and inner sources;
         * therefore it is much easier to traverse boundary sources
         * and localize them into 3d matrix;
         * this is done using macros matrix_3d and possible matrix_2d;
         * when traversing whole boudary source it is better to use direct
         * access instead of macro matrix_2d */

        /* note: I ignore that some entries of 3d matrix are rewritten several times;
         * in order to enable direct access */

        for ( j = 0; j < 4; j++ ) {
            fine_hexa = & ( fine_hexa_array [ fine_hexa_id++ ] );
            if ( ( fine_hexa->fine_id = ( int * ) calloc( ( level + 2 ) * ( level + 2 ) * ( level + 2 ), sizeof( int ) ) ) == NULL ) {
                error_message("Memory allocation error");
            }

            mesh_quad = mesh_qd [ 0 ];
            for ( k = 0; k < 4; k++ ) {
                if ( mesh_quad->node [ k ] == fe_hexa->node [ j ] ) {
                    break;
                }
            }

            if ( swap [ 0 ] == 1 ) {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, n, m, 0) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            } else {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, m, n, 0) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            }

            mesh_quad = mesh_qd [ j + 1 ];
            for ( k = 0; k < 4; k++ ) {
                if ( mesh_quad->node [ k ] == fe_hexa->node [ j ] ) {
                    break;
                }
            }

            if ( swap [ j + 1 ] == 1 ) {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, n, 0, m) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            } else {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, m, 0, n) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            }

            if ( ( jj = j ) == 0 ) {
                jj = 4;
            }

            mesh_quad = mesh_qd [ jj ];
            for ( k = 0; k < 4; k++ ) {
                if ( mesh_quad->node [ k ] == fe_hexa->node [ j ] ) {
                    break;
                }
            }

            if ( swap [ jj ] == 1 ) {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, 0, m, n) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            } else {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, 0, n, m) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            }

            for ( pos = 0, n = 0; n < level + 2; n++ ) {
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    matrix_3d(fine_hexa->fine_id, level + 1, m, n) = tmp_array [ j ] [ pos ]; /* macro */
                    matrix_3d(fine_hexa->fine_id, m, level + 1, n) = tmp_array [ jj - 1 ] [ pos ]; /* macro */
                    matrix_3d(fine_hexa->fine_id, m, n, level + 1) = tmp_array [ j + 8 ] [ pos ]; /* macro */
                }
            }

            if ( level != 0 ) {
                for ( k = 1; k < level + 1; k++ ) {
                    for ( n = 1; n < level + 1; n++ ) {
                        for ( m = 1; m < level + 1; m++ ) {
                            matrix_3d(fine_hexa->fine_id, m, n, k) = ++refine_node_id;  /* macro */
                        }
                    }
                }
            }
        }

        /* form top fine hexas */

        for ( j = 0; j < 4; j++ ) {
            fine_hexa = & ( fine_hexa_array [ fine_hexa_id++ ] );
            if ( ( fine_hexa->fine_id = ( int * ) calloc( ( level + 2 ) * ( level + 2 ) * ( level + 2 ), sizeof( int ) ) ) == NULL ) {
                error_message("Memory allocation error");
            }

            mesh_quad = mesh_qd [ 5 ];
            for ( k = 0; k < 4; k++ ) {
                if ( mesh_quad->node [ k ] == fe_hexa->node [ j + 4 ] ) {
                    break;
                }
            }

            if ( swap [ 5 ] == 1 ) {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, n, m, 0) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            } else {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, m, n, 0) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            }

            mesh_quad = mesh_qd [ j + 1 ];
            for ( k = 0; k < 4; k++ ) {
                if ( mesh_quad->node [ k ] == fe_hexa->node [ j + 4 ] ) {
                    break;
                }
            }

            if ( swap [ j + 1 ] == 1 ) {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, 0, m, n) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            } else {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, 0, n, m) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            }

            if ( ( jj = j ) == 0 ) {
                jj = 4;
            }

            mesh_quad = mesh_qd [ jj ];
            for ( k = 0; k < 4; k++ ) {
                if ( mesh_quad->node [ k ] == fe_hexa->node [ j + 4 ] ) {
                    break;
                }
            }

            if ( swap [ jj ] == 1 ) {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, n, 0, m) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            } else {
                for ( pos = 0, n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        matrix_3d(fine_hexa->fine_id, m, 0, n) = mesh_quad->fine_id [ k ] [ pos ]; /* macro */
                    }
                }
            }

            for ( pos = 0, n = 0; n < level + 2; n++ ) {
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    matrix_3d(fine_hexa->fine_id, level + 1, m, n) = tmp_array [ jj - 1 + 4 ] [ pos ]; /* macro */
                    matrix_3d(fine_hexa->fine_id, m, level + 1, n) = tmp_array [ j + 4 ] [ pos ]; /* macro */
                    matrix_3d(fine_hexa->fine_id, n, m, level + 1) = tmp_array [ j + 8 ] [ pos ]; /* macro */
                }
            }

            if ( level != 0 ) {
                for ( k = 1; k < level + 1; k++ ) {
                    for ( n = 1; n < level + 1; n++ ) {
                        for ( m = 1; m < level + 1; m++ ) {
                            matrix_3d(fine_hexa->fine_id, m, n, k) = ++refine_node_id;  /* macro */
                        }
                    }
                }
            }
        }
    }

    for ( j = 0; j < 6; j++ ) {
        free(tmp_array_cen [ j ]);
        free(tmp_array [ j ]);
        free(tmp_array [ j + 6 ]);
    }

    if ( fe_node_array != NULL ) {
        free(fe_node_array);
    }

    if ( fe_edge_array != NULL ) {
        free(fe_edge_array);
    }

    if ( fe_face_array != NULL ) {
        free(fe_face_array);
    }

    if ( fe_quad_array != NULL ) {
        free(fe_quad_array);
    }

    if ( fe_tetra_array != NULL ) {
        free(fe_tetra_array);
    }

    if ( fe_hexa_array != NULL ) {
        free(fe_hexa_array);
    }

    mesh_face = mesh_face_array;
    for ( i = 0; i < mesh_faces; i++, mesh_face++ ) {
        for ( j = 0; j < 3; j++ ) {
            if ( mesh_face->fine_id [ j ] != NULL ) {
                free(mesh_face->fine_id [ j ]);
            }
        }
    }

    mesh_quad = mesh_quad_array;
    for ( i = 0; i < mesh_quads; i++, mesh_quad++ ) {
        for ( j = 0; j < 4; j++ ) {
            if ( mesh_quad->fine_id [ j ] != NULL ) {
                free(mesh_quad->fine_id [ j ]);
            }
        }
    }

    if ( mesh_face_array != NULL ) {
        free(mesh_face_array);
    }

    if ( mesh_quad_array != NULL ) {
        free(mesh_quad_array);
    }

    /*
     * if(refine_nodes != 0){
     * if((refine_node_array = (fe_node_rec *)calloc(refine_nodes, sizeof(fe_node_rec))) == NULL)
     * error_message("Memory allocation error");
     *
     * fine_edge = fine_edge_array;
     * for(i = 0; i < fine_edges; i++, fine_edge++){
     * pos = 0;
     * for(m = 0; m < level + 2; m++, pos++){
     *  node_id = fine_edge -> fine_id[pos];
     *  if(node_id <= fe_nodes + fine_nodes)continue;
     *
     *  refine_node = &(refine_node_array[node_id - fe_nodes - fine_nodes - 1]);
     *  refine_node -> id = node_id;
     * }
     * }
     *
     * fine_quad = fine_quad_array;
     * for(i = 0; i < fine_quads; i++, fine_quad++){
     * pos = 0;
     * for(n = 0; n < level + 2; n++){
     *  for(m = 0; m < level + 2; m++, pos++){
     *   node_id = fine_quad -> fine_id[pos];
     *   if(node_id <= fe_nodes + fine_nodes)continue;
     *
     *   refine_node = &(refine_node_array[node_id - fe_nodes - fine_nodes - 1]);
     *   refine_node -> id = node_id;
     *  }
     * }
     * }
     *
     * fine_hexa = fine_hexa_array;
     * for(i = 0; i < fine_hexas; i++, fine_hexa++){
     * pos = 0;
     * for(k = 0; k < level + 2; k++){
     *  for(n = 0; n < level + 2; n++){
     *   for(m = 0; m < level + 2; m++, pos++){
     *    node_id = fine_hexa -> fine_id[pos];
     *    if(node_id <= fe_nodes + fine_nodes)continue;
     *
     *    refine_node = &(refine_node_array[node_id - fe_nodes - fine_nodes - 1]);
     *    refine_node -> id = node_id;
     *   }
     *  }
     * }
     * }
     * }
     */

    if ( fine_node_array != NULL ) {
        free(fine_node_array);
    }

    /*
     * if(refine_node_array != NULL)free(refine_node_array);
     */

    edge_id = face_id = quad_id = tetra_id = hexa_id = 0;
    for ( i = 0; i < fe_elems; i++ ) {
        element = d->giveElement(i + 1);
        RefinedElement &refinedElement = refinedElementList.at(i + 1);
        boundary = refinedElement.giveBoundaryFlagArray();

        switch ( element->giveGeometryType() ) {
        case EGT_line_1:
            pos = edge_id * 2;
            fine_edge = & ( fine_edge_array [ pos ] );
            for ( j = 0; j < 2; j++, fine_edge++ ) {
                connectivity = refinedElement.giveFineNodeArray(j + 1);
                connectivity->resize(level + 2);
                for ( k = 0; k < level + 2; k++ ) {
                    connectivity->at(k + 1) = fine_edge->fine_id [ k ];
                }
            }

            for ( j = 0; j < 2; j++ ) {
                /* note: I cannot use fe_edge (it is already released)
                 *
                 *  node = fe_edge -> node[j] -> id;
                 */

                node = ( & ( fine_edge_array [ pos + j ] ) )->fine_id [ 0 ];
                if ( node_num_elems [ node ] - node_num_elems [ node - 1 ] == 1 ) {
                    boundary->at(j + 1) = 1;
                } else {
                    boundary->at(j + 1) = 0;
                }
            }

            edge_id++;
            break;

        case EGT_triangle_1:
            pos = face_id * 3;
            fine_quad = & ( fine_quad_array [ pos ] );
            for ( j = 0; j < 3; j++, fine_quad++ ) {
                connectivity = refinedElement.giveFineNodeArray(j + 1);
                connectivity->resize( ( level + 2 ) * ( level + 2 ) );
                for ( k = 0; k < ( level + 2 ) * ( level + 2 ); k++ ) {
                    connectivity->at(k + 1) = fine_quad->fine_id [ k ];
                }
            }

            /* note: I am not using tmp_face (and its ngb_elem_id) because tmp_face exists
             * only if COMPLETE_DATA_STRUCTURE is required but that is supported for 2-manifolds only */

            elem1 = global_face_id(face_id) + 1;
            for ( j = 0; j < 3; j++ ) {
                boundary->at(j + 1) = 1;

                j1 = face_ed_nd [ j ] [ 0 ];
                j2 = face_ed_nd [ j ] [ 1 ];

                /* note: I cannot use fe_face (it is already released)
                 *
                 *  node1 = fe_face -> node[j1] -> id;
                 *  node2 = fe_face -> node[j2] -> id;
                 */

                node1 = ( & ( fine_quad_array [ pos + j1 ] ) )->fine_id [ 0 ];
                node2 = ( & ( fine_quad_array [ pos + j2 ] ) )->fine_id [ 0 ];

                flag = 1;
                for ( k = node_num_elems [ node1 - 1 ]; k < node_num_elems [ node1 ]; k++ ) {
                    elem2 = node_con_elems [ k ];
                    if ( elem2 == elem1 ) {
                        continue;
                    }

                    if ( is_edge(elem2) == 0 ) {                                /* macro */
                        for ( m = node_num_elems [ node2 - 1 ]; m < node_num_elems [ node2 ]; m++ ) {
                            if ( elem2 == node_con_elems [ m ] ) {
                                boundary->at(j + 1) = flag = 0;
                                break;
                            }
                        }

                        if ( flag == 0 ) {
                            break;
                        }
                    }
                }
            }

            face_id++;
            break;

        case EGT_quad_1:
            pos = fe_faces * 3 + quad_id * 4;
            fine_quad = & ( fine_quad_array [ pos ] );
            for ( j = 0; j < 4; j++, fine_quad++ ) {
                connectivity = refinedElement.giveFineNodeArray(j + 1);
                connectivity->resize( ( level + 2 ) * ( level + 2 ) );
                for ( k = 0; k < ( level + 2 ) * ( level + 2 ); k++ ) {
                    connectivity->at(k + 1) = fine_quad->fine_id [ k ];
                }
            }

            /* note: I am not using tmp_quad (and its ngb_elem_id) because tmp_quad exists
             * only if COMPLETE_DATA_STRUCTURE is required but that is supported for 2-manifolds only */

            elem1 = global_quad_id(quad_id) + 1;
            for ( j = 0; j < 4; j++ ) {
                boundary->at(j + 1) = 1;

                j1 = quad_ed_nd [ j ] [ 0 ];
                j2 = quad_ed_nd [ j ] [ 1 ];

                /* note: I cannot use fe_quad (it is already released)
                 *
                 *  node1 = fe_quad -> node[j1] -> id;
                 *  node2 = fe_quad -> node[j2] -> id;
                 */

                node1 = ( & ( fine_quad_array [ pos + j1 ] ) )->fine_id [ 0 ];
                node2 = ( & ( fine_quad_array [ pos + j2 ] ) )->fine_id [ 0 ];

                flag = 1;
                for ( k = node_num_elems [ node1 - 1 ]; k < node_num_elems [ node1 ]; k++ ) {
                    elem2 = node_con_elems [ k ];
                    if ( elem2 == elem1 ) {
                        continue;
                    }

                    if ( is_edge(elem2) == 0 ) {                                /* macro */
                        for ( m = node_num_elems [ node2 - 1 ]; m < node_num_elems [ node2 ]; m++ ) {
                            if ( elem2 == node_con_elems [ m ] ) {
                                boundary->at(j + 1) = flag = 0;
                                break;
                            }
                        }

                        if ( flag == 0 ) {
                            break;
                        }
                    }
                }
            }

            quad_id++;
            break;

        case EGT_tetra_1:
            pos = tetra_id * 4;
            fine_hexa = & ( fine_hexa_array [ pos ] );
            for ( j = 0; j < 4; j++, fine_hexa++ ) {
                connectivity = refinedElement.giveFineNodeArray(j + 1);
                connectivity->resize( ( level + 2 ) * ( level + 2 ) * ( level + 2 ) );
                for ( k = 0; k < ( level + 2 ) * ( level + 2 ) * ( level + 2 ); k++ ) {
                    connectivity->at(k + 1) = fine_hexa->fine_id [ k ];
                }
            }

            tmp_tetra = & ( tmp_tetra_array [ tetra_id ] );
            for ( j = 0; j < 4; j++ ) {
                if ( tmp_tetra->ngb_elem_id [ j ] != 0 ) {
                    boundary->at(j + 1) = 0;
                } else {
                    boundary->at(j + 1) = 1;
                }
            }

            tetra_id++;
            break;

        case EGT_hexa_1:
            pos = fe_tetras * 4 + hexa_id * 8;
            fine_hexa = & ( fine_hexa_array [ pos ] );
            for ( j = 0; j < 8; j++, fine_hexa++ ) {
                connectivity = refinedElement.giveFineNodeArray(j + 1);
                connectivity->resize( ( level + 2 ) * ( level + 2 ) * ( level + 2 ) );
                for ( k = 0; k < ( level + 2 ) * ( level + 2 ) * ( level + 2 ); k++ ) {
                    connectivity->at(k + 1) = fine_hexa->fine_id [ k ];
                }
            }

            tmp_hexa = & ( tmp_hexa_array [ hexa_id ] );
            for ( j = 0; j < 6; j++ ) {
                if ( tmp_hexa->ngb_elem_id [ j ] != 0 ) {
                    boundary->at(j + 1) = 0;
                } else {
                    boundary->at(j + 1) = 1;
                }
            }

            hexa_id++;
            break;

        default:
            error_message("Not supported element type");
            break;
        }
    }

    /* the following 4 arrays cannot be released earlier because they are used to recover boundary flag */

    if ( node_num_elems != NULL ) {
        free(node_num_elems);
    }

    if ( node_con_elems != NULL ) {
        free(node_con_elems);
    }

    if ( tmp_tetra_array != NULL ) {
        free(tmp_tetra_array);
    }

    if ( tmp_hexa_array != NULL ) {
        free(tmp_hexa_array);
    }

    fine_edge = fine_edge_array;
    for ( i = 0; i < fine_edges; i++, fine_edge++ ) {
        if ( fine_edge->fine_id != NULL ) {
            free(fine_edge->fine_id);
        }
    }

    fine_quad = fine_quad_array;
    for ( i = 0; i < fine_quads; i++, fine_quad++ ) {
        if ( fine_quad->fine_id != NULL ) {
            free(fine_quad->fine_id);
        }
    }

    fine_hexa = fine_hexa_array;
    for ( i = 0; i < fine_hexas; i++, fine_hexa++ ) {
        if ( fine_hexa->fine_id != NULL ) {
            free(fine_hexa->fine_id);
        }
    }

    if ( fine_edge_array != NULL ) {
        free(fine_edge_array);
    }

    if ( fine_quad_array != NULL ) {
        free(fine_quad_array);
    }

    if ( fine_hexa_array != NULL ) {
        free(fine_hexa_array);
    }

    this->nodes = fe_nodes + fine_nodes + refine_nodes;
    this->edges = fine_edges * ( level + 1 );
    this->quads = fine_quads * ( level + 1 ) * ( level + 1 );
    this->hexas = fine_hexas * ( level + 1 ) * ( level + 1 ) * ( level + 1 );
    this->elems = this->edges + this->quads + this->hexas;

    return ( 0 );
}
} // end namespace oofem
