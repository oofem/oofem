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

#ifndef refinedmesh_h
#define refinedmesh_h

#include <vector>
#include <cstdlib>

namespace oofem {
#define SAVE                       /* memory savings */
//#define COMPLETE_DATA_STRUCTURE    /* only for 2-manifolds */

#ifdef SAVE
 #undef COMPLETE_DATA_STRUCTURE
#endif

class Domain;
class RefinedElement;

class RefinedMesh
{
public:
    int completed;
    int nodes, elems;
    int edges, quads, hexas;

private:
    typedef struct fine_edge_rec {
        int *fine_id;
    } fine_edge_rec;

    typedef struct fine_quad_rec {
        int *fine_id;
    } fine_quad_rec;

    typedef struct fine_hexa_rec {
        int *fine_id;
    } fine_hexa_rec;


    typedef struct fe_node_rec {
        int id;
    } fe_node_rec;

    typedef struct fe_edge_rec {
        fe_node_rec *node [ 2 ];
    } fe_edge_rec;

    typedef struct fe_face_rec {
        fe_node_rec *node [ 3 ];
    } fe_face_rec;

    typedef struct fe_quad_rec {
        fe_node_rec *node [ 4 ];
    } fe_quad_rec;

    typedef struct fe_tetra_rec {
        fe_node_rec *node [ 4 ];
    } fe_tetra_rec;

    typedef struct fe_hexa_rec {
        fe_node_rec *node [ 8 ];
    } fe_hexa_rec;



    typedef struct mesh_edge_rec {
        fe_node_rec *node [ 2 ];
        fe_node_rec *mid_node;
        int *fine_id [ 2 ];
    } mesh_edge_rec;

    typedef struct mesh_face_rec {
        fe_node_rec *node [ 3 ];
        fe_node_rec *mid_node;
        int *fine_id [ 3 ];
    } mesh_face_rec;

    typedef struct mesh_quad_rec {
        fe_node_rec *node [ 4 ];
        fe_node_rec *mid_node;
        int *fine_id [ 4 ];
    } mesh_quad_rec;



#ifdef COMPLETE_DATA_STRUCTURE
    typedef struct tmp_face_rec {
        int ngb_elem_id [ 3 ];
        mesh_edge_rec *edge [ 3 ];
    } tmp_face_rec;

    typedef struct tmp_quad_rec {
        int ngb_elem_id [ 4 ];
        mesh_edge_rec *edge [ 4 ];
    } tmp_quad_rec;
#endif


    typedef struct tmp_tetra_rec {
        int ngb_elem_id [ 4 ];
        mesh_face_rec *face [ 4 ];
    } tmp_tetra_rec;

    typedef struct tmp_hexa_rec {
        int ngb_elem_id [ 6 ];
        mesh_quad_rec *quad [ 6 ];
    } tmp_hexa_rec;


public:
    RefinedMesh() {
        completed = 0;
        nodes = edges = quads = hexas = 0;
    };

    ~RefinedMesh() { }

    int refineMeshGlobally(Domain *d, int level, std :: vector< RefinedElement > &refinedElementList);
};
} // end namespace oofem
#endif // refinedmesh_h
