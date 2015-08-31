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

//   *******************************
//   *** DELAUNAY MESH GENERATOR ***
//   *******************************

// based on oofemlib/mesherinterface.h


#ifndef delaunaytrinagulator_h
#define delaunaytrinagulator_h

#include <list>
#include "contextioresulttype.h"
#include "timer.h"
#include "octreelocalizert.h"
#include "intarray.h"
#include "edge2d.h"

namespace oofem {
class Domain;
class TimeStep;
class IntArray;
class FloatArray;
class InsertionData;

/**
 * Mesh generator for the PFEM problem, using Bowyer-Watson algorithm of the Delaunay triangulation
 * of a set of nodes (PFEMParticle) creating TR1_2D_PFEM elements.
 *
 * @author David Krybus
 */
class DelaunayTriangulator
{
protected:
    /// Domain of the PFEM problem containing nodes to be triangulated
    Domain *domain;
    // Value of alpha for the boundary recognition via alpha shape algorithm
    double alphaValue;
    // Number of nodes
    int nnode;
    // Option 2: setting bounds for computed Alpha
    //double minAlpha;
    //double maxAlpha;

    /// Measures overall time of triangulation procedure
    Timer meshingTimer;
    /// Measures time needed by searching for non-delaunay triangles
    Timer searchingTimer;
    /// Measures time needed for identifying polygon to be retriangulated
    Timer polygonTimer;
    /// Measures time needed for creating new Delaunay triangles
    Timer creativeTimer;
    /// Contains all triangles (even not valid)
    std :: list< DelaunayTriangle * >generalTriangleList;
    std :: list< DelaunayTriangle * > :: iterator genIT;

    /// Contains resulting alpha-shape in form of a list
    std :: list< AlphaEdge2D * >alphaShapeEdgeList;

    /// contains all edges of the triangulation
    std :: list< AlphaEdge2D * >edgeList;
    std :: list< AlphaEdge2D * > :: iterator elIT;

    /// Octree with Delaunay triangles allowing fast search
    OctreeSpatialLocalizerT< DelaunayTriangle * >triangleOctree;

public:
    /// Constructor
    DelaunayTriangulator(Domain *d, double setAlpha);
    /// Destructor
    ~DelaunayTriangulator();

    /// Main call
    void generateMesh();

private:
    /// Edge is added to the polygon only if it's not contained. Otherwise both are removed (edge shared by two non-Delaunay triangles).
    void addUniqueEdgeToPolygon(Edge2D *edge, std :: list< Edge2D > &polygon);

    /// Identifies the bounding box of pfemparticles and creates initial triangulation consisting of 2 triangles conecting bounding box nodes
    void buildInitialBBXMesh(InsertTriangleBasedOnCircumcircle &tInsert);
    /// Writes the mesh into the domain by creating new tr1_2d_pfem elements and prescribes zero-pressure boundary condition on alpha-shape nodes
    void writeMesh();

    /// Reads the triangulation and fills tha edgeList container with alpha-shape edges and set their bounds
    void computeAlphaComplex();

    /**
     * Fills the edgeList with unique alphaEdges. If an edge is already contained a pointer to it is returned and inserted edge is removed.
     */
    AlphaEdge2D *giveBackEdgeIfAlreadyContainedInList(AlphaEdge2D *alphaEdge);

    /// Iterates through the edgeList container and compares alpha-value with alphaEdge bounds. Alpha shape is stored in the alphaShapeEdgeList
    void giveAlphaShape();

    /// Initializes Timers and
    void initializeTimers();

    /// Looks for non-Delaunay triangles in octree and creates a polygon
    void findNonDelaunayTriangles(int insertedNode, InsertTriangleBasedOnCircumcircle &tInsert, std :: list< Edge2D > &polygon);
    /// Retriangulates the polygon
    void meshPolygon(int insertedNode, InsertTriangleBasedOnCircumcircle &tInsert, std :: list< Edge2D > &polygon);
    /// Prints the time report
    void giveTimeReport();
    /// Iterates through generalTringleList und removes non-valid ones or those containing bounding box nodes
    void cleanUpTriangleList();
};
} // end namespace oofem
#endif // delaunaytrinagulator_h
