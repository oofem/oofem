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

#ifndef particletopologydescription_h
#define particletopologydescription_h

#include "topologydescription.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include <string>
#include <list>
#include <vector>
#include <memory>

///@todo These are general enough to be able to use outside this class.
///@name Input fields for Circle
//@{
#define _IFT_Circle_center "center"
#define _IFT_Circle_radius "radius"
#define _IFT_Circle_start "start"
#define _IFT_Circle_end "end"
//@}

///@name Input fields for Line
//@{
#define _IFT_Line_start "start"
#define _IFT_Line_end "end"
//@}

///@name Input fields for Point
//@{
#define _IFT_Point_coords "coords"
//@}

///@name Input fields for Meshing
//@{
#define _IFT_Meshing_elementType "elementtype"
#define _IFT_Meshing_set "set"
//@}

///@name Input fields for ParticleTopology
//@{
#define _IFT_ParticleTopologyDescription_Name "particletopology"
#define _IFT_ParticleTopologyDescription_nsd "nsd"
#define _IFT_ParticleTopologyDescription_baseResolution "baseresolution"
#define _IFT_ParticleTopologyDescription_tubeWidth "tubewidth"
#define _IFT_ParticleTopologyDescription_neighbors "neighbors"
#define _IFT_ParticleTopologyDescription_boundingBoxA "bboxa"
#define _IFT_ParticleTopologyDescription_boundingBoxB "bboxb"
#define _IFT_ParticleTopologyDescription_numberOfSegments "nsegments"
#define _IFT_ParticleTopologyDescription_regionOutside "regionoutside"
#define _IFT_ParticleTopologyDescription_regionInside "regioninside"
#define _IFT_ParticleTopologyDescription_identification "id"
//@}


namespace oofem
{
template< class Point >class ParticleGrid;
class Domain;
struct Triangle_PSLG;

/// Default point type for describing topology.
struct ParticlePoint {
    ParticlePoint() {};
    ParticlePoint(const FloatArray &foot, int id, const FloatArray &normal, double distance2) :
        foot(foot), id(id), normal(normal), distance2(distance2), removal(false) {};
    FloatArray foot; ///< Closest coordinate on surface.
    // Auxiliary information
    FloatArray total_displacement; ///< Total displacement since last resampling.
    int id; ///< Id of surface.
    FloatArray normal; ///< Surface normal at foot point.
    FloatArray corner; ///< Corner particle (for open surfaces).
    double distance2; ///< Squared distance stored for efficiency.
    int node; ///< Node number (for meshing)
    int c_node; ///< Corner node number (for open surfaces).
    bool removal; /// Convenience variable for clearing out nodes which have conflicting information.
};

/**
 * A grid based particle method for describing topology.
 * Based on papers:
 * Shingyu Leung and Hongkai Zhao, A Grid Based Particle Method for Evolution of Open Curves and Surfaces. Journal of Computational Physics, Volume 228, Issue 20, November 1 2009, Pages 7706-7728.
 * Shingyu Leung and Hongkai Zhao, A Grid Based Particle Method for Moving Interface Problems. UCLA-CAM 08-08. Journal of Computational Physics, Volume 228, Issue 8, May 1 2009, Pages 2993-3024.
 * with some modifications and missing implementations.
 *
 * @note Experimental! Both inconvenient to use, and not very robust.
 *
 * @todo Absolutely vital is the use of adaptivity, which is not yet implemented.
 * @todo Merging of regions is still lacking.
 * @todo Vanishing pores occasionally experience problems.
 * @todo Partial reconstruction should be possible.
 *
 * @author Mikael Öhman
 */
class ParticleTopologyDescription : public TopologyDescription
{
protected:
    /// Determines if velocity or displacements dofs should be used to update geometry.
    bool useDisplacements;
    /// Conditional for printing VTK output.
    bool writeVTK;
    /// Denotes if the active grid is newly resampled.
    bool resampled;
    /// Maximum squared displacement of any particle.
    double maxdisp2;
    /// Width of the tube around the interfaces
    double tubeWidth;
    /// Number of points to use for resampling
    int m;
    /// The grid of points, the actual topological information
    std :: unique_ptr< ParticleGrid< ParticlePoint > > grid;

    /// Corner nodes
    std :: list< ParticlePoint >corners;

    /**
     * Mapping of regions from delimited by each id.
     * This is only needed for remeshing.
     * Regions should be strictly different from any id in the domain, as the id's themselves map directly to regions.
     * First value is the region on the normals direction, the second value is the other side.
     */
    IntArray regionOutside, regionInside;
    FloatMatrix mergeID, controlID;

    /**
     * Mapping from region to FE components.
     */
    std :: vector< std :: string >regionElementType;
    IntArray regionSet;

    /**
     * Resamples the grid. Entire grid is replaced.
     * Finds new foot points for new grid points along the curve.
     */
    void resample();

    /**
     * Deactivates points with inconsistent information.
     * This is typically due to merging surfaces or vanishing pores.
     * @return TS_NeedsRemeshing if points have been deactivated.
     */
    TopologyState checkOverlap();

    /**
     * Clears all points marked for removal.
     * Works recursively on refined grids.
     * @param g Grid to clear.
     */
    void removePoints(ParticleGrid< ParticlePoint > &g) const;

    /**
     * Finds the displacement for the underlying FE-mesh.
     * @param answer The requested displacement.
     * @param id Only elements of specified region is taken into account (unless id == 0).
     * @param footpoint The spatial coordinate where to find the displacement.
     * @param tStep The time step for which to find the displacement.
     */
    bool findDisplacement(FloatArray &answer, int id, const FloatArray &footpoint, TimeStep *tStep) const;

    /**
     * Collects neighboring points according to some specification.
     * Finds a maximum of around m neighbors, could be fewer.
     * A optional distance can be added to take into account the displaced foot points.
     * @param answer The neighboring particles.
     * @param p Point to compute neighbors around.
     * @param dist Extra distance to take neighbors from.
     */
    void collectNeighbors(std :: list< ParticlePoint * > &answer, const ParticlePoint *p, double dist = 0) const;

    /**
     * Helper for common task of fetching a bounding box around a point.
     */
    static void getBoundingBox(FloatArray &x0, FloatArray &x1, const FloatArray &c, double width);

    /**
     * Shortest distance from least square fit based on 2nd order polynomial.
     */
    void calculateShortestDistance(const ParticlePoint *p, std :: list< ParticlePoint * > &points, ParticleGrid< ParticlePoint > &grid) const;

    /**
     * Helper for calculateShortestDistance
     */
    double shortestDistanceFromCurve(const FloatArray &a, double txi_min, double txi_max,
                                     const FloatArray &n0, const FloatArray &y0, const FloatArray &p, FloatArray &foot, FloatArray &normal) const;

    /**
     * Used for initialization, calculating the distance from primitives.
     * Adds line segment from p0 to p1.
     * @param id ID for segment.
     * @param p0 First corner of the edge.
     * @param p1 Second corner of the edge.
     * @param grid Grid to add points to.
     */
    void addLineSegment(int id, const FloatArray &p0, const FloatArray &p1, ParticleGrid< ParticlePoint > &grid) const;
    /**
     * Used for initialization, calculating the distance from primitives.
     * Adds a circle segment from p0 to p1, with center at c.
     * @param id ID for segment.
     * @param c Center coordinate.
     * @param r Radius.
     * @param v0 Lower angle [-pi,pi].
     * @param v1 Upper angle [-pi,pi].
     * @param grid Grid to add points to.
     */
    void addCircleSegment(int id, const FloatArray &c, double r, double v0, double v1, ParticleGrid< ParticlePoint > &grid) const;
    /**
     * Adds a corner node.
     * @param id ID for corner.
     * @param c Coordinate of corner.
     * @param grid Grid to add corner to.
     */
    void addCorner(int id, const FloatArray &c, ParticleGrid< ParticlePoint > &grid);

    /**
     * Generates the PSLG for meshing with Triangle.
     * @param PSLG Generated structure ready for meshing.
     */
    void generatePSLG(Triangle_PSLG &PSLG);

public:
    ParticleTopologyDescription(Domain *d);
    virtual ~ParticleTopologyDescription();
    virtual bool instanciateYourself(DataReader *dr);

    virtual TopologyState updateYourself(TimeStep *tStep);

    /**
     * Generates a mesh from the topology.
     * @param nodes Nodes created.
     * @param elements Bulk elements created.
     * @param segments Edge/surface elements describing the topology.
     * @param n_markers Node markers.
     * @param e_markers Element markers.
     * @param s_markers Segment markers.
     * @param e_egt Element geometry types.
     * @param s_egt Segment geometry types.
     */
    virtual void generateMesh(std :: vector< FloatArray > &nodes, std :: vector< IntArray > &elements, std :: vector< IntArray > &segments,
                              std :: vector< IntArray > &n_markers, IntArray &e_markers, IntArray &s_markers, IntArray &e_egt, IntArray &s_egt);

    virtual void replaceFEMesh();

    virtual void doOutput(TimeStep *tStep);
    virtual void writeDataToFile(const char *name) const;
    virtual void writeVTKFile(const char *name) const;

    virtual const char *giveClassName() const { return "ParticleTopologyDescription"; }
};
} // end namespace oofem
#endif // particletopologydescription_h
