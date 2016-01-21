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

#ifndef trianglemesherinterface_h
#define trianglemesherinterface_h

#include "floatarray.h"
#include "intarray.h"
#include "oofemcfg.h"

#include <vector>

namespace oofem {
/**
 * Plane straight line graph used as input for meshing with triangle.
 */
struct Triangle_PSLG
{
    FloatArray nx; ///< Nodes x coordinates.
    FloatArray ny; ///< Nodes y coordinates.
    IntArray segment_a; ///< First segment connection.
    IntArray segment_b; ///< Second segment connection.
    IntArray segment_marker; ///< Segment markers
};

/**
 * Interface to Triangle (Delaunay mesher).
 *
 * @todo{Possible reform slightly and merge this with the other meshing interfaces.}
 * @see http://www.cs.cmu.edu/~quake/triangle.html
 * @author Mikael Öhman
 */
class OOFEM_EXPORT TriangleMesherInterface
{
protected:
    double minAngle;
    double maxArea;
    bool useRegions;
    bool quadratic;

public:
    /**
     * Constructor.
     * @param minAngle Minimum angle. Should be less than 33 (20 is default).
     * @param maxArea Maximum allowed area of a triangle.
     * @param quadratic True if generated mesh should be quadratic (6 nodes).
     */
    TriangleMesherInterface(double minAngle, double maxArea, bool quadratic) :
        minAngle(minAngle), maxArea(maxArea), quadratic(quadratic) { }
    /// Destructor.
    ~TriangleMesherInterface() { }

    /**
     * Simplifies a PSLG while respecting topology, running in linear time.
     * The algorithm removes unnecessary nodes while preserving topology.
     * This means any node with other than 2 edges are kept.
     * The error measure is is the distance from any point to the simplified curve.
     * A straight line will be simplified to its end points, even for zero limit value.
     * Direction of edges is preserved.
     * @param coarse Resulting PSLG.
     * @param pslg Input PSLG.
     * @param limit Maximum difference allowed from new to old line.
     * @param minlen Minimum length of any segment (shorter segments will be forcibly removed). Will not respect the error limit.
     */
    static void simplifyPSLG(Triangle_PSLG &coarse, const Triangle_PSLG &pslg, double limit, double minlen = 0.0);

    /**
     * Constructs a mesh from a PSLG.
     * The PSLG needs to be oriented correctly, as every segment determines the bulk region number at both sides,
     * which then spreads to surrounding elements. The length of those arrays must be at least as big as the largest segment region number.
     * Inconsistent regions are automatically detected.
     * @param[in] pslg Input PSLG.
     * @param[in] outside Segment region to bulk region mapping.
     * @param[in] inside Segment region to bulk region mapping.
     * @param[out] nodes Output nodes.
     * @param[out] n_markers Node markers.
     * @param[out] triangles Output triangles (3 or 6 nodes depending on option).
     * @param[out] t_markers Triangle markers
     * @param[out] segments Output segments (2 or 3 nodes depending on option).
     * @param[out] s_markers Segment markers.
     * @return True if mesh generation was successful.
     */
    bool meshPSLG(const Triangle_PSLG &pslg,
                  const IntArray &outside, const IntArray &inside,
                  std :: vector< FloatArray > &nodes, std :: vector< IntArray > &n_markers,
                  std :: vector< IntArray > &triangles, IntArray &t_markers,
                  std :: vector< IntArray > &segments, IntArray &s_markers) const;

protected:
    /**
     * Adds all neighboring regions to every node region.
     * Necessary since triangle can only store a single node number.
     */
    static void fixNodeMarkers(const std :: vector< FloatArray > &nodes, std :: vector< IntArray > &n_markers,
                               const std :: vector< IntArray > &triangles, const IntArray &t_markers,
                               const std :: vector< IntArray > &segments, const IntArray &s_markers);
};
} // end namespace oofem
#endif // trianglemesherinterface_h
