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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef aabb_h
#define aabb_h

namespace oofem {

class Vector {
  public:
    double x, y, z;
    Vector();
    Vector(double x, double y, double z);
    double operator[](int index) const;
};

/**
 * @brief Axis-aligned bounding box.
 *
 * The AABB class represents an axis-aligned bounding box defined in global
 * coordinates. It is primarily used in contact detection and spatial search
 * algorithms to efficiently test for spatial proximity and overlap between
 * geometric entities.
 *
 * An AABB is characterized by its minimum and maximum corner coordinates
 * aligned with the global coordinate axes. This alignment enables fast
 * intersection tests and simple updates as objects move or deform.
 *
 * The class provides basic operations for initializing, updating, and
 * querying bounding boxes, and serves as a lightweight geometric utility
 * within the contact and search infrastructure.
 */
  
class AABB {
  public:
    Vector min;
    Vector max;
  /**
   * @brief Default constructor.
   *
   * Creates an empty or uninitialized axis-aligned bounding box.
   * The bounding box is typically initialized later by expanding
   * it to include points or other bounding boxes.
   */
    AABB() = default;
  /**
   * @brief Constructs an axis-aligned bounding box from corner coordinates.
   *
   * Initializes the bounding box using the given minimum and maximum
   * corner coordinates in global space.
   *
   * @param minCorner Vector of minimum coordinates along each axis.
   * @param maxCorner Vector of maximum coordinates along each axis.
   */
    AABB(const Vector& min, const Vector& max);
  /**
   * @brief Checks whether a point is contained within the bounding box.
   *
   * Tests if the given point lies inside or on the boundary of this
   * axis-aligned bounding box. The check is performed independently
   * along each coordinate axis.
   *
   * @param v Vector in global coordinates.
   * @return True if the point is inside or on the boundary of the box,
   *         false otherwise.
   */
    bool contains(const Vector& v);
    bool contains(double x, double y, double z);
  /**
   * @brief Expands the bounding box to include a given point.
   *
   * Updates the minimum and maximum corner coordinates of the bounding box
   * such that the point @p v is contained within the box. For each coordinate
   * direction, the bounds are enlarged only if the point lies outside the
   * current extent.
   *
   * @param v Point in global coordinates to be merged into the bounding box.
   */
  void merge(const Vector& v);
  /**
   * @brief Expands the bounding box to include a given point.
   *
   * Updates the minimum and maximum corner coordinates of the bounding box
   * such that the point @p v is contained within the box. For each coordinate
   * direction, the bounds are enlarged only if the point lies outside the
   * current extent.
   *
   * @param x X-coordinate of the point.
   * @param y Y-coordinate of the point.
   * @param z Z-coordinate of the point.
   */
  void merge(double x, double y, double z);
};

} // end namespace oofem
#endif // aabb_h
