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

#ifndef spatiallocalizer_h
#define spatiallocalizer_h

#include "femcmpnn.h"
#include "compiler.h"

#include "interface.h"
#include "logger.h"
#ifndef __MAKEDEPEND
 #include <set>
 #include <list>
#endif

namespace oofem {
class Domain;
class Element;
class TimeStep;
class GaussPoint;

/**
 * The spatial localizer element interface associated to spatial localizer.
 */
class SpatialLocalizerInterface : public Interface
{
public:
    SpatialLocalizerInterface() : Interface() { }

    /**
     * @name The element interface required by SpatialLocalizerInterface
     */
    //@{
    /// @return Reference to corresponding element.
    virtual Element *SpatialLocalizerI_giveElement() = 0;
    /**
     * Checks if element contains specified coordinate.
     * @param coords Global coordinate.
     * @return Nonzero if given element contains given point.
     */
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords) = 0;
    /**
     * Creates a bounding box of the nodes and checks if it includes the given coordinate.
     * @param coords Global coordinate.
     * @return Nonzero if given element bounding box contains given point.
     */
    int SpatialLocalizerI_BBoxContainsPoint(const FloatArray &coords);
    /**
     * Creates a bounding box of the nodes and checks if it includes the given coordinate.
     * @param bb0 Lower bounding box.
     * @param bb1 Upper bounding box.
     */
    virtual void SpatialLocalizerI_giveBBox(FloatArray &bb0, FloatArray &bb1);
    /**
     * Check the distance from the parametric center.
     * @param coords Global coordinate.
     * @return Distance of given point from element parametric center.
     */
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
    {
        OOFEM_ERROR2("SpatialLocalizerInterface :: SpatialLocalizerI_giveDistanceFromParametricCenter - Not implemented for %s", this->giveClassName());
        return 0.0;
    }

    /**
     * Gives the closest point on the element.
     * Default implementation uses the element interpolation.
     * @param[out] lcoords Local coordinates of closest point within the element.
     * @param[out] closest Global coordinates of closest point within the element.
     * @param gcoords Global coordinates.
     * @return Distance between answer and gcoords. Zero if gcoords is within the element, negative if point could not be found.
     */
    virtual double SpatialLocalizerI_giveClosestPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &gcoords);
    //@}
};


/**
 * The base class for all spatial localizers.
 * The basic task is to provide spatial information and localization for domain, to which receiver is associated.
 * Typical services include searching the closes node to give position, searching of an element containing given point, etc.
 * If special element algorithms required, these should be included using interface concept.
 */
class SpatialLocalizer : public FEMComponent
{
public:
    /// Typedefs to introduce the container type for element numbers, returned by some services.
    typedef std :: set< int >elementContainerType;
    /// Typedefs to introduce the container type for nodal numbers, returned by some services.
    typedef std :: list< int >nodeContainerType;

    /// Constructor
    SpatialLocalizer(int n, Domain *d) : FEMComponent(n, d) { }

    /**
     * Returns the element, containing given point and belonging to one of the region in region list.
     * @param coords Global problem coordinates of point of interest.
     * @param regionList Only elements within given regions are considered, if NULL all regions are considered.
     * @return The element belonging to associated domain, containing given point, NULL otherwise.
     */
    virtual Element *giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList = NULL) = 0;
    /**
     * Returns the element close to point
     * @param coords Global problem coordinates of point of interest.
     * @param regionList Only elements within given regions are considered, if NULL all regions are considered.
     * @return The element belonging to associated domain, close to given point, NULL otherwise.
     */
    virtual Element *giveElementCloseToPoint(const FloatArray &coords, const IntArray *regionList = NULL) = 0;
    /**
     * Returns the element closest to a given point.
     * @param[out] lcoords Local coordinates in element found.
     * @param[out] closest Global coordinates for found point.
     * @param coords Global problem coordinates of point of interest.
     * @param region Only elements within given region are considered, if 0 all regions are considered.
     * @return The element belonging to associated domain, close to given point, NULL otherwise.
     */
    virtual Element *giveElementClosestToPoint(FloatArray &lcoords, FloatArray &closest,
            const FloatArray &coords, int region = 0)
    {
        OOFEM_ERROR2("SpatialLocalizer :: giveElementClosestToPoint - Not implemented for %s", this->giveClassName());
        return NULL;
    }
    /**
     * Returns the integration point in associated domain, which is closest
     * to given point. Since IP holds the information about its element,
     * the IP reference is containing all the information.
     * @note{Only the gp belonging to the given region are taken into account.}
     * @param coords Global problem coordinates of point of interest
     * @param region If value > 0 then only closet point from given
     * region will be considered, if value < 0 all regions will be valid
     * @return The IP belonging to associated domain (only those provided by elements in default integration rule
     * are taken into account), NULL otherwise
     */
    virtual GaussPoint *giveClosestIP(const FloatArray &coords, int region) = 0;

    /**
     * Returns container (set) of all domain elements having integration point within given box.
     * @param elemSet Answer containing the list of elements meeting the criteria.
     * @param coords Center of box of interest.
     * @param radius Radius of bounding sphere.
     */
    virtual void giveAllElementsWithIpWithinBox(elementContainerType &elemSet, const FloatArray &coords,
                                                const double radius) = 0;
    /**
     * Returns container (set) of all domain elements having node within given box.
     * @param elemSet Answer containing the list of elements meeting the criteria.
     * @param coords Center of box of interest.
     * @param radius Radius of bounding sphere.
     */
    virtual void giveAllElementsWithNodesWithinBox(elementContainerType &elemSet, const FloatArray &coords,
                                                   const double radius);

    /**
     * Returns container (list) of all domain nodes within given box.
     * @param nodeList Answer containing the list of nodes meeting the criteria.
     * @param coords Center of box of interest.
     * @param radius Radius of bounding sphere.
     */
    virtual void giveAllNodesWithinBox(nodeContainerType &nodeList, const FloatArray &coords, const double radius) = 0;

    /**
     * Initialize receiver data structure if not done previously
     * If force is set to true, the initialization is enforced (useful if domain geometry has changed)
     * @return Nonzero if successful.
     */
    virtual int init(bool force = false) { return 1; }

    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    virtual const char *giveClassName() const { return "SpatialLocalizer"; }
    virtual classType giveClassID() const { return SpatialLocalizerClass; }
};
} // end namespace oofem
#endif // spatiallocalizer_h
