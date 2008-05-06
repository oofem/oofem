/* $Header: /home/cvs/bp/oofem/oofemlib/src/spatiallocalizer.h,v 1.9.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


//   *******************************
//   *** CLASS SPATIAL LOCALIZER ***
//   *******************************

#ifndef spatiallocalizer_h
#define spatiallocalizer_h

#include "femcmpnn.h"
#include "compiler.h"

#include "interface.h"
#include "logger.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <set>
#include <list>
#endif

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

    /// Returns reference to corresponding element
    virtual Element *SpatialLocalizerI_giveElement() = 0;
    /// Returns nonzero if given element contains given point
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords) = 0;
    /// Returns nonzero if given element bounding box contains given point
    int SpatialLocalizerI_BBoxContainsPoint(const FloatArray &coords);
    /// Returns distance of given point from element parametric center
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
    { OOFEM_ERROR("SpatialLocalizerI_giveDistanceFromParametricCenter: not implemented");
      return 0.0; }
};


/**
 * The base class for all spatial localizers.
 * The basic task is to provide spatial information and localization for domain, to which receiver is associated.
 * Typical services include searching the closes node to give position, serching of an element containing given point, etc.
 * If special element algorithms required, these should be included using interface concept.
 */
class SpatialLocalizer : public FEMComponent
{
protected:

public:

    /**
     * Typedefs to introduce the container type for element numbers, returned by some services
     */
    typedef std :: set< int > elementContainerType;
    /**
     * Typedefs to introduce the container type for nodal numbers, returned by some services
     */
    typedef std :: list< int > nodeContainerType;


    /// Constructor
    SpatialLocalizer(int n, Domain *d) : FEMComponent(n, d) { }

    /**
     * Returns the element, containing given point and belonging to one of the region in region list.
     * @param coords global problem coordinates of point of interest
     * @param regionList only elements within given regions are considered, if NULL all regions are considered.
     * @return the element belonging to associated domain, containing given point, NULL otherwise
     */
    virtual Element *giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList = NULL) = 0;
    /**
     * Returns the element close to point
     * @param coords global problem coordinates of point of interest
     * @param regionList only elements within given regions are considered, if NULL all regions are considered.
     * @return the element belonging to associated domain, close to given point, NULL otherwise
     */
    virtual Element *giveElementCloseToPoint(const FloatArray &coords, const IntArray *regionList = NULL) = 0;
    /**
     * Returns the integration point in associated domain, which is closest
     * to given point. Since IP holds the information about its element,
     * the IP reference is containing all the information.
     * NOTE: Only the gp belonging to the given region are taken into account.
     * @param coords global problem coordinates of point of interest
     * @param region - if  value > 0 then only closet point from given
     * region will be considered, if value < 0 all regions will be valid
     * @return the IP belonging to associated domain (only those provided by elements in default integration rule
     * are taken into acount), NULL otherwise
     */
    virtual GaussPoint *giveClosestIP(const FloatArray &coords, int region) = 0;

    /**
     * Returns container (set) of all domain elements having integration point within given box.
     * @param elemSet answer containing the list of elements meeting the criteria
     * @param coords center of box of interest
     * @param radius radius of bounding sphere
     */
    virtual void giveAllElementsWithIpWithinBox(elementContainerType &elemSet, const FloatArray &coords,
                                                const double radius) = 0;


    /** Initializes receiver acording to object description stored in input record.
     * This function is called immediately after creating object using
     * constructor. Input record can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.
     * @see readInteger, readDouble and similar functions */
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "SpatialLocalizer"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return SpatialLocalizerClass; }

protected:
};

#endif // spatiallocalizer_h






