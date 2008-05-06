/* $Header: /home/cvs/bp/oofem/oofemlib/src/dummylocalizer.h,v 1.6 2003/04/06 14:08:23 bp Exp $ */
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


//   *************************************
//   *** CLASS DUMMY SPATIAL LOCALIZER ***
//   *************************************

#ifndef dummylocalizer_h
#define dummylocalizer_h

#include "spatiallocalizer.h"
#include "compiler.h"


class Domain;
class Element;
class TimeStep;

/**
 * The dummy implementation of spatial localizer based on traversing the whole domain.
 * The basic task is to provide spatial information and localization for domain, to which receiver is associated.
 * Typical services include searching the closes node to give position, serching of an element containing given point, etc.
 * If special element algorithms required, these should be included using interface concept.
 */
class DummySpatialLocalizer : public SpatialLocalizer
{
protected:

public:
    /// Constructor
    DummySpatialLocalizer(int n, Domain *d) : SpatialLocalizer(n, d) { }

    /**
     * Returns the element, containing given point.
     * @param coords global problem coordinates of point of interest
     * @param regionList only elements within given regions are considered, if NULL all regions are considered.
     * @return the element belonging to associated domain, containing given point, NULL otherwise
     */
    Element *giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList = NULL);
    /**
     * Returns the element close to point
     * @param coords global problem coordinates of point of interest
     * @param regionList only elements within given regions are considered, if NULL all regions are considered.
     * @return the element belonging to associated domain, close to given point, NULL otherwise
     */
    Element *giveElementCloseToPoint(const FloatArray &coords, const IntArray *regionList = NULL);
    /**
     * Returns the integration point in associated domain, which is closest
     * to given point. Since IP holds the information about its element,
     * the IP reference is containing all the information.
     * @param coords global problem coordinates of point of interest
     * @return the IP belonging to associated domain (only those provided by elements in default integration rule
     * are taken into acount), NULL otherwise
     */
    GaussPoint *giveClosestIP(const FloatArray &coords, int region);
    /**
     * Returns container (set) of all domain elements having integration point within given box.
     * @param elemSet answer containing the list of elements meeting the criteria
     * @param coords center of box of interest
     * @param radius radius of bounding sphere
     */
    void giveAllElementsWithIpWithinBox(elementContainerType &elemSet, const FloatArray &coords, const double radius);

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "DummySpatialLocalizer"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return DummySpatialLocalizerClass; }

protected:
};


#endif // dummylocalizer_h






