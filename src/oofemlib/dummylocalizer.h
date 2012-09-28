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

#ifndef dummylocalizer_h
#define dummylocalizer_h

#include "spatiallocalizer.h"
#include "alist.h"

namespace oofem {
class IntArray;

/**
 * The dummy implementation of spatial localizer based on traversing the whole domain.
 * The basic task is to provide spatial information and localization for domain, to which receiver is associated.
 * Typical services include searching the closes node to give position, searching of an element containing given point, etc.
 * If special element algorithms required, these should be included using interface concept.
 */
class DummySpatialLocalizer : public SpatialLocalizer
{
protected:
    AList<IntArray> region_elements;
    bool initialized;

public:
    /// Constructor
    DummySpatialLocalizer(int n, Domain *d) : SpatialLocalizer(n, d), initialized(false) { }
    /// Destructor
    virtual ~DummySpatialLocalizer() { this->region_elements.clear(true); }

    virtual int init(bool force = false);

    virtual Element *giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList = NULL);
    virtual Element *giveElementCloseToPoint(const FloatArray &coords, const IntArray *regionList = NULL);
    virtual Element *giveElementClosestToPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &coords, int region = 0);
    virtual GaussPoint *giveClosestIP(const FloatArray &coords, int region);
    virtual void giveAllElementsWithIpWithinBox(elementContainerType &elemSet, const FloatArray &coords, const double radius);
    virtual void giveAllNodesWithinBox(nodeContainerType &nodeList, const FloatArray &coords, const double radius);

    virtual const char *giveClassName() const { return "DummySpatialLocalizer"; }
    virtual classType giveClassID() const { return DummySpatialLocalizerClass; }
};
} // end namespace oofem
#endif // dummylocalizer_h
