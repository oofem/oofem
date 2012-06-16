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

#ifndef enrichmentitem_h
#define enrichmentitem_h

#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "enrichmentfunction.h"

namespace oofem {
class BasicGeometry;

/**
 * Abstract class representing entity, which is included in the FE model using one (or more)
 * global functions. Such entity may represent crack, material interface, etc.
 * As the geometry of such entity may be represented in a number of ways, the hierarchy of classes
 * derived from base Geometry class is used to achieve flexibility of geometry representation.
 *
 * Each EnrichmentItem keeps its DOF labels (assigned/allocated by XFemManager, its geometry representation, and
 * keeps the list of its EnrichmentFunctions.
 * @author chamrova
 */
class EnrichmentItem : public FEMComponent
{
public:
    /// Constructor.
    EnrichmentItem(int n, XfemManager *xm, Domain *aDomain);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "EnrichmentItem"; }

    /// Accessor.
    BasicGeometry *giveGeometry();
    /// Checks whether EnrichmentItem interacts element.
    bool interacts(Element *element);
    /// Computes intersection points with Element.
    void computeIntersectionPoints(AList< FloatArray > *intersectionPoints, Element *element);
    /// Computes number intersection points with Element.
    int computeNumberOfIntersectionPoints(Element *element);
    /// Accessor.
    EnrichmentFunction *giveEnrichmentFunction();
    /// Gives number of dofs.
    int giveNumberOfDofs() { return this->giveEnrichmentFunction()->giveNumberOfDofs(); }
    /// Sets DofId Array of an Enrichment Item.
    void setDofIdArray(IntArray &dofId) { this->dofsId = dofId; }
    /// Accessor.
    IntArray *getDofIdArray() { return & dofsId; }
    /// Finds out whether a DofManager is enriched.
    bool isDofManEnriched(int nodeNumber);
    /// Checks whether a Geometry is inside or outside.
    bool isOutside(BasicGeometry *bg);
    virtual Material *giveMaterial() { return NULL; }
    /// Updates receiver geometry to the state reached at given time step.
    virtual void updateGeometry(TimeStep *tStep) {}
protected:
    /// Link to associated Xfem manager.
    XfemManager *xmanager;
    /// Geometry associated with EnrichmentItem.
    int geometry;
    /// EnrichmentFunction associated with the EnrichmentItem.
    int enrichmentFunction;
    /// Additional dofIds from Enrichment.
    IntArray dofsId;
};

/** Concrete representation of EnrichmentItem. */
class CrackTip : public EnrichmentItem
{
public:
    CrackTip(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain) { }
};

/** Concrete representation of EnrichmentItem. */
class CrackInterior : public EnrichmentItem
{
public:
    CrackInterior(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain) { }
};

/** Concrete representation of EnrichmentItem. */
class Inclusion : public EnrichmentItem
{
protected:
    Material *mat;

public:
    Inclusion(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain) { }
    virtual const char *giveClassName() const { return "Inclusion"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Material *giveMaterial() { return mat; }
};
} // end namespace oofem
#endif  // enrichmentitem_h
