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
#include "layeredcrosssection.h"

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
 * @author Jim Brouzoulis
 */
class EnrichmentItem : public FEMComponent
{
public:
    /// Constructor.
    EnrichmentItem(int n, XfemManager *xm, Domain *aDomain);

    virtual IRResultType initializeFrom(InputRecord *ir);
    int instanciateYourself(DataReader *dr);
    virtual const char *giveClassName() const { return "EnrichmentItem"; }

    /// Accessor. should there be support for several geom. objects describing one EI, probably yes. Ex. inclusion given as union of several geom.s?
    BasicGeometry *giveGeometry(int i);
    BasicGeometry *giveGeometry();

    // Spatial search queries
    
    /// Computes intersection points with Element. - based on the geometry of the enrichment
    void computeIntersectionPoints(AList< FloatArray > *intersectionPoints, Element *element);
    
    /// Computes number intersection points with Element.
    int computeNumberOfIntersectionPoints(Element *element);
    
    /// Checks whether a Geometry is inside or outside. - what if the geometry is not closed? like a set of segments
    bool isOutside(BasicGeometry *bg);


    /// Accessor.
    //EnrichmentFunction *giveEnrichmentFunction(){}; // but there may be several functions?
    EnrichmentFunction *giveEnrichmentFunction(int n);

    /// Gives number of dofs.
    //int giveNumberOfDofs() { return this->giveEnrichmentFunction()->giveNumberOfDofs(); }
    int giveNumberOfDofs() { return 1; } // should loop over all EF and ask them - what should be meant. active dofs? total or enriched?

    /// Sets DofId Array of an Enrichment Item.
    void setDofIdArray(IntArray &dofId) { this->dofsId = dofId; }
    /// Accessor.
    IntArray *getDofIdArray() { return & dofsId; }

    /// Finds out whether a DofManager is enriched.
    bool isDofManEnriched(int nodeNumber);
    /// Checks whether EnrichmentItem interacts element. - check if el. is enriched by this particular EI
   // bool interacts(Element *element){ return this->isElementEnriched(element); }; // old method
    bool isElementEnriched(Element *element); 


    
    virtual Material *giveMaterial() { return NULL; }
    
    /// Updates receiver geometry to the state reached at given time step.
    /// Geometry update; calls individual enrichment item updateGeometry method.
    virtual void updateGeometry(TimeStep *tStep) {}
        
protected:
    /// Link to associated Xfem manager.
    XfemManager *xmanager;
    /// Geometry associated with EnrichmentItem.
    int geometry;
    /// EnrichmentFunction associated with the EnrichmentItem. - should be a list of functions
    int enrichmentFunction;
    /// Additional dofIds from Enrichment. - depends on problem type and spatial dimension
    IntArray dofsId;

    // New -JB
    /// Geometry list.
    AList< BasicGeometry > *geometryList;
    /// Enrichment function list.
    AList< EnrichmentFunction > *enrichmentFunctionList;
    int numberOfEnrichmentFunctions;
    int numberOfGeometryItems;
};

/** Concrete representation of EnrichmentItem. */
class CrackTip : public EnrichmentItem // only for 2D. Only the tip element belong to this
{
public:
    CrackTip(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain) { }
};

/** Concrete representation of EnrichmentItem. */
class CrackInterior : public EnrichmentItem // rest of the crack el. that does not contain any tip 
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


class Delamination : public EnrichmentItem, LayeredCrossSection // rest of the crack el. that does not contain any tip 
{
    // maaybe layeredDelamination?
public:
    Delamination(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain), LayeredCrossSection(n, aDomain){}
    int numberOfDelaminations;
    IntArray delaminatedLayers; // ambigious
    FloatArray delaminationZCoords; // must they be ordered?
    void updateIntegrationRule();
    int giveDelaminationGroupAt();
    double giveDelaminationGroupMidZ();
    double giveDelaminationGroupThickness();

    /* for each delamination group = nDelam + 1
        int dGroup = giveDelaminationGroupAt(gp);
        double dMidZ = giveDelaminationGroupMidZ(dGroup);
        double dThickness = giveDelaminationGroupThickness(dGroup);
        b = dMidZ + 0.5*dThickness;
        a = dMidZ - 0.5*dThickness;
        remap the xi-coords -> 1-2(b-x)/(b-a)

    */

};



} // end namespace oofem




#endif  // enrichmentitem_h
