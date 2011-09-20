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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef xfemmanager_h
#define xfemmanager_h

#include "alist.h"

namespace oofem {
class EngngModel;
class Domain;
class BasicGeometry;
class EnrichmentItem;
class EnrichmentFunction;


/**
 * This class manages the xfem part as well as takes over some functions which would appear
 * in the Domain and Node class.
 * @author chamrova
 */
class XfemManager
{
protected:
    /// Associated EngineeringModel.
    EngngModel *emodel;
    /// Index of the associated domain.
    int domainIndex;
    /// Enrichment item list.
    AList< EnrichmentItem > *enrichmentItemList;
    /// Geometry list.
    AList< BasicGeometry > *geometryList;
    /// Enrichment function list.
    AList< EnrichmentFunction > *enrichmentFunctionList;
    /// Map giving for a node a position of its fictitious node.
    AList< IntArray > *fictPosition;
    /// Index of next available dofId from pool.
    int dofIdPos;
    int numberOfEnrichmentItems;
    int numberOfEnrichmentFunctions;
    int numberOfGeometryItems;

public:
    enum XfemType {
        SPLIT = 1, TIP = 4, STANDARD = 0
    };
    /// Constructor.
    XfemManager(EngngModel *emodel, int index);
    /// Destructor.
    ~XfemManager();
    /**
     * Gets interacted enrichment items for a particular element, the enrichment items
     * are referenced by a number from the domain
     */
    void getInteractedEI(IntArray &answer, Element *elem);
    /// checks whether an element is interacted
    bool isInteracted(Element *elem);
    /// checks whether a node is interacted
    bool isEnriched(int nodeNumber);
    /// Accessor.
    EnrichmentItem *giveEnrichmentItem(int n);
    /// Accessor.
    BasicGeometry *giveGeometry(int n);
    /// Accessor.
    EnrichmentFunction *giveEnrichmentFunction(int n);
    int giveNumberOfEnrichmentItems() { return enrichmentItemList->giveSize(); }
    /// Computes for each node position of its fictitious node.
    int computeFictPosition();
    /// Computes the type of node enrichment, returns zero if the node is not enriched.
    XfemType computeNodeEnrichmentType(int nodeNumber);
    /// Initializes receiver according to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);
    /// Instantiates the Xfem components.
    int instanciateYourself(DataReader *dr);
    const char *giveClassName() const { return "XfemManager"; }
    const char *giveInputRecordName() const { return "XfemManager"; }
    /// Wrapper for updating the integration rule.
    void updateIntegrationRule();
    /// Gives Domain.
    Domain *giveDomain();
    /// Accessor.
    IntArray *giveFictPosition(int nodeNumber) { return fictPosition->at(nodeNumber); }
    /// Stores receiver state to output stream.
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL) { return CIO_OK; }
    /// Restores the receiver state previously written in stream.
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL) { return CIO_OK; }
    /// Geometry update; calls individual enrichment item updateGeometry method.
    void updateGeometry(TimeStep *tStep);

protected:
    // Changes dofIdPos to next index.
    DofIDItem allocateNewDofID();
};
} // end namespace oofem
#endif // xfemmanager_h
