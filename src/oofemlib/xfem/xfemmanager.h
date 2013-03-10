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

#ifndef xfemmanager_h
#define xfemmanager_h

#include "alist.h"
#include "datareader.h"
#include "dofiditem.h"
#include "inputrecord.h"
#include "classtype.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "timestep.h"

namespace oofem {
class EngngModel;
class Domain;
class BasicGeometry;
class EnrichmentItem;
class EnrichmentFunction;
class IntArray;
class Element;
class DataStream;

/**
 * This class manages the xfem part as well as takes over some functions which would appear
 * in the Domain and Node class.
 *
 * @author Ruzena Chamrova
 * @author Jim Brouzoulis
 */
class XfemManager
{
protected:
    /// Associated Engineering Model.
    EngngModel *emodel;
    /// Index of the associated domain.
    int domainIndex;
    /// Enrichment item list.
    AList< EnrichmentItem > *enrichmentItemList;


    /// Index of next available dofId from pool.
    int startOfDofIdPool;
    int numberOfEnrichmentItems;

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
     * are referenced by a number from the domain - Don't like the name 'interacted' // JB
     */
    void giveActiveEIsFor(IntArray &answer, Element *elem);
    
    /// Checks whether an element is interacted.
    bool isElementEnriched(Element *elem);

    /// Checks whether a node is interacted. or 'isEnriched'
    bool isEnriched(int nodeNumber){ return isNodeEnriched(nodeNumber); };
    bool isNodeEnriched(int nodeNumber);

    /// Accessor.
    EnrichmentItem *giveEnrichmentItem(int n);
    
    /// Accessor.
    int giveNumberOfEnrichmentItems() { return enrichmentItemList->giveSize(); }
    
    /// Computes for each node position of its fictitious node. - What is this used for?
    void createEnrichedDofs();

    /// Computes the type of node enrichment, returns zero if the node is not enriched.
    XfemType computeNodeEnrichmentType(int nodeNumber); // ask node for EI and then type. but could be several?

    /// Initializes receiver according to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);

    /// Instantiates the Xfem components.
    int instanciateYourself(DataReader *dr);
   // const char *giveClassName() const { return "XfemManager"; }
     const char *giveClassName() const { return ""; }
    const char *giveInputRecordName() const { return "XfemManager"; }
    
    /// Wrapper for updating the integration rule.
    void updateIntegrationRule();

    Domain *giveDomain() { return emodel->giveDomain(domainIndex); }


    /// Clear the receiver
    void clear();


    /**
     * Stores the state of receiver to output stream.
     * @param stream Context stream.
     * @param mode Determines amount of info in stream.
     * @param obj Special parameter, used to pass optional parameters.
     * @return contextIOResultType.
     * @exception ContextIOERR If error encountered.
     */
   
    //contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the state of receiver from output stream.
     * @param stream Context file.
     * @param mode Determines amount of info in stream.
     * @param obj Special parameter for sending extra information.
     * @return contextIOResultType.
     * @exception ContextIOERR exception if error encountered.
     */
    //contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

};
} // end namespace oofem
#endif // xfemmanager_h
