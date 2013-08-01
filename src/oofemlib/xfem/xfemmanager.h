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
#include "inputrecord.h"
#include "contextioresulttype.h"
#include "contextmode.h"

///@name Input fields for XfemManager
//@{
#define _IFT_XfemManager_Name "xfemmanager"
#define _IFT_XfemManager_numberOfGeometryItems "numberofgeometryitems"  // -> numberOfEnrichmentDomains
#define _IFT_XfemManager_numberOfEnrichmentItems "numberofenrichmentitems"
#define _IFT_XfemManager_numberOfEnrichmentFunctions "numberofenrichmentfunctions"
//@}

namespace oofem {
class Domain;
class EnrichmentItem;
//class EnrichmentFunction;
class IntArray;
class Element;
class DataStream;

/**
 * This class manages the xfem part
 *
 * @author Ruzena Chamrova
 * @author Jim Brouzoulis
 */
class XfemManager
{
protected:
    Domain *domain;
    /// Enrichment item list.
    AList< EnrichmentItem > *enrichmentItemList;

    /// Index of next available dofId from pool.
    int numberOfEnrichmentItems;

public:
    enum XfemType { // not in use right now
        SPLIT = 1, TIP = 4, STANDARD = 0
    };
    /// Constructor.
    XfemManager(Domain *domain);
    /// Destructor.
    ~XfemManager();


    // Returns the active enrichment items for a particular element, the enrichment items
    // are referenced by a number from the domain
    void giveActiveEIsFor(IntArray &answer, const Element *elem);

    bool isElementEnriched(const Element *elem);

    bool isAllElNodesEnriched(const Element *elem);

    /// Accessor.
    EnrichmentItem *giveEnrichmentItem(int n);
    int giveNumberOfEnrichmentItems() { return enrichmentItemList->giveSize(); }

    void createEnrichedDofs();

    /// Computes the type of node enrichment, returns zero if the node is not enriched.
    // Old method: should instead return an array if there are several active /JB
    XfemType computeNodeEnrichmentType(int nodeNumber); 

    /// Initializes receiver according to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);

    int instanciateYourself(DataReader *dr);
    const char *giveClassName() const { return "XfemManager"; }
    const char *giveInputRecordName() const { return _IFT_XfemManager_Name; }
    
    /// Wrapper for updating the integration rule.
    void updateIntegrationRule();

    //Domain *giveDomain() { return emodel->giveDomain(domainIndex); }
    Domain *giveDomain() { return this->domain; }

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
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the state of receiver from output stream.
     * @param stream Context file.
     * @param mode Determines amount of info in stream.
     * @param obj Special parameter for sending extra information.
     * @return contextIOResultType.
     * @exception ContextIOERR exception if error encountered.
     */
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);


    // Update level sets
    void updateYourself();
};
} // end namespace oofem
#endif // xfemmanager_h
