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

#ifndef xfemmanager_h
#define xfemmanager_h

#include "oofemcfg.h"
#include "datareader.h"
#include "inputrecord.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "enrichmentitem.h"
#include "enumitem.h"
#include "internalstatevaluetype.h"

#include <unordered_map>
#include <list>
#include <vector>
#include <memory>

///@name Input fields for XfemManager
//@{
#define _IFT_XfemManager_Name "xfemmanager"
#define _IFT_XfemManager_numberOfEnrichmentItems "numberofenrichmentitems"
#define _IFT_XfemManager_numberOfNucleationCriteria "numberofnucleationcriteria"
#define _IFT_XfemManager_numberOfGpPerTri "numberofgppertri"

/// How many times a subtriangle should be refined
#define _IFT_XfemManager_numberOfTriRefs "numberoftrirefs"

#define _IFT_XfemManager_enrDofScaleFac "enrdofscalefac"

#define _IFT_XfemManager_debugVTK "debugvtk"
#define _IFT_XfemManager_VTKExport "vtkexport"
#define _IFT_XfemManager_VTKExportFields "exportfields"
//@}

//#define ENABLE_XFEM_CPP11

namespace oofem {
class Domain;
class EnrichmentItem;
class IntArray;
class Element;
class DataStream;
class DynamicInputRecord;
class NucleationCriterion;
//class InternalStateValueType;

//
// The following types determine the state types associated with xfem
//
#define XFEMStateType_DEF \
    ENUM_ITEM_WITH_VALUE(XFEMST_Undefined, 0) \
    ENUM_ITEM_WITH_VALUE(XFEMST_Enrichment, 1) \
    ENUM_ITEM_WITH_VALUE(XFEMST_LevelSetPhi, 2) \
    ENUM_ITEM_WITH_VALUE(XFEMST_LevelSetGamma, 3) \
    ENUM_ITEM_WITH_VALUE(XFEMST_NodeEnrMarker, 4) \
    ENUM_ITEM_WITH_VALUE(XFEMST_NumIntersecPoints, 5)

enum XFEMStateType {
    XFEMStateType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h

const char *__XFEMStateTypeToString(XFEMStateType _value);

/**
 * This class manages the xfem part
 *
 * @author Ruzena Chamrova
 * @author Jim Brouzoulis
 * @author Erik Svenning
 */
class OOFEM_EXPORT XfemManager
{
protected:
    Domain *domain;
    /// Enrichment item list.
    std :: vector< std :: unique_ptr< EnrichmentItem > >enrichmentItemList;

    int numberOfEnrichmentItems;

    int numberOfNucleationCriteria;

    /**
     * The number of Gauss points to be used in each sub-triangle when
     * subdividing cut elements.
     */
    int mNumGpPerTri;

    /**
     * The number of times a subtriangle should be refined.
     */
    int mNumTriRef;

    /* Scale factor for the enrichment dofs. Implies that the corresponding
       dofs must be scaled with 1/factor in the input file
     */
    double mEnrDofScaleFac;

    bool doVTKExport;

    /// If extra debug vtk files should be written.
    bool mDebugVTK;

    /**
     * Let the XfemManager keep track of enrichment items enriching each
     * node and each element, to allow more efficient computations.
     */
    std :: vector< std :: vector< int > >mNodeEnrichmentItemIndices;
    std :: unordered_map< int, std :: vector< int > >mElementEnrichmentItemIndices;

    /**
     * Keep track of enrichment items that may assign a different
     * material to some Gauss points.
     */
    std :: vector< int >mMaterialModifyingEnrItemIndices;

    /**
     * Nucleation of new enrichment items. (For example, nucleation of new cracks.)
     */
    std::vector< std :: unique_ptr< NucleationCriterion > > mNucleationCriteria;

    // IDs of all potential enriched dofs
    IntArray mXFEMPotentialDofIDs;

public:

    /**
     * List with the fields that should be exported to VTK
     */
    IntArray vtkExportFields; // make private later
    InternalStateValueType giveXFEMStateValueType(XFEMStateType type);

    /// Constructor.
    XfemManager(Domain *domain);
    /// Destructor.
    virtual ~XfemManager();
    XfemManager(const XfemManager &) = delete;

    XfemManager &operator=(const XfemManager &) = delete;

    int giveNumGpPerTri() const { return mNumGpPerTri; } /// Number of Gauss points per sub-triangle in cut elements.
    int giveNumTriRefs() const { return mNumTriRef;}
    double giveEnrDofScaleFactor() const {return mEnrDofScaleFac;}

    bool isElementEnriched(const Element *elem);

    inline EnrichmentItem *giveEnrichmentItem(int n) { return enrichmentItemList [ n - 1 ].get(); }
    int giveNumberOfEnrichmentItems() const { return ( int ) enrichmentItemList.size(); }

    inline NucleationCriterion *giveNucleationCriterion(int n) { return mNucleationCriteria [ n - 1 ].get(); }
    int giveNumberOfNucleationCriteria() const { return ( int ) mNucleationCriteria.size(); }

    void createEnrichedDofs();
    const IntArray &giveEnrichedDofIDs() const {return mXFEMPotentialDofIDs;}
    IntArray giveEnrichedDofIDs(const DofManager &iDMan) const;

    /// Initializes receiver according to object description stored in input record.
    virtual void initializeFrom(InputRecord &ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual int instanciateYourself(DataReader &dr);
    virtual const char *giveClassName() const { return "XfemManager"; }
    virtual const char *giveInputRecordName() const { return _IFT_XfemManager_Name; }

    virtual void postInitialize();

    Domain *giveDomain() { return this->domain; }
    void setDomain(Domain *ipDomain);

    /**
     * Stores the state of receiver to output stream.
     * @param stream Context stream.
     * @param mode Determines amount of info in stream.
     * @exception ContextIOERR If error encountered.
     */
    void saveContext(DataStream &stream, ContextMode mode);
    /**
     * Restores the state of receiver from output stream.
     * @param stream Context file.
     * @param mode Determines amount of info in stream.
     * @exception ContextIOERR exception if error encountered.
     */
    void restoreContext(DataStream &stream, ContextMode mode);


    /**
     * Update enrichment items (level sets).
     */
    virtual void updateYourself(TimeStep *tStep);

    virtual void propagateFronts(bool &oAnyFronHasPropagated);
    void initiateFronts(bool &oAnyFronHasPropagated, TimeStep *tStep);
    bool hasPropagatingFronts();
    bool hasInitiationCriteria();

    /// Remove all enrichment items
    void clearEnrichmentItems();

    void appendEnrichmentItems(std :: vector< std :: unique_ptr< EnrichmentItem > > &iEIlist);

    void nucleateEnrichmentItems(bool &oNewItemsWereNucleated);
    bool hasNucleationCriteria();

    bool giveVtkDebug() const { return mDebugVTK; }
    void setVtkDebug(bool iDebug) { mDebugVTK = iDebug; }

    void updateNodeEnrichmentItemMap();

    const std :: vector< int > &giveNodeEnrichmentItemIndices(int iNodeIndex) const { return mNodeEnrichmentItemIndices [ iNodeIndex - 1 ]; }
    void giveElementEnrichmentItemIndices(std :: vector< int > &oElemEnrInd, int iElementIndex) const;

    const std :: vector< int > &giveMaterialModifyingEnrItemIndices() const { return mMaterialModifyingEnrItemIndices; }
};
} // end namespace oofem
#endif // xfemmanager_h
