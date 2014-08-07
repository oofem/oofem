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

#ifndef enrichmentitem_h
#define enrichmentitem_h

#include "femcmpnn.h"
#include "tipinfo.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "dofmanager.h"
#include "xfem/enrichmentfronts/enrichmentfront.h"

#include <memory>
#include <algorithm>
#include <unordered_map>


///@name Input fields for EnrichmentItem
//@{
#define _IFT_EnrichmentItem_domains "enrichmentdomains"
#define _IFT_EnrichmentItem_domain "enrichmentdomain"
#define _IFT_EnrichmentItem_function "enrichmentfunction"
#define _IFT_EnrichmentItem_front "enrichmentfront"
#define _IFT_EnrichmentItem_propagationlaw "propagationlaw"
#define _IFT_EnrichmentItem_inheritbc "inheritbc"
//@}

#define _IFT_Inclusion_Name "inclusion"
#define _IFT_Inclusion_CrossSection "crosssection"


namespace oofem {
class BasicGeometry;
class EnrichmentFunction;
class EnrichmentDomain;
class FractureManager;
class FailureCriteriaStatus;
class EnrichmentDomain_BG;
class DofManList;
class WholeDomain;
class EnrichmentFront;
class LinElBranchFunction;
class PropagationLaw;
class DynamicDataReader;
class Triangle;
class GnuplotExportModule;
class GaussPoint;
class CrossSection;
class Element;


enum NodeEnrichmentType : int {
    NodeEnr_NONE = 0,
    NodeEnr_BULK = 1,
    NodeEnr_START_TIP = 2,
    NodeEnr_END_TIP = 3,
    NodeEnr_START_AND_END_TIP = 4
};

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
 * @author Erik Svenning
 */
class OOFEM_EXPORT EnrichmentItem : public FEMComponent
{
public:
    /// Constructor / destructor
    EnrichmentItem(int n, XfemManager *xm, Domain *aDomain);
    virtual ~EnrichmentItem();

    virtual IRResultType initializeFrom(InputRecord *ir);

    /**
     * Note the special treatment here, the "normal" syntax
     * would be giveInputRecord(DynamicInputRecord &input).
     * Passing the entire DataReader instead allows us
     * to have separate InputRecords for the
     * EnrichmentDomain, EnrichmentFront and PropagationLaw
     * without have to keep track of them globally.
     */
    virtual void giveInputRecord(DynamicInputRecord &input) { OOFEM_ERROR("This function must be called with DynamicDataReader as input."); }
    virtual void appendInputRecords(DynamicDataReader &oDR);

    int instanciateYourself(DataReader *dr);
    virtual const char *giveClassName() const = 0;
    const IntArray *giveEnrichesDofsWithIdArray() const { return & mpEnrichesDofsWithIdArray; }
    int giveNumberOfEnrDofs() const;

    void writeVtkDebug() const;

    // Spatial query
    bool isElementEnriched(const Element *element) const;
    inline bool isDofManEnriched(const DofManager &iDMan) const;
    int  giveNumDofManEnrichments(const DofManager &iDMan) const;
    int giveNumEnrichedDofs(const DofManager &iDMan) const;

    // Returns true if the enrichment item can assign
    // a different material to any Gauss point.
    inline virtual bool canModifyMaterial() const { return false; }

    // Returns true if the enrichment item assigns a different material to the Gauss point
    virtual bool isMaterialModified(GaussPoint &iGP, Element &iEl, CrossSection * &opCS) const;


    // Should update receiver geometry to the state reached at given time step.
    virtual void updateGeometry(FailureCriteriaStatus *fc, TimeStep *tStep) { };
    virtual void updateGeometry();
    virtual void propagateFronts();

    virtual bool hasPropagatingFronts() const;


    int giveStartOfDofIdPool() const { return this->startOfDofIdPool; };
    int giveEndOfDofIdPool() const { return this->endOfDofIdPool; };

    /**
     * Compute Id's of enriched dofs for a given DofManager.
     */
    virtual void computeEnrichedDofManDofIdArray(IntArray &oDofIdArray, DofManager &iDMan);

    void giveEIDofIdArray(IntArray &answer) const; // list of id's for the enrichment dofs


    void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd = -1) const;
    void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const;
    void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, int iNodeInd, GaussPoint &iGP, bool iGPLivesOnCurrentCrack) const;

    bool evalLevelSetNormalInNode(double &oLevelSet, int iNodeInd) const;
    bool evalLevelSetTangInNode(double &oLevelSet, int iNodeInd) const;
    bool evalNodeEnrMarkerInNode(double &oNodeEnrMarker, int iNodeInd) const;

    bool levelSetChangesSignInEl(const IntArray &iElNodes) const;

    void interpLevelSet(double &oLevelSet, const FloatArray &iGlobalCoord) const;

    // By templating the function this way, we may choose if we want to pass iNodeInd as
    // an IntArray, a std::vector<int> or something else.
    // Any container that contains int and implements [] is legal.
    template< typename T >
    inline void interpLevelSet(double &oLevelSet, const FloatArray &iN, const T &iNodeInd) const;

    template< typename T >
    void interpLevelSetTangential(double &oLevelSet, const FloatArray &iN, const T &iNodeInd) const;

    template< typename T >
    void interpGradLevelSet(FloatArray &oGradLevelSet, const FloatMatrix &idNdX, const T &iNodeInd) const;

    // JB - temporary
    template< typename T >
    void interpSurfaceLevelSet(double &oLevelSet, const FloatArray &iN, const T &iNodeInd, double iXi) const;
    void interpSurfaceLevelSet(double &oLevelSet, double iXi) const;

    // Level set routines
    bool giveLevelSetsNeedUpdate() const { return mLevelSetsNeedUpdate; }
    virtual void updateLevelSets(XfemManager &ixFemMan);
    virtual void updateNodeEnrMarker(XfemManager &ixFemMan, const EnrichmentDomain_BG &iEnrichmentDomain_BG);
    virtual void updateNodeEnrMarker(XfemManager &ixFemMan, const DofManList &iDofManList);
    virtual void updateNodeEnrMarker(XfemManager &ixFemMan, const WholeDomain &iWholeDomain);

    virtual void createEnrichedDofs();

    virtual void computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, std :: vector< int > &oIntersectedEdgeInd, Element *element, std :: vector< double > &oMinDistArcPos) const;
    virtual void computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, std :: vector< int > &oIntersectedEdgeInd, Element *element, const Triangle &iTri, std :: vector< double > &oMinDistArcPos) const;


    // Return the coordinates of the tip in element iElIndex,
    // if the element contains a tip.
    bool giveElementTipCoord(FloatArray &oCoord, double &oArcPos, int iElIndex, const FloatArray &iElCenter) const;
    bool giveElementTipCoord(FloatArray &oCoord, double &oArcPos, int iElIndex, const Triangle &iTri, const FloatArray &iElCenter) const;

    // Help functions
    static double calcXiZeroLevel(const double &iQ1, const double &iQ2);
    static void calcPolarCoord(double &oR, double &oTheta, const FloatArray &iOrigin, const FloatArray &iPos, const FloatArray &iN, const FloatArray &iT, const EfInput &iEfInput, bool iFlipTangent);

    PropagationLaw *givePropagationLaw() { return this->mpPropagationLaw; };
    bool hasPropagationLaw() { return this->mPropLawIndex != 0; };

    void giveSubPolygon(std :: vector< FloatArray > &oPoints, const double &iXiStart, const double &iXiEnd) const;

    virtual void callGnuplotExportModule(GnuplotExportModule &iExpMod);

    EnrichmentDomain *giveEnrichmentDomain() const { return mpEnrichmentDomain; }

    const std :: unordered_map< int, NodeEnrichmentType > &giveEnrNodeMap() const { return mNodeEnrMarkerMap; }

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius);

    EnrichmentFront *giveEnrichmentFrontStart() {return mpEnrichmentFrontStart;}
    void setEnrichmentFrontStart(EnrichmentFront *ipEnrichmentFrontStart);

    EnrichmentFront *giveEnrichmentFrontEnd() {return mpEnrichmentFrontEnd;}
    void setEnrichmentFrontEnd(EnrichmentFront *ipEnrichmentFrontEnd);

protected:

    EnrichmentDomain *mpEnrichmentDomain;

    EnrichmentFunction *mpEnrichmentFunc;

    EnrichmentFront *mpEnrichmentFrontStart, *mpEnrichmentFrontEnd;

    /// mEnrFrontIndex: nonzero if an enrichment front is present, zero otherwise.
    int mEnrFrontIndex;

    PropagationLaw *mpPropagationLaw;

    /// mPropLawIndex: nonzero if a propagation law is present, zero otherwise.
    int mPropLawIndex;

    /**
     * If newly created enriched dofs should inherit boundary conditions
     * from the node they are introduced in. Default is false, i.e.
     * XFEM dofs are free by default. Note: the routine takes the first
     * Dirichlet BC it finds in the node. Therefore, we may get in trouble
     * if the node has different Dirichlet BCs for different dofs.
     */
    bool mInheritBoundaryConditions;

    int startOfDofIdPool; // points to the first available dofId number associated with the ei
    int endOfDofIdPool;

    /// Geometry associated with EnrichmentItem.
    IntArray mpEnrichesDofsWithIdArray;


    // Level set for signed distance to the interface.
    //	The sign is determined by the interface normal direction.
    // This level set function is relevant for both open and closed interfaces.
    std :: unordered_map< int, double >mLevelSetNormalDirMap;

    // Level set for signed distance along the interface.
    // Only relevant for open interfaces.
    std :: unordered_map< int, double >mLevelSetTangDirMap;

    //	The sign is determined by the surface normal direction. Currently used
    //  to keep track of a delamination surface in a shell element
    std :: vector< double >mLevelSetSurfaceNormalDir;

    // Field with desired node enrichment types
    std :: unordered_map< int, NodeEnrichmentType >mNodeEnrMarkerMap;

    bool mLevelSetsNeedUpdate;

    static const double mLevelSetTol;
    static const double mLevelSetRelTol;
    const double mLevelSetTol2;
};

inline bool EnrichmentItem :: isDofManEnriched(const DofManager &iDMan) const
{
    auto res = mNodeEnrMarkerMap.find( iDMan.giveGlobalNumber() );
    return !( res == mNodeEnrMarkerMap.end() );
}

/** Inclusion. */
class OOFEM_EXPORT Inclusion : public EnrichmentItem
{
protected:
    CrossSection *mpCrossSection;
public:
    Inclusion(int n, XfemManager *xm, Domain *aDomain);
    virtual ~Inclusion();

    // Returns true if the enrichment item can assign
    // a different material to any Gauss point.
    inline virtual bool canModifyMaterial() const { return true; }

    // Returns true if the enrichment item assigns a different material to the Gauss point
    virtual bool isMaterialModified(GaussPoint &iGP, Element &iEl, CrossSection * &opCS) const;


    virtual const char *giveClassName() const { return "Inclusion"; }
    virtual const char *giveInputRecordName() const { return _IFT_Inclusion_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    CrossSection *giveCrossSection() { return mpCrossSection; }
};


/////////////////////////////////////////////////
// Function implementations

template< typename T >
inline void EnrichmentItem :: interpLevelSet(double &oLevelSet, const FloatArray &iN, const T &iNodeInd) const
{
    oLevelSet = 0.0;
    for ( int i = 1; i <= iN.giveSize(); i++ ) {
        double levelSetNode = 0.0;
        if ( evalLevelSetNormalInNode(levelSetNode, iNodeInd [ i - 1 ]) ) {
            oLevelSet += iN.at(i) * levelSetNode;
        }
    }
}

template< typename T >
void EnrichmentItem :: interpLevelSetTangential(double &oLevelSet, const FloatArray &iN, const T &iNodeInd) const
{
    oLevelSet = 0.0;
    for ( int i = 1; i <= iN.giveSize(); i++ ) {
        double levelSetNode = 0.0;
        if ( evalLevelSetTangInNode(levelSetNode, iNodeInd [ i - 1 ]) ) {
            oLevelSet += iN.at(i) * levelSetNode;
        }
    }
}

template< typename T >
void EnrichmentItem :: interpGradLevelSet(FloatArray &oGradLevelSet, const FloatMatrix &idNdX, const T &iNodeInd) const
{
    int dim = idNdX.giveNumberOfColumns();

    if ( oGradLevelSet.giveSize() != dim ) {
        oGradLevelSet.resize(dim);
    }

    oGradLevelSet.zero();

    for ( int i = 1; i <= idNdX.giveNumberOfRows(); i++ ) {
        for ( int j = 1; j <= dim; j++ ) {
            double levelSetNode = 0.0;
            if ( evalLevelSetNormalInNode(levelSetNode, iNodeInd [ i - 1 ]) ) {
                oGradLevelSet.at(j) += idNdX.at(i, j) * levelSetNode;
            }
        }
    }
}


// should be generalised - JB
template< typename T >
void EnrichmentItem :: interpSurfaceLevelSet(double &oLevelSet, const FloatArray &iN, const T &iNodeInd, double iXi) const
{
    //oLevelSet = 0.0;
    oLevelSet = iXi;
    for ( int i = 1; i <= iN.giveSize(); i++ ) {
        oLevelSet -= iN.at(i) * mLevelSetSurfaceNormalDir [ iNodeInd [ i - 1 ] - 1 ];
    }
}
} // end namespace oofem



#endif  // enrichmentitem_h
