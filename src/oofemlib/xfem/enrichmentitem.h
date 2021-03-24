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
#include "dofiditem.h"
#include "tipinfo.h"
#include "intarray.h"
#include "dofmanager.h"
#include "xfem/enrichmentfronts/enrichmentfront.h"
#include "xfem/enrichmentfunction.h"
#include "error.h"

#include <vector>
#include <unordered_map>

///@name Input fields for XFEM
//@{

#define _IFT_EnrichmentItem_domains "enrichmentdomains"
#define _IFT_EnrichmentItem_domain "enrichmentdomain"
#define _IFT_EnrichmentItem_function "enrichmentfunction"
#define _IFT_EnrichmentItem_front "enrichmentfront"
#define _IFT_EnrichmentItem_propagationlaw "propagationlaw"

#define _IFT_EnrichmentItem_inheritbc "inheritbc"
#define _IFT_EnrichmentItem_inheritorderedbc "inheritorderedbc"

//@}

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
class Element;
class CrossSection;
class Node;
class FloatMatrix;
class DofManager;
class DataReader;

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

    void initializeFrom(InputRecord &ir) override;
    virtual int instanciateYourself(DataReader &dr) = 0;
    virtual void postInitialize() = 0;

    /**
     * Note the special treatment here, the "normal" syntax
     * would be giveInputRecord(DynamicInputRecord &input).
     * Passing the entire DataReader instead allows us
     * to have separate InputRecords for the
     * EnrichmentDomain, EnrichmentFront and PropagationLaw
     * without have to keep track of them globally.
     */
    void giveInputRecord(DynamicInputRecord &input) override
    { OOFEM_ERROR("This function must be called with DynamicDataReader as input."); }
    virtual void appendInputRecords(DynamicDataReader &oDR) = 0;

    const IntArray *giveEnrichesDofsWithIdArray() const { return & mpEnrichesDofsWithIdArray; }
    int giveNumberOfEnrDofs() const;

    virtual void writeVtkDebug() const {};

    // Spatial query
    bool isElementEnriched(const Element *element) const;
    inline bool isDofManEnriched(const DofManager &iDMan) const;
    int  giveNumDofManEnrichments(const DofManager &iDMan) const;

    // Returns true if the enrichment item can assign
    // a different material to any Gauss point.
    inline virtual bool canModifyMaterial() const { return false; }

    // Returns true if the enrichment item assigns a different material to the Gauss point
    virtual bool isMaterialModified(GaussPoint &iGP, Element &iEl, CrossSection * &opCS) const;


    // Should update receiver geometry to the state reached at given time step.
    virtual void updateGeometry(FailureCriteriaStatus *fc, TimeStep *tStep) { }
    virtual void updateGeometry(TimeStep *tStep) { }
    virtual void updateGeometry() = 0;
    virtual void propagateFronts(bool &oFrontsHavePropagated) = 0;

    virtual bool hasPropagatingFronts() const;
    virtual bool hasInitiationCriteria() { return false; }


    int giveStartOfDofIdPool() const { return this->startOfDofIdPool; }
    int giveEndOfDofIdPool() const { return this->endOfDofIdPool; }
    virtual int giveDofPoolSize() const;
    /**
     * Compute Id's of enriched dofs for a given DofManager.
     */
    virtual void computeEnrichedDofManDofIdArray(IntArray &oDofIdArray, DofManager &iDMan);

    virtual void giveEIDofIdArray(IntArray &answer) const; // list of id's for the enrichment dofs
    virtual void givePotentialEIDofIdArray(IntArray &answer) const; // List of potential IDs for enrichment

    virtual void evaluateEnrFuncInNode(std :: vector< double > &oEnrFunc, const Node &iNode) const = 0;

    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl) const = 0;
    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl, const FloatArray &iN, const IntArray &iElNodes) const = 0;

    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl) const = 0;
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl, const FloatArray &iN, const FloatMatrix &idNdX, const IntArray &iElNodes) const = 0;

    bool evalLevelSetNormalInNode(double &oLevelSet, int iNodeInd, const FloatArray &iGlobalCoord) const;
    bool evalLevelSetTangInNode(double &oLevelSet, int iNodeInd, const FloatArray &iGlobalCoord) const;
    bool evalNodeEnrMarkerInNode(double &oNodeEnrMarker, int iNodeInd) const;

protected:
    /**
     * Evaluate the normal direction level set in the
     * point iGlobalCoord. To improve performance, basis
     * function values corresponding to that point should also be provided.
     */
    virtual void evalLevelSetNormal(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const = 0;

    /**
     * Evaluate the tangential direction level set in the
     * point iGlobalCoord. To improve performance, basis
     * function values corresponding to that point should also be provided.
     */
    virtual void evalLevelSetTangential(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const = 0;

    /**
     * Evaluate the gradient of the normal direction level set in the
     * point iGlobalCoord. To improve performance, basis
     * function values corresponding to that point should also be provided.
     */
    virtual void evalGradLevelSetNormal(FloatArray &oGradLevelSet, const FloatArray &iGlobalCoord, const FloatMatrix &idNdX, const IntArray &iNodeInd) const = 0;


    // Level set routines
    virtual void updateNodeEnrMarker(XfemManager &ixFemMan) = 0;

public:

    virtual void createEnrichedDofs();

    // Return the coordinates of the tip in element iElIndex,
    // if the element contains a tip.
    virtual bool giveElementTipCoord(FloatArray &oCoord, double &oArcPos, Element &iEl, const FloatArray &iElCenter) const = 0;

    // Help functions
    static double calcXiZeroLevel(const double &iQ1, const double &iQ2);
    static void calcPolarCoord(double &oR, double &oTheta, const FloatArray &iOrigin, const FloatArray &iPos, const FloatArray &iN, const FloatArray &iT, const EfInput &iEfInput, bool iFlipTangent);

    PropagationLaw *givePropagationLaw() { return this->mpPropagationLaw.get(); }
    void setPropagationLaw(std::unique_ptr<PropagationLaw> ipPropagationLaw);
    bool hasPropagationLaw() { return this->mPropLawIndex != 0; }


    virtual void callGnuplotExportModule(GnuplotExportModule &iExpMod, TimeStep *tStep);


    const std :: unordered_map< int, NodeEnrichmentType > &giveEnrNodeMap() const { return mNodeEnrMarkerMap; }

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius) = 0;

    EnrichmentFront *giveEnrichmentFrontStart() { return mpEnrichmentFrontStart.get(); }
    void setEnrichmentFrontStart(std::unique_ptr<EnrichmentFront> ipEnrichmentFrontStart, bool iDeleteOld = true);

    EnrichmentFront *giveEnrichmentFrontEnd() { return mpEnrichmentFrontEnd.get(); }
    void setEnrichmentFrontEnd(std::unique_ptr<EnrichmentFront> ipEnrichmentFrontEnd, bool iDeleteOld = true);

    bool tipIsTouchingEI(const TipInfo &iTipInfo);

    void setEnrichmentFunction(std::unique_ptr<EnrichmentFunction> ipEnrichmentFunc) { mpEnrichmentFunc = std::move(ipEnrichmentFunc); }

protected:

    std::unique_ptr<EnrichmentFunction> mpEnrichmentFunc;

    std::unique_ptr<EnrichmentFront> mpEnrichmentFrontStart, mpEnrichmentFrontEnd;

    /// mEnrFrontIndex: nonzero if an enrichment front is present, zero otherwise.
    int mEnrFrontIndex;

    std::unique_ptr<PropagationLaw> mpPropagationLaw;

    /// mPropLawIndex: nonzero if a propagation law is present, zero otherwise.
    int mPropLawIndex;

    /**
     * If newly created enriched dofs should inherit boundary conditions
     * from the node they are introduced in. Default is false, i.e.
     * XFEM dofs are free by default. Two alternatives exists:
     * mInheritBoundaryConditions takes the first Dirichlet BC it finds in the 
     * node. Therefore, we may get in trouble if the node has different Dirichlet 
     * BCs for different dofs.
     * mInheritOrderedBoundaryConditions assumes that the enriched dofs are of the
     * same type and in the same order as the original dofs and uses the same BC
     * on the enriched dofs
     * NB: These routines basically only works for zero Dirichlet BC, since enriched
     * dofs often are related to the existing dofs. If nonzero BC are prescibed this 
     * will be a problem.
     */
    bool mInheritBoundaryConditions;
    bool mInheritOrderedBoundaryConditions;

    int startOfDofIdPool; // points to the first available dofId number associated with the ei
    int endOfDofIdPool;

    /// Geometry associated with EnrichmentItem.
    IntArray mpEnrichesDofsWithIdArray;


    // Level set for signed distance to the interface.
    // The sign is determined by the interface normal direction.
    // This level set function is relevant for both open and closed interfaces.
    std :: unordered_map< int, double >mLevelSetNormalDirMap;

    // Level set for signed distance along the interface.
    // Only relevant for open interfaces.
    std :: unordered_map< int, double >mLevelSetTangDirMap;


    // Field with desired node enrichment types
    std :: unordered_map< int, NodeEnrichmentType >mNodeEnrMarkerMap;

    // Enrichment dof IDs used by the enrichment item.
    IntArray mEIDofIdArray;

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
} // end namespace oofem



#endif  // enrichmentitem_h
