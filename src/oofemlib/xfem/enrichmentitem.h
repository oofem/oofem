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
#include "domain.h"
#include "floatmatrix.h"
#include "layeredcrosssection.h"
#include "dofiditem.h"
#include "tipinfo.h"

#include <memory>

#include "dofmanager.h"
#include <algorithm>

///@name Input fields for XFEM
//@{
#define _IFT_Inclusion_Name "inclusion"
#define _IFT_Inclusion_CrossSection "crosssection"

#define _IFT_EnrichmentItem_domains "enrichmentdomains"
#define _IFT_EnrichmentItem_domain "enrichmentdomain"
#define _IFT_EnrichmentItem_function "enrichmentfunction"
#define _IFT_EnrichmentItem_front "enrichmentfront"
#define _IFT_EnrichmentItem_propagationlaw "propagationlaw"

#define _IFT_Delamination_Name "delamination"
#define _IFT_Delamination_xiCoord "delaminationxicoord"
#define _IFT_Delamination_interfacenum "interfacenum"
#define _IFT_Delamination_csnum "csnum"
#define _IFT_Delamination_CohesiveZoneMaterial "czmaterial"
//#define _IFT_MultipleDelamination_Name "multipledelamination"
//@}

#define _IFT_Crack_Name "crack"

#define _IFT_EnrFrontDoNothing_Name "enrfrontdonothing"
#define _IFT_EnrFrontExtend_Name "enrfrontextend"
#define _IFT_EnrFrontLinearBranchFuncRadius_Name "enrfrontlinearbranchfuncradius"
#define _IFT_EnrFrontLinearBranchFuncRadius_Radius "radius"


namespace oofem {
template< class T >class AList;
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
    virtual void giveInputRecord(DynamicInputRecord &input) { OOFEM_ERROR("In EnrichmentItem::giveInputRecord(): This function must be called with DynamicDataReader as input.\n"); }
    virtual void appendInputRecords(DynamicDataReader &oDR);

    int instanciateYourself(DataReader *dr);
    virtual const char *giveClassName() const = 0;
    const IntArray *giveEnrichesDofsWithIdArray() const { return mpEnrichesDofsWithIdArray; }
    int giveNumberOfEnrDofs() const;

    // Spatial query
    bool isElementEnriched(const Element *element) const;
    inline bool isDofManEnriched(const DofManager &iDMan) const;
    int  giveNumDofManEnrichments(const DofManager &iDMan) const;

    // Returns true if the enrichment item assigns a different material to the Gauss point
    virtual bool isMaterialModified(GaussPoint &iGP, Element &iEl, CrossSection * &opCS) const;

    // Should update receiver geometry to the state reached at given time step.
    virtual void updateGeometry(FailureCriteriaStatus *fc, TimeStep *tStep) {};
    virtual void updateGeometry();
    virtual void propagateFronts();

    virtual bool hasPropagatingFronts() const { return mPropLawIndex != 0; }


    int giveStartOfDofIdPool() const { return this->startOfDofIdPool; };
    int giveEndOfDofIdPool() const { return this->endOfDofIdPool; };
    void computeDofManDofIdArray(IntArray &DofIdArray, DofManager *dMan); // list of id's a particular dof manager supports
    void giveEIDofIdArray(IntArray &answer) const; // list of id's for the enrichment dofs


    void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd = -1) const;
    void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const;
    void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, int iNodeInd) const;

    void evalLevelSetNormalInNode(double &oLevelSet, int iNodeInd) const { oLevelSet = mLevelSetNormalDir [ iNodeInd - 1 ]; }
    void evalLevelSetTangInNode(double &oLevelSet, int iNodeInd) const { oLevelSet = mLevelSetTangDir [ iNodeInd - 1 ]; }
    void evalNodeEnrMarkerInNode(double &oLevelSet, int iNodeInd) const { oLevelSet = mNodeEnrMarker [ iNodeInd - 1 ]; }

    bool levelSetChangesSignInEl(const IntArray &iElNodes) const;

    // By templating the function this way, we may choose if we want to pass iNodeInd as
    // an IntArray, a std::vector<int> or something else.
    // Any container that contains int and implements [] is legal.
    template< typename T >
    void interpLevelSet(double &oLevelSet, const FloatArray &iN, const T &iNodeInd) const;

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

    void createEnrichedDofs();

    virtual void computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, std :: vector< int > &oIntersectedEdgeInd, Element *element);


    // Return the coordinates of the tip in element iElIndex,
    // if the element contains a tip.
    bool giveElementTipCoord(FloatArray &oCoord, int iElIndex) const;

    // Help functions
    double calcXiZeroLevel(const double &iQ1, const double &iQ2) const;
    static void calcPolarCoord(double &oR, double &oTheta, const FloatArray &iOrigin, const FloatArray &iPos, const FloatArray &iN, const FloatArray &iT);

    PropagationLaw *givePropagationLaw() { return this->mpPropagationLaw; };
    bool hasPropagationLaw() { return this->mPropLawIndex != 0; };

    virtual void postInitialize();

protected:

    EnrichmentDomain *mpEnrichmentDomain;

    EnrichmentFunction *mpEnrichmentFunc;

    EnrichmentFront *mpEnrichmentFront;

    /// mEnrFrontIndex: nonzero if an enrichment front is present, zero otherwise.
    int mEnrFrontIndex;

    PropagationLaw *mpPropagationLaw;

    /// mPropLawIndex: nonzero if a propagation law is present, zero otherwise.
    int mPropLawIndex;

    int startOfDofIdPool; // points to the first available dofId number associated with the ei
    int endOfDofIdPool;

    /// Geometry associated with EnrichmentItem.
    IntArray *mpEnrichesDofsWithIdArray;


    // Level set for signed distance to the interface.
    //	The sign is determined by the interface normal direction.
    // This level set function is relevant for both open and closed interfaces.
    std :: vector< double >mLevelSetNormalDir;

    // Level set for signed distance along the interface.
    // Only relevant for open interfaces.
    std :: vector< double >mLevelSetTangDir;

    //	The sign is determined by the surface normal direction. Currently used
    //  to keep track of a delamination surface in a shell element
    std :: vector< double >mLevelSetSurfaceNormalDir;

    // Field with desired node enrichment types
    std :: vector< int >mNodeEnrMarker;

    // Indices of enriched nodes: this list is used to tell
    // if a given node is enriched.
    std :: vector< int >mEnrNodeIndices;

    bool mLevelSetsNeedUpdate;

    const double mLevelSetTol, mLevelSetTol2;
};

inline bool EnrichmentItem :: isDofManEnriched(const DofManager &iDMan) const
{
    return std :: binary_search( mEnrNodeIndices.begin(), mEnrNodeIndices.end(), iDMan.giveGlobalNumber() );
}

/** Inclusion. */
class OOFEM_EXPORT Inclusion : public EnrichmentItem
{
protected:
    CrossSection *mpCrossSection;
public:
    Inclusion(int n, XfemManager *xm, Domain *aDomain);
    virtual ~Inclusion();

    // Returns true if the enrichment item assigns a different material to the Gauss point
    virtual bool isMaterialModified(GaussPoint &iGP, Element &iEl, CrossSection * &opCS) const;


    virtual const char *giveClassName() const { return "Inclusion"; }
    virtual const char *giveInputRecordName() const { return _IFT_Inclusion_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    CrossSection *giveCrossSection() { return mpCrossSection; }
};


/** Delamination. */
class OOFEM_EXPORT Delamination : public EnrichmentItem
{
protected:
    Material *mat;  // Material for cohesive zone model
    int interfaceNum;
    int crossSectionNum;
    int matNum;
    double delamXiCoord;    // defines at what local xi-coord the delamination is defined
public:
    Delamination(int n, XfemManager *xm, Domain *aDomain);

    virtual const char *giveClassName() const { return "Delamination"; }
    virtual const char *giveInputRecordName() const { return _IFT_Delamination_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void appendInputRecords(DynamicDataReader &oDR);

    double giveDelamXiCoord() { return delamXiCoord; };
    //virtual Material *giveMaterial() { return mat; }
    virtual void updateGeometry(FailureCriteriaStatus *fc, TimeStep *tStep);
};





/**
 * Crack.
 * @author Erik Svenning
 */
class OOFEM_EXPORT Crack : public EnrichmentItem
{
public:
    Crack(int n, XfemManager *xm, Domain *aDomain);

    virtual const char *giveClassName() const { return "Crack"; }
    virtual const char *giveInputRecordName() const { return _IFT_Crack_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
};


/////////////////////////////////////////////////
// Function implementations

template< typename T >
void EnrichmentItem :: interpLevelSet(double &oLevelSet, const FloatArray &iN, const T &iNodeInd) const
{
    oLevelSet = 0.0;
    for ( int i = 1; i <= iN.giveSize(); i++ ) {
        oLevelSet += iN.at(i) * mLevelSetNormalDir [ iNodeInd [ i - 1 ] - 1 ];
    }
}

template< typename T >
void EnrichmentItem :: interpLevelSetTangential(double &oLevelSet, const FloatArray &iN, const T &iNodeInd) const
{
    oLevelSet = 0.0;
    for ( int i = 1; i <= iN.giveSize(); i++ ) {
        oLevelSet += iN.at(i) * mLevelSetTangDir [ iNodeInd [ i - 1 ] - 1 ];
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
            oGradLevelSet.at(j) += idNdX.at(i, j) * mLevelSetNormalDir [ iNodeInd [ i - 1 ] - 1 ];
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




/*
 * Class EnrichmentFront: describes the edge or tip of an XFEM enrichment.
 * The purpose is to add a different treatment of the front than the "interior"
 * enrichments. We may, e.g.
 *  - Apply branch functions at a crack tip for the element containing the crack tip.
 *  - Apply branch functions on nodes within a certain radius from the crack tip.
 *  - Exclude nodes touched by the front.
 *
 *  The desired behavior is obtained by choosing a suitable EnrichmentFront.
 *
 * @author Erik Svenning
 * August 2013
 */
class OOFEM_EXPORT EnrichmentFront
{
public:
    EnrichmentFront() {};
    virtual ~EnrichmentFront() {};

    /*
     *  MarkNodesAsFront:
     *  Intput:
     *  -ioNodeEnrMarker:   A vector with the same size as the number of nodes in the mesh
     *                      where the nodes corresponding to interior XFEM enrichments are
     *                      marked with 1, other entries are zero.
     *
     *  Output:
     *  -ioNodeEnrMarker:	Modifies the vector by marking tip nodes as 2, meaning that they
     *                      should get special treatment. May also modify the set of nodes
     *                      enriched by the interior enrichment.
     */
    virtual void MarkNodesAsFront(std :: vector< int > &ioNodeEnrMarker, XfemManager &ixFemMan, const std :: vector< double > &iLevelSetNormalDir, const std :: vector< double > &iLevelSetTangDir, const std :: vector< TipInfo > &iTipInfo) = 0;

    // The number of enrichment functions applied to tip nodes.
    virtual int  giveNumEnrichments(const DofManager &iDMan) const = 0;
    virtual int  giveMaxNumEnrichments() const = 0;


    // Evaluate the enrichment function and its derivative in front nodes.
    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const = 0;
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const = 0;
    virtual void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps) const = 0;

    virtual const char *giveClassName() const = 0;
    virtual const char *giveInputRecordName() const = 0;

    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
    virtual void giveInputRecord(DynamicInputRecord &input) = 0;

    virtual bool giveElementTipCoord(FloatArray &oCoord, int iElIndex) const;

protected:
    std :: vector< TipInfo >mTipInfo;

    /**
     * Keep record of the tips associated with an enriched node:
     * pair.first -> node index
     * pair.second-> tip indices
     */
    std :: vector< std :: pair< int, std :: vector< int > > >mNodeTipIndices;

    void addTipIndexToNode(int iNodeInd, int iTipInd); // Help function for updating mNodeTipIndices
    void giveNodeTipIndices(int iNodeInd, std :: vector< int > &oTipIndices) const;
};

class OOFEM_EXPORT EnrFrontDoNothing : public EnrichmentFront
{
public:
    EnrFrontDoNothing() {};
    virtual ~EnrFrontDoNothing() {};

    virtual void MarkNodesAsFront(std :: vector< int > &ioNodeEnrMarker, XfemManager &ixFemMan, const std :: vector< double > &iLevelSetNormalDir, const std :: vector< double > &iLevelSetTangDir, const std :: vector< TipInfo > &iTipInfo) { /*printf("Entering EnrFrontDoNothing::MarkNodesAsFront().\n");*/ }

    // No special tip enrichments are applied with this model.
    virtual int  giveNumEnrichments(const DofManager &iDMan) const { return 0; }
    virtual int  giveMaxNumEnrichments() const { return 0; }

    // Evaluate the enrichment function and its derivative in front nodes.
    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const {};
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const {};
    virtual void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps) const {};

    virtual const char *giveClassName() const { return "EnrFrontDoNothing"; }
    virtual const char *giveInputRecordName() const { return _IFT_EnrFrontDoNothing_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    virtual void giveInputRecord(DynamicInputRecord &input);
};

class OOFEM_EXPORT EnrFrontExtend : public EnrichmentFront
{
public:
    EnrFrontExtend() {};
    virtual ~EnrFrontExtend() {};

    virtual void MarkNodesAsFront(std :: vector< int > &ioNodeEnrMarker, XfemManager &ixFemMan, const std :: vector< double > &iLevelSetNormalDir, const std :: vector< double > &iLevelSetTangDir, const std :: vector< TipInfo > &iTipInfo);

    // No special tip enrichments are applied with this model,
    // it only modifies the set of nodes subject to bulk enrichment.
    virtual int  giveNumEnrichments(const DofManager &iDMan) const { return 0; }
    virtual int  giveMaxNumEnrichments() const { return 0; }

    // Evaluate the enrichment function and its derivative in front nodes.
    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const {};
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const {};
    virtual void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps) const {};


    virtual const char *giveClassName() const { return "EnrFrontExtend"; }
    virtual const char *giveInputRecordName() const { return _IFT_EnrFrontExtend_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    virtual void giveInputRecord(DynamicInputRecord &input);
};

class OOFEM_EXPORT EnrFrontLinearBranchFuncRadius : public EnrichmentFront
{
public:
    EnrFrontLinearBranchFuncRadius();
    virtual ~EnrFrontLinearBranchFuncRadius();

    virtual void MarkNodesAsFront(std :: vector< int > &ioNodeEnrMarker, XfemManager &ixFemMan, const std :: vector< double > &iLevelSetNormalDir, const std :: vector< double > &iLevelSetTangDir, const std :: vector< TipInfo > &iTipInfo);

    virtual int  giveNumEnrichments(const DofManager &iDMan) const;
    virtual int  giveMaxNumEnrichments() const { return 4; }

    // Evaluate the enrichment function and its derivative in front nodes.
    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const;
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const;
    virtual void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps) const;

    virtual const char *giveClassName() const { return "EnrFrontLinearBranchFuncRadius"; }
    virtual const char *giveInputRecordName() const { return _IFT_EnrFrontLinearBranchFuncRadius_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

private:
    double mEnrichmentRadius;
    LinElBranchFunction *mpBranchFunc;
};
} // end namespace oofem



#endif  // enrichmentitem_h
