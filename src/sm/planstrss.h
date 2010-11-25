/* $Header: /home/cvs/bp/oofem/sm/src/planstrss.h,v 1.6.4.1 2004/04/05 15:19:47 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   **************************
//   *** CLASS PLANE STRAIN ***
//   **************************

#ifndef planstrss_h
#define planstrss_h

#include "nlstructuralelement.h"
#include "fei2dquadlin.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"

#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"
#include "huertaerrorestimator.h"

namespace oofem {
/// Comment or uncomment the following line to force full or reduced integration
#define PlaneStress2d_reducedShearIntegration

class PlaneStress2d : public NLStructuralElement, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
    public SpatialLocalizerInterface
    , public DirectErrorIndicatorRCInterface, public EIPrimaryUnknownMapperInterface,
    public HuertaErrorEstimatorInterface, public HuertaRemeshingCriteriaInterface
{
    /*
     * This class implements an isoparametric four-node quadrilateral plane-
     * stress elasticity finite element. Each node has 2 degrees of freedom.
     *
     * DESCRIPTION :
     *
     * TASKS :
     *
     * - calculating its Gauss points ;
     * - calculating its B,D,N matrices and dV.
     */

protected:
    static FEI2dQuadLin interpolation;
    int numberOfGaussPoints;

public:
    PlaneStress2d(int, Domain *); // constructor
    ~PlaneStress2d();           // destructor

    virtual int            computeNumberOfDofs(EquationID ut) { return 8; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // characteristic length in gp (for some material models)
    double        giveCharacteristicLenght(GaussPoint *, const FloatArray &);


    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);

    //int    hasEdgeLoadSupport () {return 1;}
    double                computeVolumeAround(GaussPoint *);
    /**
     * Computes the global coordinates from given element's local coordinates.
     * Required by nonlocal material models.
     * @returns nonzero if successful
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    /**
     * Computes the element local coordinates from given global coordinates.
     * @returns nonzero if successful (point inside); zero otherwise
     */
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    // returns interpolation type
    FEInterpolation *giveInterpolation() { return & interpolation; } // rch
    /**
     * @name The element interface required by ZZNodalRecoveryModel
     */
    //@{
    /**
     * Returns the size of DofManger record required to hold recovered values for given mode.
     * @param type determines the type of internal variable to be recovered
     * @return size of DofManger record required to hold recovered values
     */
    int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    /**
     * Returns the corresponding element to interface
     */
    Element *ZZNodalRecoveryMI_giveElement() { return this; }
    /**
     * Evaluates N matrix (interpolation estimated stress matrix).
     */
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint,
                                                             InternalStateType type);
    //@}
    /**
     * @name The element interface required by SPRNodalRecoveryModelInterface
     */
    //@{
    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    int SPRNodalRecoveryMI_giveNumberOfIP();
    //void SPRNodalRecoveryMI_giveIPValue (FloatArray& answer, int ipNum, InternalStateType type);
    void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    SPRPatchType SPRNodalRecoveryMI_givePatchType();
    //@}
    /**
     * @name The element interface required by SpatialLocalizerInterface
     */
    //@{
    /// Returns reference to corresponding element
    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    /// Returns nonzero if given element contains given point
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    /// Returns distance of given point from element parametric center
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);
    //@}


    /**
     * @name The element interface required by HuertaErrorEstimatorInterface and HuertaRemeshingCriteriaInterface
     */
    //@{
    virtual Element *HuertaErrorEstimatorI_giveElement() { return this; }

    virtual void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode);
    void HuertaErrorEstimatorI_computeLocalCoords(FloatArray &answer, const FloatArray &coords)
    { computeLocalCoordinates(answer, coords); }
    void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
    { computeNmatrixAt(aGaussPoint, answer); }

    // HuertaRemeshingCriteriaInterface
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize() { return DirectErrorIndicatorRCI_giveCharacteristicSize(); }
    virtual int HuertaRemeshingCriteriaI_givePolynOrder() { return 1; };
    //@}


    /**
     * @name The element interface required by SDirectErrorIndicatorRCInterface
     */
    //@{
    /*
     * Determines the characteristic size of element. This quantity is defined as follows:
     * For 1D it is the element length, for 2D it is the square root of element area.
     */
    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize();
    //@}
    /**
     * @name The element interface required by EIPrimaryUnknownMapperInterface
     */
    //@{
    /**
     * Computes the element vector of primary unknowns at given point. Similar to computeVectorOf,
     * but the interpolation from element DOFs to given point using element shape function is done.
     * The method should work also for point outside the volume of element (adaptivity mapping).
     * @param u    Identifies mode of unknown (eg. total value or velocity of unknown).
     * @param stepN Time step, when vector of unknowns is requested.
     * @param coords global coordinates of point of interest
     * @param answer vector of unknowns.
     */
    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType u,
                                                                 TimeStep *stepN, const FloatArray &coords,
                                                                 FloatArray &answer);
    /**
     * Returns the dof meaning of element vector of primary unknowns.
     * @param answer contains values of DofIDItem type that identify physical meaning of DOFs
     */
    virtual void EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer);

    //@}


#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void  drawScalar(oofegGraphicContext &context);
    virtual void  drawSpecial(oofegGraphicContext &);
    //     void          drawInternalState (oofegGraphicContext&);
#endif

    //
    // definition & identification
    //
    const char *giveClassName() const { return "PlaneStress2d"; }
    classType             giveClassID()          const { return PlaneStress2dClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_quad_1; }
    integrationDomain  giveIntegrationDomain() { return _Square; }
    MaterialMode          giveMaterialMode()  { return _PlaneStress; }

protected:
    // edge load support
    void  computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *);
    void  giveEdgeDofMapping(IntArray &answer, int) const;
    double        computeEdgeVolumeAround(GaussPoint *, int);
    void          computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    int   computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *);

    void                  computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void                  computeNLBMatrixAt(FloatMatrix &, GaussPoint *, int i);
    void                  computeNmatrixAt(GaussPoint *, FloatMatrix &);
    // give Transformation matrix from global coord. syst. to local coordinate system in nodes.
    // i.e. r(n)=T r(g)
    // int   computeGtoNRotationMatrix (FloatMatrix&);
    void                  computeGaussPoints();

    int           giveApproxOrder() { return 1; }
    int           giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // planstrss_h
