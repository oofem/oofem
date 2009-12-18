/* $Header: /home/cvs/bp/oofem/sm/src/trplanestrain.h,v 1.8 2003/04/06 14:08:32 bp Exp $ */
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

#ifndef trplanstrain_h
#define trplanstrain_h

#include "nlstructuralelement.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"

#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"
#include "zzerrorestimator.h"
#include "mmashapefunctprojection.h"
#include "huertaerrorestimator.h"

namespace oofem {

class TrPlaneStrain : public StructuralElement, public ZZNodalRecoveryModelInterface,
    public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
    public SpatialLocalizerInterface
    , public DirectErrorIndicatorRCInterface,
    public EIPrimaryUnknownMapperInterface,
    public ZZErrorEstimatorInterface, public ZZRemeshingCriteriaInterface,
    public MMAShapeFunctProjectionInterface,
    public HuertaErrorEstimatorInterface, public HuertaRemeshingCriteriaInterface
{
    /*
     * This class implements an triangular three-node  plane-
     * strain elasticity finite element. Each node has 2 degrees of freedom.
     * DESCRIPTION :
     *
     * TASKS :
     *
     * - calculating its B,D,N matrices and dV.
     */

protected:
    double area;
    int numberOfGaussPoints;

public:
    TrPlaneStrain(int, Domain *);                        // constructor
    ~TrPlaneStrain()  { }                                // destructor

    virtual int            computeNumberOfDofs(EquationID ut) { return 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;


    // characteristic length in gp (for some material models)
    double        giveCharacteristicLenght(GaussPoint *, const FloatArray &);

    /**
     * Computes the global coordinates from given element's local coordinates.
     * Required by nonlocal material models.
     * @returns nonzero if successful
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    /**
     * Computes the element local coordinates from given global coordinates.
     * @returns nonzero if successful (if point is inside element); zero otherwise
     */
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    double             computeVolumeAround(GaussPoint *);

    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);

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
     * @name The element interface required by NodalAveragingRecoveryModel
     */
    //@{
    /**
     * Computes the element value in given node.
     * @param answer contains the result
     * @param node element node number
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    /**
     * Computes the element value in given side.
     * @param answer contains the result
     * @param node element side number
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    /**
     * Returns the size of DofManger record required to hold recovered values for given mode.
     * @param type determines the type of internal variable to be recovered
     * @return size of DofManger record required to hold recovered values
     */
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
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
     * @name The element interface required by DirectErrorIndicatorRC Interface
     */
    //@{
    /*
     * Determines the characteristic size of element. This quantity is defined as follows:
     * For 1D it is the element length, for 2D it is the square root of element area.
     */
    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize();
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

    // ZZErrorEstimatorInterface
    virtual Element *ZZErrorEstimatorI_giveElement() { return this; }
    virtual void ZZErrorEstimatorI_computeEstimatedStressInterpolationMtrx(FloatMatrix &answer, GaussPoint *gp,
                                                                           InternalStateType type)
    { ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(answer, gp, type); }

    // HuertaErrorEstimatorInterface
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

    // ZZRemeshingCriteriaInterface
    virtual double ZZRemeshingCriteriaI_giveCharacteristicSize() { return DirectErrorIndicatorRCI_giveCharacteristicSize(); }
    virtual int ZZRemeshingCriteriaI_givePolynOrder() { return 1; };

    // HuertaRemeshingCriteriaInterface
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize() { return DirectErrorIndicatorRCI_giveCharacteristicSize(); };
    virtual int HuertaRemeshingCriteriaI_givePolynOrder() { return 1; };

    /**
     * @name The element interface required by MMAShapeFunctProjectionInterface
     */
    //@{
    /**
     * Interpolates the internal variables, using given nodal values
     * to given point (given by global coordinates) using element shape functions.
     * @param answer computed internal variable
     * @param coords global/local coordinates of point of interest, see param ct
     * @param ct determines type of coordinates system used
     * @param list container of nodal values
     * @param type internal variable type
     * @param tStep solution step
     */
    virtual void MMAShapeFunctProjectionInterface_interpolateIntVarAt(FloatArray &answer, FloatArray &coords,
                                                                      coordType ct, nodalValContainerType &list,
                                                                      InternalStateType type, TimeStep *tStep);
    //@}

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void  drawScalar(oofegGraphicContext &context);
    //void          drawInternalState (oofegGraphicContext&);
#endif
    //
    // definition & identification
    //
    const char *giveClassName() const { return "TrPlaneStrain"; }
    classType            giveClassID() const { return TrPlaneStrainClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }


protected:
    // edge load support
    void  computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *);
    void  giveEdgeDofMapping(IntArray &answer, int) const;
    double        computeEdgeVolumeAround(GaussPoint *, int);
    void          computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    int   computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *);
    void          computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    //  void          computeNLBMatrixAt(FloatMatrix& answer, GaussPoint*, int);

    void               computeNmatrixAt(GaussPoint *, FloatMatrix &);
    // give Transformation matrix from global coord. syst. to local coordinate system in nodes.
    // i.e. r(n)=T r(g).
    //  int                computeGtoNRotationMatrix (FloatMatrix&);
    void               computeGaussPoints();
    integrationDomain  giveIntegrationDomain() { return _Triangle; }

    virtual double     giveArea();
    virtual FloatArray *GivebCoeff();
    virtual FloatArray *GivecCoeff();
    int           giveApproxOrder() { return 1; }
    /**
     * Return desired number of integration points for consistent mass matrix
     * computation, if required.
     */
    int  giveNumberOfIPForMassMtrxIntegration() { return 1; }
};

} // end namespace oofem
#endif // trplanstrain_h
