/* $Header: /home/cvs/bp/oofem/sm/src/lspace.h,v 1.6 2003/05/19 13:04:00 bp Exp $ */
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

//   ************************************
//   ***      CLASS 3D Element        ***
//   ************************************

#ifndef lspace_h
#define lspace_h

#include "nlstructuralelement.h"
#include "fei3dhexalin.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"
#include "huertaerrorestimator.h"


class LSpace  : public NLStructuralElement, public ZZNodalRecoveryModelInterface,
    public SPRNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface,
    public SpatialLocalizerInterface
    , public EIPrimaryUnknownMapperInterface,
    public HuertaErrorEstimatorInterface, public HuertaRemeshingCriteriaInterface
{
    /*
     * This class implements a Linear 3d 8-node finite element for stress analysis.
     * Each node has 3 degrees of freedom.
     * DESCRIPTION :
     * One single additional attribute is needed for Gauss integration purpose :
     * 'jacobianMatrix'. This 3x3 matrix contains polynomials.
     * TASKS :
     * - calculating its Gauss points ;
     * - calculating its B,D,N matrices and dV.
     */

protected:
    //FloatMatrix*  jacobianMatrix ;
    static FEI3dHexaLin interpolation;
    int numberOfGaussPoints;

public:
    LSpace(int, Domain *);                       // constructor
    ~LSpace()  { }                                // destructor

    virtual int            computeNumberOfDofs(EquationID ut) { return 24; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    double             computeVolumeAround(GaussPoint *);
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

    /// characteristic length in gp (for some material models)
    double        giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);

    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);
    virtual int testElementExtension(ElementExtension ext)
    { return ( ( ( ext == Element_EdgeLoadSupport ) || ( ext == Element_SurfaceLoadSupport ) ) ? 1 : 0 ); }

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
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize();
    virtual int HuertaRemeshingCriteriaI_givePolynOrder() { return 1; };
    //@}

    //
    // definition & identification
    //
    const char *giveClassName() const { return "LSpace"; }
    classType             giveClassID()  const { return LSpaceClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_hexa_1; }


#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void          drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    virtual void  drawScalar(oofegGraphicContext &context);
    virtual void  drawSpecial(oofegGraphicContext &);
    void           drawTriad(FloatArray &, int isurf);
#endif



protected:
    void               computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void               computeNLBMatrixAt(FloatMatrix &, GaussPoint *, int i);
    void               computeNmatrixAt(GaussPoint *, FloatMatrix &);
    void       computeGaussPoints();
    integrationDomain  giveIntegrationDomain() { return _Cube; }
    int           giveApproxOrder() { return 1; }
    int           giveNumberOfIPForMassMtrxIntegration() { return 8; }

    /**
     * @name Edge load support
     */
    //@{
    void  computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *);
    void  giveEdgeDofMapping(IntArray &answer, int) const;
    double        computeEdgeVolumeAround(GaussPoint *, int);
    void          computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    int   computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *);

    //@}
    /**
     * @name Surface load support
     */
    //@{
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *);
    virtual void  giveSurfaceDofMapping(IntArray &answer, int) const;
    virtual IntegrationRule *GetSurfaceIntegrationRule(int);
    virtual double        computeSurfaceVolumeAround(GaussPoint *, int);
    virtual void        computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *, int);
    virtual int  computeLoadLSToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *);
    //@}
};

#endif // lspace_h
