/* $Header: /home/cvs/bp/oofem/sm/src/truss1d.h,v 1.8 2003/04/06 14:08:32 bp Exp $ */
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

//   *********************
//   *** CLASS TRUSS2D ***
//   *********************


#ifndef truss1d_h
#define truss1d_h


#include "structuralelement.h"
#include "gaussintegrationrule.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"
#include "zzerrorestimator.h"
#include "mmashapefunctprojection.h"
#include "huertaerrorestimator.h"

namespace oofem {

class Truss1d : public StructuralElement, public ZZNodalRecoveryModelInterface,
    public NodalAveragingRecoveryModelInterface, public SpatialLocalizerInterface,
    public DirectErrorIndicatorRCInterface,
    public EIPrimaryUnknownMapperInterface,
    public ZZErrorEstimatorInterface, public ZZRemeshingCriteriaInterface, public MMAShapeFunctProjectionInterface,
    public HuertaErrorEstimatorInterface, public HuertaRemeshingCriteriaInterface
{
    /*
     * This class implements a two-node truss bar element for two-dimensional
     * analysis.
     * DESCRIPTION :
     * A truss bar element is characterized by its 'length' and its 'pitch'. The
     * pitch is the angle in radians between the X-axis anf the axis of the
     * element (oriented node1 to node2).
     * The 'rotationMatrix' R is such that u{loc}=R*u{glob}.
     * Note: element is formulated in global c.s.
     * TASKS :
     * - calculating its Gauss points ;
     * - calculating its B,D,N matrices and dV ;
     * - expressing M,K,f,etc, in global axes. Methods like 'computeStiffness-
     *   Matrix' of class Element are here overloaded in order to account for
     *   rotational effects.
     */

protected:
    double length;
    // FloatMatrix*  rotationMatrix ;

public:
    Truss1d(int, Domain *);                         // constructor
    ~Truss1d()   { }                                // destructor

    // FloatArray*   ComputeBodyLoadVectorAt (TimeStep*) ;
    void          computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void          computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
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

    virtual int            computeNumberOfDofs(EquationID ut) { return 2; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // characteristic length in gp (for some material models)
    double        giveCharacteristicLenght(GaussPoint *, const FloatArray &)
    { return this->giveLength(); }

    double        computeVolumeAround(GaussPoint *);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void          drawScalar(oofegGraphicContext &context);
#endif
    //
    // definition & identification
    //
    const char *giveClassName() const { return "Truss1d"; }
    classType            giveClassID() const { return Truss1dClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

    /**
     * @name The element interface required by ZZNodalRecoveryModel
     */
    //@{
    /**
     * Computes the element contribution to \f$\int_\Omega N^TDBr\;d\Omega\f$.
     * The size of answer should be recordSize*numberofDofManagers.
     * @param answer contains the result
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    //void ZZNodalRecoveryMI_computeNValProduct (FloatArray& answer, InternalStateType type, TimeStep* tStep);
    /**
     * Computes the element contribution to \f$\int_\Omega N^TN\;d\Omega\f$ term.
     * The size of answer should be [recordSize*numberofDofManagers].
     * @param answer contain diagonalized result
     */
    //void ZZNodalRecoveryMI_computeNNMatrix (FloatArray& answer, InternalStateType type);
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
     * @name The element interface required by SPRNodalRecoveryModelInterface
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
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize() { return DirectErrorIndicatorRCI_giveCharacteristicSize(); }
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


protected:
    void          computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void          computeNmatrixAt(GaussPoint *, FloatMatrix &);
    //  int           computeGtoNRotationMatrix (FloatMatrix&);
    void          computeGaussPoints();
    integrationDomain  giveIntegrationDomain() { return _Line; }

    double        giveLength();
    double        givePitch();
    int           giveApproxOrder() { return 1; }
};

} // end namespace oofem
#endif // truss1d_h
