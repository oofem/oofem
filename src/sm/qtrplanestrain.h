/* $Header: /home/cvs/bp/oofem/sm/src/qtrplanestrain.h,v 1.5.4.1 2011/05/06 12:03:47 bp Exp $ */
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
//   *** CLASS QTrPlaneStrain ***
//   **************************

#ifndef qtrplanestrain_h
#define qtrplanestrain_h


#include "structuralelement.h"
#include "fei2dtrquad.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"

#include "mathfem.h"

namespace oofem {
class QTrPlaneStrain : public StructuralElement, public SpatialLocalizerInterface,
    public SPRNodalRecoveryModelInterface, public ZZNodalRecoveryModelInterface,
    public DirectErrorIndicatorRCInterface, public EIPrimaryUnknownMapperInterface
{
    /*
     * This class implements an triangular three-node  plane-
     * stress elasticity finite element. Each node has 2 degrees of freedom.
     *
     * DESCRIPTION :
     *
     * TASKS :
     *
     * - calculating its B,D,N matrices and dV.
     */

protected:
    static FEI2dTrQuad interpolation;
    int numberOfGaussPoints;

public:
    QTrPlaneStrain(int, Domain *);     // constructor
    ~QTrPlaneStrain()  { }             // destructor

    // FloatMatrix*    ComputeConstitutiveMatrixAt (GaussPoint *) ;
    // FloatArray*     ComputeBodyLoadVectorAt (TimeStep*) ;

    virtual int  computeNumberOfDofs(EquationID ut) { return 12; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    // characteristic length in gp (for some material models)
    double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane) {
        return this->giveLenghtInDir(normalToCrackPlane) / sqrt( ( double ) this->numberOfGaussPoints );
    }
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

    double             computeVolumeAround(GaussPoint *);

    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 0 : 0 ); }
    //int    hasEdgeLoadSupport () {return 0;}

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void  drawScalar(oofegGraphicContext &context);
    virtual void  drawSpecial(oofegGraphicContext &);
#endif
    //
    // definition & identification
    //
    const char *giveClassName() const { return "QTrPlaneStrain"; }
    classType       giveClassID()          const { return QTrPlaneStrainClass; }
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }
    IRResultType initializeFrom(InputRecord *ir);

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
     * @name The element interface required by ZZNodalRecoveryModel
     * currently used but problems in lumping procedure
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
    Element* ZZNodalRecoveryMI_giveElement () {return this;}
    /**
     * Evaluates N matrix (interpolation estimated stress matrix).
     */
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx (FloatMatrix& answer, GaussPoint* aGaussPoint, InternalStateType type);
    //@}

    /**
     * @name The element interface required by SPRNodalRecoveryModelInterface
     */
    //@{
    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    int SPRNodalRecoveryMI_giveNumberOfIP();
    //void SPRNodalRecoveryMI_giveIPValue (FloatArray& answer, int ipNum, InternalStateType type);
    void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    SPRPatchType SPRNodalRecoveryMI_givePatchType();
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

    integrationDomain  giveIntegrationDomain() { return _Triangle; }
    MaterialMode          giveMaterialMode()  { return _PlaneStrain; }


protected:
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *, FloatMatrix &);
    void computeGaussPoints();
    int giveApproxOrder() { return 2; }
    int giveNumberOfIPForMassMtrxIntegration() { return 4; }

    void computeDerivativeKsi(FloatArray &, double, double);
    void computeDerivativeEta(FloatArray &, double, double);
    void computeJacobianMatrixAt(FloatMatrix &, GaussPoint *);
};
} // end namespace oofem
#endif // qtrplanestrain_h
