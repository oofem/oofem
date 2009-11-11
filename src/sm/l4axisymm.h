/* $Header: /home/cvs/bp/oofem/sm/src/l4axisymm.h,v 1.4.4.1 2004/04/05 15:19:47 bp Exp $ */
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
// 9.5.1995
//
#ifndef l4axisymm_h
#define l4axisymm_h

#include "nlstructuralelement.h"
#include "fei2dquadlin.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"

class L4Axisymm : public NLStructuralElement, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
    public SpatialLocalizerInterface
{
    /*
     * This class implements an isoparametric four-node quadrilateral axisymme-
     * tric finite element. Each node has 2 degrees of freedom.
     *
     * DESCRIPTION :
     *
     * One single additional attribute is needed for Gauss integration purpose :
     * 'jacobianMatrix'. This 2x2 matrix contains polynomials.
     *
     * TASKS :
     *
     * - calculating its Gauss points ;
     * - calculating its B,D,N matrices and dV.
     */

protected:

    static FEI2dQuadLin interpolation;
    int numberOfGaussPoints, numberOfFiAndShGaussPoints;

public:

    L4Axisymm(int, Domain *);                    // constructor
    ~L4Axisymm();                                // destructor


    virtual int            computeNumberOfDofs(EquationID ut) { return 8; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    double                computeVolumeAround(GaussPoint *);
    void                  computeStrainVector(FloatArray &answer, GaussPoint *, TimeStep *);

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
    // characteristic length in gp (for some material models)
    double        giveCharacteristicLenght(GaussPoint *, const FloatArray &);

    //
    // definition & identification
    //
    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);
    MaterialMode          giveMaterialMode() { return _3dMat; }
    const char *giveClassName() const { return "L4Axisymm"; }
    classType     giveClassID() const { return L4AxisymmClass; }
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_quad_1; }

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
    //@}

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void          drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    virtual void  drawScalar(oofegGraphicContext &context);
    // void          drawInternalState (oofegGraphicContext&);

#endif

protected:
    void                  computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void                  computeNmatrixAt(GaussPoint *, FloatMatrix &);
    void                  computeGaussPoints();
    integrationDomain  giveIntegrationDomain() { return _Square; }

    // edge load support
    void  computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *);
    void  giveEdgeDofMapping(IntArray &answer, int) const;
    double        computeEdgeVolumeAround(GaussPoint *, int);
    void          computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    int   computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *);
};

#endif // l4axisymm_h
