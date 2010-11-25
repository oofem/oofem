/* $Header: /home/cvs/bp/oofem/oofemlib/src/element.h,v 1.27 2003/04/06 14:08:24 bp Exp $ */
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

//   **********************
//   *** CLASS CCTPlate ***
//   **********************

#ifndef cct_h
#define cct_h

#include "nlstructuralelement.h"
#include "layeredcrosssection.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

namespace oofem {
class CCTPlate : public NLStructuralElement,
    public LayeredCrossSectionInterface, public ZZNodalRecoveryModelInterface,
    public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface
{
    /*
     * This class implements an triangular three-node  plate CCT
     * finite element. Each node has 3 degrees of freedom.
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
    CCTPlate(int, Domain *);                      // constructor
    ~CCTPlate() { }                               // destructor

protected:
    integrationDomain giveIntegrationDomain() { return _Triangle; }
    MaterialMode          giveMaterialMode()  { return _2dPlate; }
    void computeGaussPoints();
    void computeBodyLoadVectorAt(FloatArray &answer, Load *, TimeStep *, ValueModeType mode);
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *, FloatMatrix &);
    //int  computeGtoNRotationMatrix (FloatMatrix&);
    //void computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);

    virtual double giveArea();
    virtual void   giveNodeCoordinates(double &x1, double &x2, double &x3,
                                       double &y1, double &y2, double &y3,
                                       double *z = NULL);

public:
    //
    // definition & identification
    //
    const char *giveClassName() const { return "CCTPlate"; }
    classType    giveClassID()   const { return CCTPlateClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

    virtual int  computeNumberOfDofs(EquationID ut) { return 9; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // midPlaneNormal computation
    virtual FloatArray *ComputeMidPlaneNormal(GaussPoint *);

    // characteristic length in gp (for some material models)
    double giveCharacteristicLenght(GaussPoint *, const FloatArray &);
    double computeVolumeAround(GaussPoint *);

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }

    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);

    /**
     * Computes the global coordinates from given element's local coordinates.
     * Required by nonlocal material models.
     * @returns nonzero if successful
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    /**
     * Computes the element local (iso) coordinates from given global coordinates.
     * @returns nonzero if successful (if point is inside element); zero otherwise
     */
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    /**
     * Returns the integration point corresponding value in FULL form.
     * @param answer contain corresponding ip value, zero sized if not available.
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);


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
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type);
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
    int  SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    int SPRNodalRecoveryMI_giveNumberOfIP() { return 1; }
    //void SPRNodalRecoveryMI_giveIPValue (FloatArray& answer, int ipNum, InternalStateType type);
    void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    SPRPatchType SPRNodalRecoveryMI_givePatchType();
    //@}

    //
    // layered cross section support functions
    //
    void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                    GaussPoint *slaveGp, TimeStep *tStep);

    //
    // io routines
    //
#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void          drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    virtual void  drawScalar(oofegGraphicContext &context);
    //void          drawInternalState (oofegGraphicContext&);
#endif
};
} // end namespace oofem
#endif // cct_h
