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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef cct_h
#define cct_h

#include "nlstructuralelement.h"
#include "layeredcrosssection.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

namespace oofem {

/**
 * This class implements an triangular three-node plate CCT finite element.
 * Each node has 3 degrees of freedom.
 *
 * Tasks:
 * - calculating its B,D,N matrices and dV.
 */
class CCTPlate : public NLStructuralElement,
    public LayeredCrossSectionInterface, public ZZNodalRecoveryModelInterface,
    public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface
{
protected:
    double area;
    int numberOfGaussPoints;

public:
    CCTPlate(int n, Domain *d);
    ~CCTPlate() { }

protected:
    integrationDomain giveIntegrationDomain() { return _Triangle; }
    MaterialMode giveMaterialMode()  { return _2dPlate; }
    void computeGaussPoints();
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    //int computeGtoNRotationMatrix (FloatMatrix&);
    //void computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);

    virtual double giveArea();
    virtual void giveNodeCoordinates(double &x1, double &x2, double &x3,
                                     double &y1, double &y2, double &y3,
                                     double *z = NULL);

public:
    // definition & identification
    const char *giveClassName() const { return "CCTPlate"; }
    classType giveClassID() const { return CCTPlateClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

    virtual int computeNumberOfDofs(EquationID ut) { return 9; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    virtual void computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp);

    double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);
    double computeVolumeAround(GaussPoint *gp);

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }

    Interface *giveInterface(InterfaceType it);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);

    int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    Element *ZZNodalRecoveryMI_giveElement() { return this; }
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type);

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    int  SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    int SPRNodalRecoveryMI_giveNumberOfIP() { return 1; }
    //void SPRNodalRecoveryMI_giveIPValue(FloatArray& answer, int ipNum, InternalStateType type);
    void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    SPRPatchType SPRNodalRecoveryMI_givePatchType();

    // layered cross section support functions
    void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                    GaussPoint *slaveGp, TimeStep *tStep);

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    virtual void drawScalar(oofegGraphicContext &context);
    //void drawInternalState(oofegGraphicContext&);
#endif
};
} // end namespace oofem
#endif // cct_h
