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

#ifndef qtrplstr_h
#define qtrplstr_h

#include "structuralelement.h"
#include "fei2dtrquad.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"

#include "mathfem.h"

namespace oofem {
/**
 * This class implements a quadratic triangular 6-node plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
class QTrPlaneStress2d : public StructuralElement, public SpatialLocalizerInterface,
    public SPRNodalRecoveryModelInterface,
    public DirectErrorIndicatorRCInterface, public EIPrimaryUnknownMapperInterface
{
protected:
    static FEI2dTrQuad interpolation;
    int numberOfGaussPoints;

public:
    QTrPlaneStress2d(int n, Domain *d);
    ~QTrPlaneStress2d() { }

    virtual int  computeNumberOfDofs(EquationID ut) { return 12; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    double computeVolumeAround(GaussPoint *gp);

    Interface *giveInterface(InterfaceType it);

    FEInterpolation *giveInterpolation() { return & interpolation; }

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 0 : 0 ); }

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void  drawScalar(oofegGraphicContext &context);
    virtual void  drawSpecial(oofegGraphicContext &);
#endif

    // definition & identification
    const char *giveClassName() const { return "QTrPlaneStress2d"; }
    classType giveClassID() const { return QTrPlaneStress2dClass; }
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }
    IRResultType initializeFrom(InputRecord *ir);

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    int SPRNodalRecoveryMI_giveNumberOfIP();
    //void SPRNodalRecoveryMI_giveIPValue(FloatArray& answer, int ipNum, InternalStateType type);
    void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize();

    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType u,
                                                                 TimeStep *stepN, const FloatArray &coords,
                                                                 FloatArray &answer);
    virtual void EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer);

    integrationDomain  giveIntegrationDomain() { return _Triangle; }
    MaterialMode giveMaterialMode() { return _PlaneStress; }


protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeGaussPoints();
    int giveApproxOrder() { return 2; }
    int giveNumberOfIPForMassMtrxIntegration() { return 4; }

    void computeDerivativeKsi(FloatArray &, double, double);
    void computeDerivativeEta(FloatArray &, double, double);
    void computeJacobianMatrixAt(FloatMatrix &, GaussPoint *);
};
} // end namespace oofem
#endif // qtrplstr_h
