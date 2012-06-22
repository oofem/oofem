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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef axisymm3d_h
#define axisymm3d_h

#include "nlstructuralelement.h"
#include "fei2dtrlin.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"

namespace oofem {
/**
 * This class implements an triangular three-node finite element for axisymmetric continuum.
 * Each node has 2 degrees of freedom.
 *
 * Tasks:
 * - calculating its B,D,N matrices and dV.
 */
class Axisymm3d : public NLStructuralElement, public ZZNodalRecoveryModelInterface,
    public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
    public SpatialLocalizerInterface
{

protected:
    static FEI2dTrLin interpolation;

    int numberOfGaussPoints, numberOfFiAndShGaussPoints;
    double area;

public:
    Axisymm3d(int n, Domain *d);
    virtual ~Axisymm3d();

    virtual int computeNumberOfDofs(EquationID ut) { return 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;

    // characteristic length in gp (for some material models)
    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);
    virtual double giveArea();
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);
    virtual FEInterpolation *giveInterpolation() { return & interpolation; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    virtual void drawScalar(oofegGraphicContext &context);
    //void drawInternalState(oofegGraphicContext&);
#endif

    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);


    virtual const char *giveClassName() const { return "Axisymm3d"; }
    virtual classType giveClassID() const { return Axisymm3dClass; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual integrationDomain giveIntegrationDomain() { return _Triangle; }
    virtual MaterialMode giveMaterialMode() { return _3dMat; }

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();

    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &, int iEdge, GaussPoint *gp);
};
} // end namespace oofem
#endif // axisymm3d_h
