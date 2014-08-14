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
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef l4axisymm_h
#define l4axisymm_h

#include "Elements/nlstructuralelement.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"

#define _IFT_L4Axisymm_Name "l4axisymm"

namespace oofem {
class FEI2dQuadLin;

/**
 * This class implements an isoparametric four-node quadrilateral axisymmetric
 * finite element. Each node has 2 degrees of freedom.
 */
class L4Axisymm : public NLStructuralElement, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
public SpatialLocalizerInterface
{
protected:
    static FEI2dQuadLin interpolation;
    int numberOfGaussPoints, numberOfFiAndShGaussPoints;

public:
    L4Axisymm(int n, Domain * d);
    virtual ~L4Axisymm();


    virtual int computeNumberOfDofs() { return 8; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &crackToNormalPlane);

    virtual FEInterpolation *giveInterpolation() const;

    virtual Interface *giveInterface(InterfaceType it);

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_L4Axisymm_Name; }
    virtual const char *giveClassName() const { return "L4Axisymm"; }

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();


#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif

    virtual MaterialMode giveMaterialMode() { return _3dMat; }

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();

    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
};
} // end namespace oofem
#endif // l4axisymm_h
