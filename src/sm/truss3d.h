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

#ifndef truss3d_h
#define truss3d_h

#include "nlstructuralelement.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "directerrorindicatorrc.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "fei3dlinelin.h"

#define _IFT_Truss3d_Name "truss3d"

namespace oofem {

/**
 * This class implements a two-node truss bar element for three-dimensional
 * analysis.
 */
class Truss3d : public NLStructuralElement,
    public DirectErrorIndicatorRCInterface,
    public ZZNodalRecoveryModelInterface,
    public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI3dLineLin interp;

public:
    Truss3d(int n, Domain *d);
    virtual ~Truss3d() { }

    virtual FEInterpolation *giveInterpolation() const { return &interp; }

    virtual double computeLength();

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { this->computeLumpedMassMatrix(answer, tStep); }
    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);

    virtual int computeNumberOfDofs() { return 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
    { return this->computeLength(); }

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int testElementExtension(ElementExtension ext) { return ( ext == Element_EdgeLoadSupport ); }

    virtual Interface *giveInterface(InterfaceType it);

    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);

    virtual double DirectErrorIndicatorRCI_giveCharacteristicSize() { return this->computeLength(); }

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Truss3d_Name; }
    virtual const char *giveClassName() const { return "Truss3d"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _1dMat; }

protected:
    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
    { computeGlobalCoordinates( answer, * ( gp->giveCoordinates() ) ); }
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual void computeGaussPoints();

    virtual int giveApproxOrder() { return 1; }
};
} // end namespace oofem
#endif // truss3d_h
