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

/*
 * 	File: q9planstrss.h
 *
 *	Description: 9-node element.
 *
 *  Created on: May 22, 2013
 *  Author: Erik Svenning
 */

#ifndef Q9PLANSTRSS_H_
#define Q9PLANSTRSS_H_

#include "nlstructuralelement.h"
#include "fei2dquadbiquad.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_Q9PlaneStress2d_Name "q9planestress2d"

namespace oofem {

class Q9PlaneStress2d : public NLStructuralElement, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    int numberOfGaussPoints;
    static FEI2dQuadBiQuad interpolation;

public:
    Q9PlaneStress2d(int n, Domain *d);
    virtual ~Q9PlaneStress2d() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int computeNumberOfDofs(EquationID ut) { return 18; }

    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Q9PlaneStress2d_Name; }
    virtual const char *giveClassName() const { return "Q9PlaneStress2d"; }
    virtual classType giveClassID() const { return Q9PlaneStress2dClass; }
    virtual FEInterpolation *giveInterpolation() const { return & interpolation; }
    virtual MaterialMode giveMaterialMode() { return _PlaneStress; }

    virtual int testElementExtension(ElementExtension ext) { return ( ext == Element_EdgeLoadSupport ); }

    virtual Interface *giveInterface(InterfaceType it);

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);

    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }
    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void drawScalar(oofegGraphicContext &context);
    //void drawInternalState(DrawMode mode);
#endif

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();

    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);

};
} // end namespace oofem
#endif // qplanstrss_h

