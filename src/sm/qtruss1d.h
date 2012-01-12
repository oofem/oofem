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

#ifndef qtruss1d_h
#define qtruss1d_h

#include "fei1dquad.h"
#include "structuralelement.h"
#include "gaussintegrationrule.h"

namespace oofem {
/**
 * This class implements a two-node truss bar element for two-dimensional
 * analysis.
 *
 * A truss bar element is characterized by its 'length' and its 'pitch'. The
 * pitch is the angle in radians between the X-axis and the axis of the
 * element (oriented node1 to node2).
 * Note: element is formulated in global c.s.
 * Tasks:
 * - calculating its Gauss points ;
 * - calculating its B,D,N matrices and dV ;
 * - expressing M,K,f,etc, in global axes. Methods like 'computeStiffness-
 *   Matrix' of class Element are here overloaded in order to account for
 *   rotational effects.
 */
class QTruss1d : public StructuralElement
{
protected:
    double length;
    static FEI1dQuad interpolation;
    int numberOfGaussPoints;

public:
    QTruss1d(int n, Domain *d);
    ~QTruss1d() { }

    IRResultType initializeFrom(InputRecord *ir);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    virtual int computeNumberOfDofs(EquationID ut) { return 3; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
    { return this->giveLength(); }

    double computeVolumeAround(GaussPoint *gp);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    // definition & identification
    const char *giveClassName() const { return "QTruss1d"; }
    classType giveClassID() const { return QTruss1dClass; }

    Element_Geometry_Type giveGeometryType() const { return EGT_line_2; }
    integrationDomain giveIntegrationDomain() { return _Line; }
    MaterialMode giveMaterialMode() { return _1dMat; }

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeGaussPoints();

    double giveLength();
    double givePitch();
    int giveApproxOrder() { return 2; }
};
} // end namespace oofem
#endif // truss1d_h
