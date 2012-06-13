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

#ifndef qtruss1d_h
#define qtruss1d_h

#include "fei1dquad.h"
#include "structuralelement.h"
#include "gaussintegrationrule.h"

namespace oofem {
/**
 * This class implements a three-node truss bar element for one-dimensional
 * analysis.
 */
class QTruss1d : public StructuralElement
{
protected:
    static FEI1dQuad interpolation;
    int numberOfGaussPoints;

public:
    QTruss1d(int n, Domain *d);
    virtual ~QTruss1d() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int computeNumberOfDofs(EquationID ut) { return 3; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
    { return this->computeLength(); }

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    // definition & identification
    virtual const char *giveClassName() const { return "QTruss1d"; }
    virtual classType giveClassID() const { return QTruss1dClass; }

    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_2; }
    virtual integrationDomain giveIntegrationDomain() { return _Line; }
    virtual MaterialMode giveMaterialMode() { return _1dMat; }

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();

    virtual int giveApproxOrder() { return 2; }
};
} // end namespace oofem
#endif // truss1d_h
