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

#ifndef qtruss1d_h
#define qtruss1d_h

#include "../sm/Elements/nlstructuralelement.h"

#define _IFT_QTruss1d_Name "qtruss1d"

namespace oofem {
class FEI1dQuad;

/**
 * This class implements a three-node truss bar element for one-dimensional
 * analysis.
 */
class QTruss1d : public NLStructuralElement
{
protected:
    static FEI1dQuad interpolation;

public:
    QTruss1d(int n, Domain * d);
    virtual ~QTruss1d() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int computeNumberOfDofs() { return 3; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    virtual double giveCharacteristicLength(const FloatArray &normalToCrackPlane)
    { return this->computeLength(); }

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QTruss1d_Name; }
    virtual const char *giveClassName() const { return "QTruss1d"; }

    virtual MaterialMode giveMaterialMode() { return _1dMat; }
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();
    
};
} // end namespace oofem
#endif // truss1d_h
