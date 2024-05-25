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

#include "sm/Elements/nlstructuralelement.h"

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
    QTruss1d(int n, Domain *d);
    virtual ~QTruss1d() { }

    void initializeFrom(InputRecord &ir) override;

    int computeNumberOfDofs() override { return 3; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;

    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    double giveCharacteristicLength(const FloatArray &normalToCrackPlane) override
    { return this->computeLength(); }

    double computeVolumeAround(GaussPoint *gp) override;

    int testElementExtension(ElementExtension ext) override { return 0; }

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QTruss1d_Name; }
    const char *giveClassName() const override { return "QTruss1d"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_line_2;}


    MaterialMode giveMaterialMode() override { return _1dMat; }
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeGaussPoints() override;
};
} // end namespace oofem
#endif // truss1d_h
