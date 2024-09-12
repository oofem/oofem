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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef interfaceelem2dlin_h
#define interfaceelem2dlin_h

#include "sm/Elements/structuralelement.h"

#define _IFT_InterfaceElem2dLin_Name "interface2dlin"
#define _IFT_InterfaceElem2dLin_axisymmode "axisymmode"

namespace oofem {
class FEI2dLineLin;

/**
 * This class implements a two dimensional interface element.
 * The approximatuion of geometry and unknowns is linear.
 */
class InterfaceElem2dLin : public StructuralElement
{
protected:
    static FEI2dLineLin interp;
    /// Flag controlling axisymmetric mode (integration over unit circumferential angle)
    bool axisymmode;

public:
    InterfaceElem2dLin(int n, Domain * d);
    virtual ~InterfaceElem2dLin() { }

    FEInterpolation *giveInterpolation() const override;

    int computeNumberOfDofs() override { return 8; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    Element_Geometry_Type giveGeometryType() const override {return EGT_line_1;}


    double computeVolumeAround(GaussPoint *gp) override;


    int testElementExtension(ElementExtension ext) override { return 0; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_InterfaceElem2dLin_Name; }
    const char *giveClassName() const override { return "InterfaceElem2dLin"; }
    void initializeFrom(InputRecord &ir) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    MaterialMode giveMaterialMode() override { return _2dInterface; }

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override { }
    void computeGaussPoints() override;

    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
};
} // end namespace oofem
#endif // interfaceelem2dlin_h
