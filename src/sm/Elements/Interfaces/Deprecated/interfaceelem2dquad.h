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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef interfaceelem2dquad_h
#define interfaceelem2dquad_h

#include "sm/Elements/structuralelement.h"

#define _IFT_InterfaceElem2dQuad_Name "interface2dquad"
#define _IFT_InterfaceElem2dQuad_axisymmode "axisymmode"

namespace oofem {
class ParamKey;
class FEI2dLineQuad;

/**
 * This class implements a two dimensional interface element.
 * Even if geometry approx is quadratic, the element is assumed straight
 * If not straight, the rotation matrix depends on actual integration point
 * and stiffness and strain computations should be modified.
 */
class InterfaceElem2dQuad : public StructuralElement
{
protected:
    static FEI2dLineQuad interp;
    /// Flag controlling axisymmetric mode (integration over unit circumferential angle)
    bool axisymmode;

    static ParamKey IPK_InterfaceElem2dQuad_axisymmode;

public:
    InterfaceElem2dQuad(int n, Domain * d);
    virtual ~InterfaceElem2dQuad() { }

    FEInterpolation *giveInterpolation() const override;

    int computeNumberOfDofs() override { return 12; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    double computeVolumeAround(GaussPoint *gp) override;


    int testElementExtension(ElementExtension ext) override { return 0; }

    Interface *giveInterface(InterfaceType) override { return NULL; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_line_2;}


#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_InterfaceElem2dQuad_Name; }
    const char *giveClassName() const override { return "InterfaceElem2dQuad"; }
    void initializeFrom(InputRecord &ir, int priority) override;
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
#endif // interfaceelem2dquad_h
