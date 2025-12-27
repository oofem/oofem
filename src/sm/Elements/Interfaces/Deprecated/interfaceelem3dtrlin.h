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

#ifndef interfaceelem3dtrlin_h
#define interfaceelem3dtrlin_h

#include "sm/Elements/structuralelement.h"

#define _IFT_InterfaceElement3dTrLin_Name "interface3dtrlin"

namespace oofem {
class FEI2dTrLin;

/**
 * This class implements 3d triangular surface interface element with linear interpolation.
 * @author Jim Brouzoulis
 * @author Borek Patzak
 */
class InterfaceElement3dTrLin : public StructuralElement
{
protected:
    ///@todo Implement FEI3dTrLin, then remove giveIntegrationDomain and giveElementGeometry and overload giveInterpolation()
    static FEI2dTrLin interpolation;

public:
    InterfaceElement3dTrLin(int n, Domain * d);
    virtual ~InterfaceElement3dTrLin() { }

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;
    bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) override;

    int computeNumberOfDofs() override { return 18; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    double computeVolumeAround(GaussPoint *gp) override;

    int testElementExtension(ElementExtension ext) override { return 0; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_InterfaceElement3dTrLin_Name; }
    const char *giveClassName() const override { return "InterfaceElement3dTrLin"; }
    Element_Geometry_Type giveGeometryType() const override { return EGT_wedge_1; }
    integrationDomain giveIntegrationDomain() const override { return _Triangle; }

    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    MaterialMode giveMaterialMode() override { return _3dInterface; }

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override { }
    void computeGaussPoints() override;

    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    void computeLCS(FloatMatrix &answer);
};
} // end namespace oofem
#endif // interfaceelem3dtrlin_h
