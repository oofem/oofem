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

#ifndef qplanstrain_h
#define qplanstrain_h

#include "sm/Elements/structural2delement.h"
#include "zznodalrecoverymodel.h"

#define _IFT_QPlaneStrain_Name "qplanestrain"

namespace oofem {
class FEI2dQuadQuad;

/**
 * This class implements an Quadratic isoparametric eight-node quadrilateral plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
class QPlaneStrain : public PlaneStrainElement, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI2dQuadQuad interpolation;

public:
    QPlaneStrain(int N, Domain * d);
    virtual ~QPlaneStrain() { }

    FEInterpolation *giveInterpolation() const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QPlaneStrain_Name; }
    const char *giveClassName() const override { return "QPlaneStrain"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_quad_2;}


    int testElementExtension(ElementExtension ext) override { return 0; } ///@todo //check this probably ok now when derived from PE-element

    Interface *giveInterface(InterfaceType it) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

};
} // end namespace oofem
#endif // qplanstrain_h
