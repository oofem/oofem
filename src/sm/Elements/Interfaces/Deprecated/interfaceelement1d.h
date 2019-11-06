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

#ifndef interfaceelement1d_h
#define interfaceelement1d_h

#include "sm/Elements/structuralelement.h"

///@name Input fields for Material
//@{
#define _IFT_InterfaceElem1d_Name "interface1d"
#define _IFT_InterfaceElem1d_refnode "refnode"
#define _IFT_InterfaceElem1d_normal "normal"
#define _IFT_InterfaceElem1d_dofIDs "dofids"
//@}

namespace oofem {
/**
 * This class implements a one-dimensional interface element connecting two nodes (with the same position)
 * In order to compute normal and tangential direction of the slip plane, a reference node or specific direction is needed.
 */
class InterfaceElem1d : public StructuralElement
{
protected:
    enum cmode { ie1d_1d, ie1d_2d, ie1d_3d } mode;
    int referenceNode;
    FloatArray normal;
    IntArray dofids;

public:
    InterfaceElem1d(int n, Domain * d);
    virtual ~InterfaceElem1d() { }

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override  { computeLumpedMassMatrix(answer, tStep); }

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;
    bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) override;

    int computeNumberOfDofs() override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    double computeVolumeAround(GaussPoint *gp) override;

    int testElementExtension(ElementExtension ext) override { return 0; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_InterfaceElem1d_Name; }
    const char *giveClassName() const override { return "InterfaceElem1d"; }
    void initializeFrom(InputRecord &ir) override;
    Element_Geometry_Type giveGeometryType() const override { return EGT_point; }
    integrationDomain giveIntegrationDomain() const override { return _Point; }
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    MaterialMode giveMaterialMode() override;

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override { }
    void computeGaussPoints() override;
    void evaluateLocalCoordinateSystem(FloatMatrix &);

    void computeLocalSlipDir(FloatArray &normal);
    cmode giveCoordMode() const { return this->mode; }
    void setCoordMode();
};
} // end namespace oofem
#endif // interfaceelement1d_h
