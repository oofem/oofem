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

#ifndef cohsur3d_h
#define cohsur3d_h

#include "sm/Elements/structuralelement.h"

///@name Input fields for CohSur3d
//@{
#define _IFT_CohesiveSurface3d_Name "cohsur3d"
#define _IFT_CohSur3d_kx "kx"
#define _IFT_CohSur3d_ky "ky"
#define _IFT_CohSur3d_kz "kz"
#define _IFT_CohSur3d_area "area"
//@}

namespace oofem {
/**
 * This class implements a cohesive surface element used by the
 * cohesive particle model.
 *
 * @author Milan Jirasek
 */
class CohesiveSurface3d : public StructuralElement
{
protected:
    double area, length;
    FloatArray center; ///< Coordinates of the center of the cohesive surface.
    FloatMatrix lcs; ///< Matrix defining the local coordinate system.

    ///@name Shift constants of periodic particles (near boundary of periodic cell).
    //@{
    int kx, ky, kz;
    double kxa, kyb, kzc;
    //@}

public:
    CohesiveSurface3d(int n, Domain * d);
    virtual ~CohesiveSurface3d() { }

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui) override;
    double computeVolumeAround(GaussPoint *gp) override;
    int computeNumberOfDofs() override { return 6 * giveNumberOfNodes(); }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    double giveLength();
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override { }
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    // definition & identification
    const char *giveClassName() const override { return "CohesiveSurface3d"; }
    const char *giveInputRecordName() const override { return _IFT_CohesiveSurface3d_Name; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_unknown;}


    // input and output
    void initializeFrom(InputRecord &ir) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

protected:
    void computeGaussPoints() override;
    void evaluateCenter();
    void evaluateLocalCoordinateSystem();

    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    MaterialMode giveMaterialMode() override { return _3dInterface; }
};
} // namespace oofem
#endif
