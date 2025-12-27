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

#ifndef libeam2d_h
#define libeam2d_h

#include "sm/Elements/structuralelement.h"
#include "sm/CrossSections/layeredcrosssection.h"

#define _IFT_LIBeam2d_Name "libeam2d"
#define _IFT_LIBeam2d_XZ "xz"
#define _IFT_LIBeam2d_XY "xy"

namespace oofem {
class FEI2dLineLin;
class ParamKey;

/**
 * A 2-dimensional Linear Isoparametric
 * Mindlin theory beam element, with reduced integration.
 */
class LIBeam2d : public StructuralElement, public LayeredCrossSectionInterface
{
protected:
    /// Interpolation
    static FEI2dLineLin interpolationXZ;
    static FEI2dLineLin interpolationXY;

    static ParamKey IPK_LIBeam2d_XZ;
    static ParamKey IPK_LIBeam2d_XY;

public:
    double pitch, length;
    bool xz = true;
    bool xy = false; //by default element in xz plane

    LIBeam2d(int n, Domain * aDomain);
    virtual ~LIBeam2d() { }

    FEInterpolation *giveInterpolation() const override;

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override { computeLumpedMassMatrix(answer, tStep); }
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;

    // layered cross section support functions
    void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                    GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    Interface *giveInterface(InterfaceType it) override;

    int computeNumberOfDofs() override { return 6; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    double computeVolumeAround(GaussPoint *gp) override;
    Element_Geometry_Type giveGeometryType() const override {return EGT_line_1;}


    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LIBeam2d_Name; }
    const char *giveClassName() const override { return "LIBeam2d"; }
    void initializeFrom(InputRecord &ir, int priority) override;
    void initializeFinish() override;

    MaterialMode giveMaterialMode() override { return _2dBeam; }

protected:
    // edge load support
    void giveEdgeDofMapping(IntArray &answer, int) const override;
    double computeEdgeVolumeAround(GaussPoint *, int) override;
    int computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *) override;
    int computeLoadGToLRotationMtrx(FloatMatrix &answer) override;
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode) override;
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &) override;
    void computeGaussPoints() override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;

    double computeLength() override;
    double givePitch();
};
} // end namespace oofem
#endif // libeam2d_h
