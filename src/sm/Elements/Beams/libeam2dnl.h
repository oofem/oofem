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

#ifndef libeam2dnl_h
#define libeam2dnl_h

#include "sm/Elements/nlstructuralelement.h"
#include "sm/CrossSections/layeredcrosssection.h"

#define _IFT_LIBeam2dNL_Name "libeam2dNL"

namespace oofem {
/**
 * This class implements a 2-dimensional Linear Isoparametric
 * Mindlin theory beam element, with reduced integration.
 * Geometric nonlinearities are taken into account.
 */
class LIBeam2dNL : public NLStructuralElement, public LayeredCrossSectionInterface
{
protected:
    double pitch, length;

public:
    LIBeam2dNL(int n, Domain *d);
    virtual ~LIBeam2dNL() { }

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { computeLumpedMassMatrix(answer, tStep); }
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override;

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    // layered cross section support functions
    void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                    GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep) override;

    Interface *giveInterface(InterfaceType it) override;

    int computeNumberOfDofs() override { return 6; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    double computeVolumeAround(GaussPoint *gp) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LIBeam2dNL_Name; }
    const char *giveClassName() const override { return "LIBeam2dNL"; }
    void initializeFrom(InputRecord &ir, int priority) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
#endif

    integrationDomain giveIntegrationDomain() const override { return _Line; }
    MaterialMode giveMaterialMode() override { return _2dBeam; }
    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

protected:
    // edge load support
    void giveEdgeDofMapping(IntArray &answer, int) const override;
    double computeEdgeVolumeAround(GaussPoint *, int) override;
    int computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *) override;
    int computeLoadGToLRotationMtrx(FloatMatrix &answer) override;
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode) override;

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    // nonlinear part of geometrical eqs. for i-th component of strain vector.
    void computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int);
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;
    void computeGaussPoints() override;
    double computeLength() override;
    double givePitch();
};
} // end namespace oofem
#endif // libeam2dnl_h
