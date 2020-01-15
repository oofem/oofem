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

#ifndef tet1_3d_supg_h
#define tet1_3d_supg_h

#include "supgelement2.h"
#include "levelsetpcs.h"

#define _IFT_Tet1_3D_SUPG_Name "tet1supg"

namespace oofem {
class FEI3dTetLin;

/**
 * Class representing 3d linear tetrahedral element
 * for solving incompressible fluid with SUPG solver
 */
class Tet1_3D_SUPG : public SUPGElement2, public LevelSetPCSElementInterface
{
protected:
    static FEI3dTetLin interpolation;

public:
    Tet1_3D_SUPG(int n, Domain * d);

    // definition
    const char *giveClassName() const override { return "Tet1_3D_SUPG"; }
    const char *giveInputRecordName() const override { return _IFT_Tet1_3D_SUPG_Name; }
    MaterialMode giveMaterialMode() override { return _3dFlow; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    int computeNumberOfDofs() override;

    Interface *giveInterface(InterfaceType t) override;

    double computeCriticalTimeStep(TimeStep *tStep) override;
    double computeVolumeAround(GaussPoint *gp) override;

    void updateStabilizationCoeffs(TimeStep *tStep) override;

    double LS_PCS_computeF(LevelSetPCS *, TimeStep *tStep) override;
    void LS_PCS_computedN(FloatMatrix &answer) override;
    double LS_PCS_computeVolume() override;
    double LS_PCS_computeS(LevelSetPCS *, TimeStep *tStep) override;
    void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

protected:
    void computeGaussPoints() override;
    void computeDeviatoricStress(FloatArray &answer, const FloatArray &eps, GaussPoint *gp, TimeStep *tStep) override;
    void computeTangent(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void computeNuMatrix(FloatMatrix &answer, GaussPoint *gp) override;
    void computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) override;
    void computeBMatrix(FloatMatrix &anwer, GaussPoint *gp) override;
    void computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp) override;
    void computeNpMatrix(FloatMatrix &answer, GaussPoint *gp) override;
    void computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp) override;
    void computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) override;
    void computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) override;
    int  giveNumberOfSpatialDimensions() override;
};
} // end namespace oofem
#endif // tet1_3d_supg_h
