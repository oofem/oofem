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

#ifndef tr1_2d_supg_axi_h
#define tr1_2d_supg_axi_h

#include "tr1_2d_supg.h"

#define _IFT_TR1_2D_SUPG_AXI_Name "tr1supgaxi"

namespace oofem {
/**
 * Class representing 2d linear axisymmetric triangular element
 * for solving incompressible fluid with SUPG solver
 *
 * This class is similar to TR1_2D_SUPG2_AXI, but difference is in handling
 * multiple fluids. This class uses rule of mixture which interpolates the
 * properties using VOF value, requiring the use of twofluidmaterial class
 * as material model for this situation.
 */
class TR1_2D_SUPG_AXI : public TR1_2D_SUPG
{
protected:
    /// Radius at element center.
    double rc;

public:
    TR1_2D_SUPG_AXI(int n, Domain * d);

    void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *tStep) override;
    void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *tStep) override;
    void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep) override;
    void computePressureTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *tStep) override;
    void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *tStep) override;
    void computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computePressureTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;

    void updateStabilizationCoeffs(TimeStep *tStep) override;
    void computeBCRhsTerm_MB(FloatArray &answer, TimeStep *tStep) override;
    void computeBCRhsTerm_MC(FloatArray &answer, TimeStep *tStep) override;
    void computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep) override;

    void computeSlipWithFrictionBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *tStep) override;
    void computePenetrationWithResistanceBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *tStep) override;
    void computeOutFlowBCTerm_MB(FloatMatrix &answer, int side, TimeStep *tStep) override;
    void computeBCLhsTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeBCLhsPressureTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;

    double computeVolumeAround(GaussPoint *gp) override;
    void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi) override;

    // definition
    const char *giveClassName() const override { return "TR1_2D_SUPG_AXI"; }
    const char *giveInputRecordName() const override { return _IFT_TR1_2D_SUPG_AXI_Name; }
    MaterialMode giveMaterialMode() override { return _2dAxiFlow; }

protected:
    void computeGaussPoints() override;
    void computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;
    void computeDeviatoricStress(FloatArray &answer, const FloatArray &eps, GaussPoint *gp, TimeStep *tStep) override;
    void computeTangent(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void initGeometry() override;
    virtual double computeRadiusAt(GaussPoint *gp);
    virtual void computeBMtrx(FloatMatrix &answer, GaussPoint *gp);
    void computeNVector(FloatArray &answer, GaussPoint *gp);
};
} // end namespace oofem
#endif // tr1_2d_supg_axi_h
