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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef tr1_2d_supg_axi_h
#define tr1_2d_supg_axi_h

#include "tr1_2d_supg.h"

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
    TR1_2D_SUPG_AXI(int n, Domain *d);
    virtual ~TR1_2D_SUPG_AXI();

    virtual void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime);
    virtual void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime);
    virtual void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime);
    virtual void computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime);
    virtual void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *atTime);
    virtual void computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime);

    virtual void updateStabilizationCoeffs(TimeStep *tStep);
    virtual void computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime);
    virtual void computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime);

    virtual void computeSlipWithFrictionBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *atTime);
    virtual void computePenetrationWithResistanceBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *atTime);
    virtual void computeOutFlowBCTerm_MB(FloatMatrix &answer, int side, TimeStep *atTime);
    virtual void computeBCLhsTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeBCLhsPressureTerm_MB(FloatMatrix &answer, TimeStep *atTime);

    virtual double computeVolumeAround(GaussPoint *gp);
    virtual void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi);

    // definition
    virtual const char *giveClassName() const { return "TR1_2D_SUPG_AXI"; }
    virtual classType giveClassID() const { return TR1_2D_SUPG_AXIClass; }
    virtual MaterialMode giveMaterialMode() { return _2dAxiFlow; }

protected:
    virtual void computeGaussPoints();
    virtual void computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    virtual void initGeometry();
    virtual double computeRadiusAt(GaussPoint *gp);
    virtual void computeBMtrx(FloatMatrix &answer, GaussPoint *gp);
    void computeNVector(FloatArray &answer, GaussPoint *gp);
};
} // end namespace oofem
#endif // tr1_2d_supg_axi_h
