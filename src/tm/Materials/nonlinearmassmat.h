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

#ifndef nonlinearheatmat_h
#define nonlinearheatmat_h

#include "tm/Materials/transportmaterial.h"

///@name Input fields for NonlinearMassTransferMaterial
//@{
#define _IFT_NonlinearMassTransferMaterial_Name "nonlinmass"
#define _IFT_NonlinearMassTransferMaterial_c "c"
#define _IFT_NonlinearMassTransferMaterial_alpha "alpha"
//@}

namespace oofem {
/**
 *
 * Class for a nonlinear fictitious transport material. The constitutive equation is given as
 * @f[
 * \mathbf{w}=-\left( 1+C \mid \mid \mathbf{\nabla} p \mid\mid^{\alpha}\right) \mathbf{\nabla} p
 * @f]
 * where @f$ \boldsymbol{w} @f$ is the seepage velocity, @f$ \alpha @f$ and @f$ C @f$ are constants and @f$ p @f$ is the pressure.
 *
 * @author Carl Sandström
 *
 */
class NonlinearMassTransferMaterial : public TransportMaterial
{
protected:
    /// Indicates the level of nonlinearity in the model
    double C = 0.;
    /// Indicates the level of nonlinearity in the model
    double alpha = 0.;

public:
    NonlinearMassTransferMaterial(int n, Domain * d) : TransportMaterial(n, d) { }

    void initializeFrom(InputRecord &ir) override;

    double giveCharacteristicValue(MatResponseMode mode,
                                   GaussPoint *gp,
                                   TimeStep *tStep) const override;

    FloatArrayF<3> computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    const char *giveInputRecordName() const override { return _IFT_NonlinearMassTransferMaterial_Name; }
    const char *giveClassName() const override { return "NonlinearMassTransferMaterial"; }
};
} // end namespace oofem
#endif // nonlinearheatmat_h
