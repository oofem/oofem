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

#ifndef nonlinearheatmat_h
#define nonlinearheatmat_h

#include "transportmaterial.h"

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
    double C;
    /// Indicates the level of nonlinearity in the model
    double alpha;

public:

    NonlinearMassTransferMaterial(int n, Domain * d) : TransportMaterial(n, d) { }
    virtual ~NonlinearMassTransferMaterial() { }

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);

    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);

    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual const char *giveInputRecordName() const { return _IFT_NonlinearMassTransferMaterial_Name; }
    virtual const char *giveClassName() const { return "NonlinearMassTransferMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);
};
} // end namespace oofem
#endif // nonlinearheatmat_h
