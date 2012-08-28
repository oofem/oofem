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

#ifndef nonlinearheatmat_h
#define nonlinearheatmat_h

#include "transportmaterial.h"
#include "anisomassmat.h"

namespace oofem {

/**
 *
 * Class for a nonlinear fictitious transport material. The constitutive equation is given as
 *
 * @f[ \mathbf{w}=-\left( 1+C \mid \mid \mathbf{\nabla} p \mid\mid^{\alpha}\right) \mathbf{\nabla} p @f]
 *
 * where @f$ \mathbm{w} @f$ is the seepage velocity, @f$ \alpha @f$ and @f$ C @f$ are constants and @f$ p @f$ is the pressure.
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

    NonlinearMassTransferMaterial(int n, Domain *d) : TransportMaterial(n, d) { };
    virtual ~NonlinearMassTransferMaterial() { };

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);

    virtual const char *giveClassName() const { return "NonlinearMassTransferMaterial"; };
    virtual classType giveClassID() const { return NonlinearMassTransferMaterialClass; };

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new AnisotropicMassTransferMaterialStatus(1, domain, gp);  };
};
} // end namespace oofem
#endif // nonlinearheatmat_h
