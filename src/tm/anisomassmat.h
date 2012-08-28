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

#ifndef anisomassmat_h
#define anisomassmat_h

#include "transportmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
class GaussPoint;

/**
 * Material status class for a linear anisotropic mass transfer material. The material is used with the Tr1Darcy element.
 *
 * @author Carl Sandström
 */
class AnisotropicMassTransferMaterialStatus : public TransportMaterialStatus
{
private:
    FloatArray pressureGradient;     ///< Vector containing the last used pressure gradient
    FloatArray seepageVelocity;      ///< Vector containing the last computed velocity

public:
    AnisotropicMassTransferMaterialStatus(int n, Domain *d, GaussPoint *g);

    /// Set pressureGradient
    void setPressureGradient(const FloatArray &gradP);

    /// Set seepageVelocity
    void setSeepageValocity(const FloatArray &w);

    /// Return last pressure gradient
    const FloatArray &giveGradP() { return pressureGradient; };

    /// Returns last seepage velocity vector
    const FloatArray &giveSeepageVelocity() { return seepageVelocity; };

    virtual const char *giveClassName() const { return "AnisotropicMassTransferMaterialStatus"; }
    virtual classType giveClassID() const { return AnisotropicMassTransferMaterialStatusClass; }
};

/**
 *
 * Class for an anisotropic linear transport material. The constitutive equation is given as
 * @f[ \mathbf{w}=-\mathbf{K} \mathbf{\nabla} p @f]
 * where @f$ \mathbm{w} @f$ is the seepage velocity, @f$ \mathbm{K} @f$ is the permeability which is given in the input file
 * and @f$ p @f$ is the pressure.
 *
 * @author Carl Sandström
 *
 */
class AnisotropicMassTransferMaterial : public TransportMaterial
{
protected:
    FloatMatrix k; ///< Conductivity/permeability matrix. This matrix is read from the input file and should be given row-wise as a vector of 4, eg "C 4 1 0 0 1".

public:

    AnisotropicMassTransferMaterial(int n, Domain *d) : TransportMaterial(n, d) { };
    virtual ~AnisotropicMassTransferMaterial() { };

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);

    virtual const char *giveClassName() const { return "AnisotropicMassTransferMaterial"; }
    virtual classType giveClassID() const { return AnisotropicMassTransferMaterialClass; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new AnisotropicMassTransferMaterialStatus(1, domain, gp);  }
};
} // end namespace oofem
#endif // anisomassmat_h
