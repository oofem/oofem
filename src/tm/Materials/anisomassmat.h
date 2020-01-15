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

#ifndef anisomassmat_h
#define anisomassmat_h

#include "tm/Materials/transportmaterial.h"
#include "floatmatrixf.h"

///@name Input fields for AnisotropicMassTransferMaterial
//@{
#define _IFT_AnisotropicMassTransferMaterial_Name "anisomass"
#define _IFT_AnisotropicMassTransferMaterial_c "c"
//@}

namespace oofem {
/**
 *
 * Class for an anisotropic linear transport material. The constitutive equation is given as
 * @f[ \mathbf{w}=-\mathbf{K} \mathbf{\nabla} p @f]
 * where @f$ \mathbf{w} @f$ is the seepage velocity, @f$ \mathbf{K} @f$ is the permeability which is given in the input file
 * and @f$ p @f$ is the pressure.
 *
 * @author Carl Sandström
 *
 */
class AnisotropicMassTransferMaterial : public TransportMaterial
{
protected:
    FloatMatrixF<3,3> k; ///< Conductivity/permeability matrix.

public:
    AnisotropicMassTransferMaterial(int n, Domain * d) : TransportMaterial(n, d) { }

    void initializeFrom(InputRecord &ir) override;

    FloatArrayF<3> computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    double giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveInputRecordName() const override { return _IFT_AnisotropicMassTransferMaterial_Name; }
    const char *giveClassName() const override { return "AnisotropicMassTransferMaterial"; }
};
} // end namespace oofem
#endif // anisomassmat_h
