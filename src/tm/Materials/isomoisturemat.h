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

#ifndef isomoisturemat_h
#define isomoisturemat_h

#include "tm/Materials/transportmaterial.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

namespace oofem {
/**
 * This class implements a isotropic moisture transport material. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 */
class IsotropicMoistureTransferMaterial : public TransportMaterial
{
public:
    IsotropicMoistureTransferMaterial(int n, Domain * d) : TransportMaterial(n, d) { }

    void initializeFrom(InputRecord &ir) override;

    FloatArrayF<3> computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    double giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    virtual double givePermeability(GaussPoint *gp, TimeStep *tStep) const = 0;
    virtual double giveMoistureCapacity(GaussPoint *gp, TimeStep *tStep) const = 0;
    //sorption isotherm, return total water mass [kg/m3]
    virtual double giveMoistureContent(double humidity) const { return 0.; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    const char *giveClassName() const override { return "IsotropicMoistureTransferMaterial"; }

};
} // end namespace oofem
#endif // isomoisturemat_h
