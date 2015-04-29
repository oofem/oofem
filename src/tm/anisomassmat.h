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

#include "transportmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

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
    FloatMatrix k; ///< Conductivity/permeability matrix. This matrix is read from the input file and should be given row-wise as a vector of 4, eg "C 4 1 0 0 1".

public:
    AnisotropicMassTransferMaterial(int n, Domain * d) : TransportMaterial(n, d) { }
    virtual ~AnisotropicMassTransferMaterial() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep);

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual double giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual const char *giveInputRecordName() const { return _IFT_AnisotropicMassTransferMaterial_Name; }
    virtual const char *giveClassName() const { return "AnisotropicMassTransferMaterial"; }
};
} // end namespace oofem
#endif // anisomassmat_h
