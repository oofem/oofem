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

#ifndef hyperelasticmaterial_h
#define hyperelasticmaterial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"

///@name Input fields for HyperElasticMaterial
//@{
#define _IFT_HyperElasticMaterial_Name "hyperelmat"
#define _IFT_HyperElasticMaterial_k "k"
#define _IFT_HyperElasticMaterial_g "g"
//@}

namespace oofem {
/**
 * Saint Venant–Kirchhoff model defined by shear and bulk modulus.
 * @todo Should we even have this? Isn't this identical to the isotropic linear elastic model? / Mikael
 */
class HyperElasticMaterial : public StructuralMaterial
{
protected:
    double K; ///< Bulk modulus.
    double G; ///< Shear modulus.

public:
    HyperElasticMaterial(int n, Domain * d);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep);


    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &reducedStrain, TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual const char *giveInputRecordName() const { return _IFT_HyperElasticMaterial_Name; }
    virtual const char *giveClassName() const { return "HyperElasticMaterial"; }
};
} // end namespace oofem
#endif
