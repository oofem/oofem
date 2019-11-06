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

#ifndef twophasemat_h
#define twophasemat_h

#include "tm/Materials/transportmaterial.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "scalarfunction.h"

///@name Input fields for IsotropicHeatTransferMaterial
//@{
#define _IFT_TwoPhaseMaterial_Name "twophasemat"
#define _IFT_TwoPhaseMaterial_mat "mat" /// Array of individual materials (void material, reference material)
//@}

namespace oofem {

/**
 * This class implements an isotropic linear heat  material. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 */
class TwoPhaseMaterial : public TransportMaterial
{
protected:
    IntArray slaveMaterial; 
 
public:
    TwoPhaseMaterial(int n, Domain *d) : TransportMaterial(n, d) { }

    void initializeFrom(InputRecord &ir) override;

    FloatArrayF<3> computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

     double giveCharacteristicValue(MatResponseMode mode,
                                   GaussPoint *gp,
                                   TimeStep *tStep) const override;

    //virtual double giveMaturityT0() const { return maturityT0; }
protected:
    TransportMaterial *giveMaterial(int i) const;
    double giveVof (GaussPoint* gp, TimeStep* tStep) const;
    const char *giveClassName() const override { return "TwoPhaseMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_TwoPhaseMaterial_mat; }

};

} // end namespace oofem
#endif // twophasemat_h
