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

#ifndef isoheatmat_h
#define isoheatmat_h

#include "transportmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
/**
 * This class implements a isotropic linear heat  material. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 */
class IsotropicHeatTransferMaterial : public TransportMaterial
{
protected:
    double conductivity; ///< Conductivity (k in input file).
    double capacity;     ///< Capacity (c in input file).

public:
    IsotropicHeatTransferMaterial(int n, Domain *d) : TransportMaterial(n, d) { }
    virtual ~IsotropicHeatTransferMaterial() { }

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime);

    virtual double giveIsotropicConductivity(GaussPoint *gp) {return conductivity;};
                                          
    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);


    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual const char *giveClassName() const { return "IsotropicHeatTransferMaterial"; }
    virtual classType giveClassID() const { return IsotropicHeatTransferMaterialClass; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double give(int aProperty, GaussPoint *gp);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TransportMaterialStatus(1, domain, gp);  }
};
} // end namespace oofem
#endif // isoheatmat_h
