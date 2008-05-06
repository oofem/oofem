/* $Header: /home/cvs/bp/oofem/tm/src/isoheatmat.h,v 1.1 2003/04/14 16:01:39 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


//   ************************************************
//   *** CLASS ISOTROPIC MATERIAL FOR HEAT
//   ************************************************

#ifndef isoheatmat_h
#define isoheatmat_h

#include "transportmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"


class GaussPoint;


class IsotropicHeatTransferMaterial : public TransportMaterial
{
    /*
     * This class implements a isotropic linear heat  material in a finite
     * element problem. A material
     * is an attribute of a domain. It is usually also attribute of many elements.
     *
     * DESCRIPTION
     * ISOTROPIC Linear Heat Material
     *
     * TASK
     *
     */

protected:

    double k, c;

public:

    IsotropicHeatTransferMaterial(int n, Domain *d) : TransportMaterial(n, d) { }
    ~IsotropicHeatTransferMaterial() { }

    void  giveCharacteristicMatrix(FloatMatrix &answer,
                                   MatResponseForm form,
                                   MatResponseMode mode,
                                   GaussPoint *gp,
                                   TimeStep *atTime); // identification and auxiliary functions

    virtual double  giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime);


    const char *giveClassName() const { return "IsotropicHeatTransferMaterial"; }
    classType giveClassID()         const { return IsotropicHeatTransferMaterialClass; }

    IRResultType initializeFrom(InputRecord *ir);

    // non-standart - returns time independent material constant
    double   give(int);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return NULL; }

protected:
};


#endif // isoheatmat_h
