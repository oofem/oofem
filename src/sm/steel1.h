/* $Header: /home/cvs/bp/oofem/sm/src/steel1.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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

//   *********************************************************************
//   *** CLASS STEEL1 - Lineary elastic - perfectly plastic Hook  Material
//   *********************************************************************

#ifndef steel1_h
#define steel1_h

#include "perfectlyplasticmaterial.h"
#include "isolinearelasticmaterial.h"

namespace oofem {
class Domain;

class Steel1 : public PerfectlyPlasticMaterial
{
    /*
     * This class implements a isotropic perfectly plastic linear material in a finite
     * element problem. A material
     * is an attribute of a domain. It is usually also attribute of many elements.
     *
     * DESCRIPTION
     * ISOTROPIC PERFECTLY PLASTIC Material with HMH plastic condition
     *
     * TASK
     * - Returning standard material stiffness marix for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint.
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     */

protected:
public:

    Steel1(int n, Domain *d);
    ~Steel1() { }

    IRResultType initializeFrom(InputRecord *ir);
    const char *giveClassName() const { return "Steel1MaterialClass"; }
    classType giveClassID()         const { return Steel1MaterialClass; }
    void  updateIfFailure(GaussPoint *, FloatArray *, FloatArray *) { }
protected:

    //
    // yield(YC-like functions) and loading(LC-like functions) criteria specific section
    //

    virtual double      computeYCValueAt(GaussPoint *, FloatArray *, FloatArray *);
    virtual FloatArray *GiveYCStressGradient(GaussPoint *, FloatArray *, FloatArray *);
    virtual FloatArray *GiveLCStressGradient(GaussPoint *, FloatArray *, FloatArray *);
    virtual FloatArray *GiveYCPlasticStrainGradient(GaussPoint *, FloatArray *, FloatArray *);
    virtual FloatArray *GiveLCPlasticStrainGradient(GaussPoint *, FloatArray *, FloatArray *);
    virtual void        updateTempYC(GaussPoint *, FloatArray *, FloatArray *) { }
    virtual void        updateTempLC(GaussPoint *, FloatArray *, FloatArray *) { }
    // update during computation

    // auxiliary function
    double      computeJ2InvariantAt(FloatArray *);
};
} // end namespace oofem
#endif // steel1_h
