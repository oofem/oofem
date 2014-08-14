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

#ifndef steel1_h
#define steel1_h

#include "perfectlyplasticmaterial.h"
#include "Materials/isolinearelasticmaterial.h"

///@name Input fields for Steel1
//@{
#define _IFT_Steel1_Name "steel1"
#define _IFT_Steel1_ry "ry"
//@}

namespace oofem {
class Domain;
/**
 * This class implements a isotropic perfectly plastic linear material in a finite
 * element problem.
 */
class Steel1 : public PerfectlyPlasticMaterial
{
public:
    Steel1(int n, Domain * d);
    virtual ~Steel1() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveInputRecordName() const { return _IFT_Steel1_Name; }
    virtual const char *giveClassName() const { return "Steel1MaterialClass"; }
    virtual void updateIfFailure(GaussPoint *gp, FloatArray *, FloatArray *) { }

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
    double computeJ2InvariantAt(FloatArray *);
};
} // end namespace oofem
#endif // steel1_h
