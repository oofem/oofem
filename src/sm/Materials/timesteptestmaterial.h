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

#ifndef timesteptestmaterial_h
#define timesteptestmaterial_h

#include "isolinearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"


///@name Input fields for IsotropicLinearElasticMaterial
//@{
#define _IFT_TimeStepTestMaterial_Name "tstmat"
#define _IFT_TimeStepTestMaterial_rtst "rtst"
#define _IFT_TimeStepTestMaterial_rtsf "rtsf"
//@}

namespace oofem {
class GaussPoint;

  class TimeStepTestMaterial: public IsotropicLinearElasticMaterial
  {
protected:

    double reduceTimeStepTime;
    double reduceTimeStepFactor;

public:

    TimeStepTestMaterial(int n, Domain *d);

    const char *giveClassName() const override { return "TimeStepTestMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_TimeStepTestMaterial_Name; }

    void initializeFrom(InputRecord &ir) override;
    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;


};
} // end namespace oofem
#endif // isolinearelasticmaterial_h
