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

#include "timesteptestmaterial.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "convergenceexception.h"

namespace oofem {
REGISTER_Material(TimeStepTestMaterial);

TimeStepTestMaterial :: TimeStepTestMaterial(int n, Domain *d) : IsotropicLinearElasticMaterial(n,d)
{ }



void
TimeStepTestMaterial :: initializeFrom(InputRecord &ir)
{
    IsotropicLinearElasticMaterial:: initializeFrom(ir);
    IR_GIVE_FIELD(ir, reduceTimeStepTime, _IFT_TimeStepTestMaterial_rtst);
    IR_GIVE_FIELD(ir, reduceTimeStepFactor, _IFT_TimeStepTestMaterial_rtsf);
   
}


FloatArrayF<6>
TimeStepTestMaterial :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{

    auto Sigma = LinearElasticMaterial :: giveRealStressVector_3d(strain, gp, tStep);
    if(reduceTimeStepTime <= tStep->giveTargetTime() ) {
      if(reduceTimeStepFactor < 0) {
	throw ConvergenceException("Testing time step reduction in case of convergence Error");
      } else {
	tStep->setTimeStepReductionFactor(reduceTimeStepFactor);
      }
    }
    return Sigma;
}



} // end namespace oofem
