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

#include "domain.h"
#include "nonlocalbarrier.h"
#include "graddamagematerialextensioninterface.h"
#include "inputrecord.h"
#include "floatmatrix.h"

#include <list>


namespace oofem {
// constructor
GradientDamageMaterialExtensionInterface :: GradientDamageMaterialExtensionInterface(Domain *d)  : 
    Interface(),
    dom(d)
{
}

void
GradientDamageMaterialExtensionInterface :: giveGradientDamageStiffnessMatrix_dd_NN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.clear();
}

void
GradientDamageMaterialExtensionInterface :: giveGradientDamageStiffnessMatrix_dd_BN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.clear();
}
  

void
GradientDamageMaterialExtensionInterface :: initializeFrom(InputRecord &ir)
{
    // read the characteristic length
    IR_GIVE_FIELD(ir, internalLength, _IFT_GradientDamageMaterialExtensionInterface_l);
}

GradientDamageMaterialStatusExtensionInterface :: GradientDamageMaterialStatusExtensionInterface() : Interface()
{
}


void
GradientDamageMaterialStatusExtensionInterface :: initTempStatus()
{
    tempLocalDamageDrivingVariable = localDamageDrivingVariable;
    tempNonlocalDamageDrivingVariable = nonlocalDamageDrivingVariable;
    tempNonlocalDamageDrivingVariableGrad = nonlocalDamageDrivingVariableGrad;
}


void
GradientDamageMaterialStatusExtensionInterface :: updateYourself(TimeStep *tStep)
{
    localDamageDrivingVariable = tempLocalDamageDrivingVariable;
    nonlocalDamageDrivingVariable = tempNonlocalDamageDrivingVariable;
    nonlocalDamageDrivingVariableGrad = tempNonlocalDamageDrivingVariableGrad;
}
} // end namespace oofem
