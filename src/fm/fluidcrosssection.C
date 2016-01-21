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

#include "fluidcrosssection.h"
#include "fluiddynamicmaterial.h"
#include "dynamicinputrecord.h"
#include "domain.h"
#include "classfactory.h"

namespace oofem {
REGISTER_CrossSection(FluidCrossSection);

FluidCrossSection :: FluidCrossSection(int n, Domain *d) : CrossSection(n, d), matNumber(0) { }

FluidCrossSection :: ~FluidCrossSection() { }


IRResultType
FluidCrossSection :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->matNumber, _IFT_FluidCrossSection_material);

    return CrossSection :: initializeFrom(ir);
}


void
FluidCrossSection :: giveInputRecord(DynamicInputRecord &input)
{
    CrossSection :: giveInputRecord(input);

    input.setField(this->matNumber, _IFT_FluidCrossSection_material);
}


int
FluidCrossSection :: checkConsistency()
{
    CrossSection :: checkConsistency();
    return dynamic_cast< FluidDynamicMaterial * >( this->domain->giveMaterial(this->matNumber) ) != NULL;
}


double
FluidCrossSection :: giveDensity(GaussPoint *gp)
{
    return this->domain->giveMaterial(this->matNumber)->give('d', gp);
}


int
FluidCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep)
{
    return this->domain->giveMaterial(this->matNumber)->giveIPValue(answer, ip, type, tStep);
}


bool
FluidCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode)
{
    return this->domain->giveMaterial(this->matNumber)->isCharacteristicMtrxSymmetric(rMode);
}


FluidDynamicMaterial *
FluidCrossSection :: giveFluidMaterial()
{
    return static_cast< FluidDynamicMaterial * >( this->domain->giveMaterial(this->matNumber) );
}


int
FluidCrossSection :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->domain->giveMaterial(this->matNumber)->packUnknowns(buff, tStep, gp);
}

int
FluidCrossSection :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->domain->giveMaterial(this->matNumber)->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
FluidCrossSection :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->domain->giveMaterial(this->matNumber)->estimatePackSize(buff, gp);
}

} // end namespace oofem
