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

#include "dummycrosssection.h"
#include "dynamicinputrecord.h"
#include "material.h"
#include "gaussintegrationrule.h"
#include "classfactory.h"


namespace oofem {
REGISTER_CrossSection(DummyCrossSection);

DummyCrossSection :: DummyCrossSection(int n, Domain* d) : CrossSection(n, d)
{
}

void
DummyCrossSection :: initializeFrom(InputRecord &ir)
{
    CrossSection :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->matNumber, _IFT_DummyCrossSection_material);
}


void
DummyCrossSection :: giveInputRecord(DynamicInputRecord &input)
{
    CrossSection :: giveInputRecord(input);

    input.setField(this->matNumber, _IFT_DummyCrossSection_material);
}


bool
DummyCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode) const
{
    return this->domain->giveMaterial(this->matNumber)->isCharacteristicMtrxSymmetric(rMode);
}

Material *
DummyCrossSection :: giveMaterial(IntegrationPoint *ip) const
{
    return this->domain->giveMaterial(this->matNumber);
}
int
DummyCrossSection :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->domain->giveMaterial(this->matNumber)->packUnknowns(buff, tStep, gp);
}

int
DummyCrossSection :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->domain->giveMaterial(this->matNumber)->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
DummyCrossSection :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->domain->giveMaterial(this->matNumber)->estimatePackSize(buff, gp);
}
} // end namespace oofem
