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

#include "simpletransportcrosssection.h"
#include "transportmaterial.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

namespace oofem {
REGISTER_CrossSection(SimpleTransportCrossSection);

SimpleTransportCrossSection :: SimpleTransportCrossSection(int n, Domain *d) : TransportCrossSection(n, d) { }

SimpleTransportCrossSection :: ~SimpleTransportCrossSection() { }


IRResultType
SimpleTransportCrossSection :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->matNumber, _IFT_SimpleTransportCrossSection_material);
    this->propertyDictionary.clear();
    if ( ir->hasField(_IFT_SimpleTransportCrossSection_thickness) ) {
        double thickness;
        IR_GIVE_FIELD(ir, thickness, _IFT_SimpleTransportCrossSection_thickness);
        this->propertyDictionary.add(CS_Thickness, thickness);
    }

    return TransportCrossSection :: initializeFrom(ir);
}


void
SimpleTransportCrossSection :: giveInputRecord(DynamicInputRecord &input)
{
    TransportCrossSection :: giveInputRecord(input);

    input.setField(this->matNumber, _IFT_SimpleTransportCrossSection_material);
}


int
SimpleTransportCrossSection :: checkConsistency()
{
    Material *mat = dynamic_cast< TransportMaterial * >( this->domain->giveMaterial(this->matNumber) );
    if ( !mat ) {
        return 0;
    }
    return TransportCrossSection :: checkConsistency();
}


TransportMaterial *
SimpleTransportCrossSection :: giveMaterial()
{
    return dynamic_cast< TransportMaterial * >( this->domain->giveMaterial(this->matNumber) );
}


int
SimpleTransportCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep)
{
    return this->domain->giveMaterial(this->matNumber)->giveIPValue(answer, ip, type, tStep);
}


bool
SimpleTransportCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode)
{
    return this->domain->giveMaterial(this->matNumber)->isCharacteristicMtrxSymmetric(rMode);
}


int
SimpleTransportCrossSection :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->domain->giveMaterial(this->matNumber)->packUnknowns(buff, tStep, gp);
}

int
SimpleTransportCrossSection :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->domain->giveMaterial(this->matNumber)->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
SimpleTransportCrossSection :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->domain->giveMaterial(this->matNumber)->estimatePackSize(buff, gp);
}

} // end namespace oofem
