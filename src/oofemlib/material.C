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

#include "material.h"
#include "verbose.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "contextioerr.h"

namespace oofem {

Material :: Material(int n, Domain* d) : FEMComponent(n, d), propertyDictionary(), castingTime ( -1. ) { }


Material :: ~Material()
{
}


double
Material :: give(int aProperty, GaussPoint *gp)
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
// tStep allows time dependent behavior to be taken into account
{
    double value = 0.0;

    if ( propertyDictionary.includes(aProperty) ) {
        value = propertyDictionary.at(aProperty);
    } else {
        OOFEM_ERROR("property #%d on element %d and GP %d not defined", aProperty, gp->giveElement()->giveNumber(), gp->giveNumber() );
    }

    return value;
}


bool
Material :: hasProperty(int aProperty, GaussPoint *gp)
// Returns true if the aProperty is defined on a material
{
    return propertyDictionary.includes(aProperty);
}


void
Material :: modifyProperty(int aProperty, double value, GaussPoint *gp)
{
    if ( propertyDictionary.includes(aProperty) ) {
        propertyDictionary.at(aProperty) = value;
    } else {
        OOFEM_ERROR("property #%d on element %d and GP %d not defined", aProperty, gp->giveElement()->giveNumber(), gp->giveNumber() );
    }
}


IRResultType
Material :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    double value;

#  ifdef VERBOSE
    // VERBOSE_PRINT1 ("Instanciating material ",this->giveNumber())
#  endif

    IR_GIVE_FIELD(ir, value, _IFT_Material_density);
    propertyDictionary.add('d', value);

    this->castingTime = -1.e10;
    IR_GIVE_OPTIONAL_FIELD(ir, castingTime, _IFT_Material_castingtime);

    return IRRT_OK;
}


void
Material :: giveInputRecord(DynamicInputRecord &input)
{
    FEMComponent :: giveInputRecord(input);
    input.setField(this->propertyDictionary.at('d'), _IFT_Material_density);
    input.setField(this->castingTime, _IFT_Material_castingtime);
}


int
Material :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return 0;
}


int
Material :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_MaterialNumber ) {
        answer.resize(1);
        answer.at(1) = this->giveNumber();
        return 1;
    }

    answer.clear();
    return 0;
}


void
Material :: printYourself()
// Prints the receiver on screen.
{
    printf("Material with properties : \n");
    propertyDictionary.printYourself();
}


//
// store & restore context - material info in gp not saved now!
//

contextIOResultType
Material :: saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp)
//
// saves full material status (saves state variables, that completely describe
// current state) stored in gp->matstatusDict with key = this->giveNumber()
// storing of corresponding context if it is defined for current material in
// gp status dictionary should be performed here by overloading this function.
// (such code should invoke also corresponding function for yield conditions,
// submaterials and so on)
//

//
{
    contextIOResultType iores;

    if ( gp == NULL ) {
        THROW_CIOERR(CIO_BADOBJ);
    }

    // write raw data - we save status there for this
    MaterialStatus *status = this->giveStatus(gp);

    if ( status ) {
        if ( ( iores = status->saveContext(stream, mode, gp) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}

contextIOResultType
Material :: restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp)
//
// restores full material status (saves state variables, that completely describe
// current state) stored in gp->matstatusDict with key = this->giveNumber()
// restoring of corresponding context if it is defined for current material in
// gp status dictionary should be performed here by overloading this function.
// (such code should invoke also corresponding function for yield conditions,
//  submaterials and so on)
//

//
{
    contextIOResultType iores;
    if ( gp == NULL ) {
        THROW_CIOERR(CIO_BADOBJ);
    }

    // read raw data - context
    MaterialStatus *status =  this->giveStatus(gp);
    if ( status ) {
        if ( ( iores = status->restoreContext(stream, mode, gp) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}



MaterialStatus *
Material :: giveStatus(GaussPoint *gp) const
/*
 * returns material status in gp corresponding to specific material class
 */
{
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);

        // if newly created status is null
        // dont include it. specific instance
        // does not have status.
        if ( status != NULL ) {
            gp->setMaterialStatus( status, this->giveNumber() );
        }
    }

    return status;
}


void
Material :: initTempStatus(GaussPoint *gp)
//
// Initialize MatStatus (respective it's temporary variables at the begining
// of integrating incremental constitutive relations) to correct values
//
{
    MaterialStatus *status = this->giveStatus(gp);
    if ( status ) {
        status->initTempStatus();
    }
}


int
Material :: initMaterial(Element *element)
{
    return 0;
}
} // end namespace oofem
