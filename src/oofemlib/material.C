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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "material.h"
#include "verbose.h"
#include "gausspnt.h"
#include "flotarry.h"
#include "mathfem.h"
#include "contextioerr.h"

namespace oofem {
void
Material :: giveCharacteristicMatrix(FloatMatrix &answer,
                                     MatResponseForm form, MatResponseMode rMode,
                                     GaussPoint *gp,
                                     TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    _error("Material::giveCharacteristicMatrix is fully abstract, no implementation");
}


double
Material :: giveCharacteristicValue(MatResponseMode rMode,
                                    GaussPoint *gp,
                                    TimeStep *atTime)
//
// Returns characteristic value of the receiver
//
{
    _error("Material :: giveCharacteristicValue is purely abstract");
    return 0.0;
}


double
Material :: give(int aProperty, GaussPoint *gp)
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
// atTime allows time dependent behavior to be taken into account
{
    double value = 0.0;

    if ( propertyDictionary->includes(aProperty) ) {
        value = propertyDictionary->at(aProperty);
    } else {
        OOFEM_ERROR4( "give: property #%d on element %d and GP %d not defined", aProperty, gp->giveElement()->giveNumber(), gp->giveNumber() );
    }

    return value;
}


bool
Material :: hasProperty(int aProperty, GaussPoint *gp)
// Returns true if the aProperty is defined on a material
{
    return propertyDictionary->includes(aProperty);
}


void
Material :: modifyProperty(int aProperty, double value, GaussPoint *gp)
{
    if ( propertyDictionary->includes(aProperty) ) {
        propertyDictionary->at(aProperty) = value;
    } else {
        OOFEM_ERROR4( "modifyProperty: property #%d on element %d and GP %d not defined", aProperty, gp->giveElement()->giveNumber(), gp->giveNumber() );
    }
}


IRResultType
Material :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    double value;

#  ifdef VERBOSE
    // VERBOSE_PRINT1 ("Instanciating material ",this->giveNumber())
#  endif

    IR_GIVE_FIELD(ir, value, IFT_Material_density, "d"); // Macro
    propertyDictionary->add('d', value);

    this->castingTime = -1.e6;
    IR_GIVE_OPTIONAL_FIELD(ir, castingTime, IFT_Material_castingtime, "castingtime"); // Macro

    return IRRT_OK;
}



int
Material :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    FEMComponent :: giveInputRecordString(str, keyword);
    sprintf( buff, " d %e", propertyDictionary->at('d') );
    str += buff;

    return 1;
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
Material :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime) {
    answer.resize(0);
    return 0;
}

int
Material :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint) {
    return 0;
}

int
Material :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode) {
    answer.resize(0);
    return 0;
}


void
Material :: printYourself()
// Prints the receiver on screen.
{
    printf("Material with properties : \n");
    propertyDictionary->printYourself();
}


//
// store & restore context - material info in gp not saved now!
//

contextIOResultType
Material :: saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp)
//
// saves full material status (saves state variables, that completely describe
// current state) stored in gp->matstatusDict with key =  (int)this->giveClassID()
// storing of corresponding context if it is defined for current material in
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

    // write raw data - we save status there for this
    MaterialStatus *status =  this->giveStatus(gp);

    if ( status ) {
        if ( ( iores = status->saveContext(stream, mode, gp) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}

contextIOResultType
Material :: restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp)
//
// restores full material status (saves state variables, that completely describe
// current state) stored in gp->matstatusDict with key =  (int)this->giveClassID()
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
    MaterialStatus *status;
    status = (MaterialStatus*) gp->giveMaterialStatus(this->giveClassID());
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);

        // if newly created status is null
        // dont include it. specific instance
        // does not have status.
        if ( status != NULL ) {
	  gp->setMaterialStatus(status,this->giveClassID());
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


void
Material :: initGpForNewStep(GaussPoint *gp)
//
// initialize gp record at the beginning of new load Increment
// initialize gp status using this->initTempStatus(gp);
//
// this means: keeping Stress, Strain PlasticStrain Vectors.
// zeroing Stress, Strain PlasticStrain Increment Vectors.
//
{
    // initialize status
    this->initTempStatus(gp);
}

int
Material :: initMaterial(Element *element) {
    return 0;
}

void
Material :: updateYourself(GaussPoint *gp, TimeStep *atTime)
//
//
// We call MaterialStatus->updateYourself()
//
{
    MaterialStatus *status = this->giveStatus(gp);
    if ( status ) {
        status->updateYourself(atTime);
    }
}
} // end namespace oofem
