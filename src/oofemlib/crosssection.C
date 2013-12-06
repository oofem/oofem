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

#include "crosssection.h"
#include "dictionary.h"
#include "gausspoint.h"
#include "material.h"
#include "contextioerr.h"
#include "gaussintegrationrule.h"



namespace oofem {

int
CrossSection :: setupIntegrationPoints(IntegrationRule &irule, int npoints, Element *element)
{
    return irule.setUpIntegrationPoints(element->giveIntegrationDomain(), npoints, element->giveMaterialMode());
}


IRResultType
CrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    // Read set number the cross section is applied to
    IR_GIVE_OPTIONAL_FIELD(ir, this->setNumber, _IFT_CrossSection_SetNumber );

    return IRRT_OK;
}


void
CrossSection :: printYourself()
// Prints the receiver on screen.
{
    printf("Cross Section with properties : \n");
    propertyDictionary->printYourself();
}


contextIOResultType
CrossSection :: saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp)
//
// saves full material context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    Material *mat = gp->giveMaterial();

    if ( ( iores = mat->saveIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
CrossSection :: restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    Material *mat = gp->giveMaterial();

    if ( ( iores = mat->restoreIPContext(stream, mode, gp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


double
CrossSection :: give(CrossSectionProperty aProperty, GaussPoint* gp)
// Returns the value of the property aProperty of the receiver.
{
    if ( propertyDictionary->includes(aProperty) ) {
        return propertyDictionary->at(aProperty);
    } else {
        OOFEM_ERROR3("Cross-section Number %d has undefined property ID %d", this->giveNumber(), aProperty);
    }

    return 0.0;
}

double
CrossSection :: give(CrossSectionProperty aProperty, const FloatArray* lc, const FloatArray* gc, Element* elem)
// Returns the value of the property aProperty of the receiver.
{
    if ( propertyDictionary->includes(aProperty) ) {
        return propertyDictionary->at(aProperty);
    } else {
        OOFEM_ERROR3("Cross-section Number %d has undefined property ID %d", this->giveNumber(), aProperty);
    }

    return 0.0;
}


#ifdef __PARALLEL_MODE
double
CrossSection :: predictRelativeComputationalCost(GaussPoint *gp)
{
    return this->giveRelativeSelfComputationalCost() * gp->giveMaterial()->predictRelativeComputationalCost(gp);
}
#endif
} // end namespace oofem
