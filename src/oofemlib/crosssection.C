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
#include "dynamicinputrecord.h"
#include "gausspoint.h"
#include "material.h"
#include "contextioerr.h"
#include "datastream.h"
#include "gaussintegrationrule.h"

namespace oofem {
    
CrossSection :: CrossSection(int n, Domain* d) : FEMComponent(n, d), propertyDictionary(), setNumber(0)
{
}

CrossSection :: ~CrossSection()
{
}

int
CrossSection :: setupIntegrationPoints(IntegrationRule &irule, int npoints, Element *element)
{
    return irule.setUpIntegrationPoints( element->giveIntegrationDomain(), npoints, element->giveMaterialMode() );
}


int
CrossSection :: setupIntegrationPoints(IntegrationRule &irule, int npointsXY, int npointsZ, Element *element)
{
    return irule.setUpIntegrationPoints( element->giveIntegrationDomain(), npointsXY, npointsZ, element->giveMaterialMode() );
}

void
CrossSection :: initializeFrom(InputRecord &ir)
{
    // Read set number the cross section is applied to
    this->setNumber = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->setNumber, _IFT_CrossSection_SetNumber);
}

void
CrossSection :: giveInputRecord(DynamicInputRecord &input)
{
	FEMComponent :: giveInputRecord(input);
    input.setField(this->setNumber, _IFT_CrossSection_SetNumber);
}

int
CrossSection :: giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_CrossSectionNumber ) {
        answer.resize(1);
        answer.at(1) = this->giveNumber();
        return 1;
    }
    return ip->giveMaterial()->giveIPValue(answer, ip, type, tStep);
}


void
CrossSection :: printYourself()
// Prints the receiver on screen.
{
    printf("Cross Section with properties : \n");
    propertyDictionary.printYourself();
}


void
CrossSection :: saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp)
{
    Material *mat = this->giveMaterial(gp);
    mat->saveIPContext(stream, mode, gp);
}


void
CrossSection :: restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp)
{
    Material *mat = this->giveMaterial(gp);
    mat->restoreIPContext(stream, mode, gp);
}

bool
CrossSection :: hasProperty(CrossSectionProperty aProperty)
{
    return propertyDictionary.includes(aProperty);
}

double
CrossSection :: give(CrossSectionProperty aProperty, GaussPoint *gp) const
// Returns the value of the property aProperty of the receiver.
{
    if ( propertyDictionary.includes(aProperty) ) {
        return propertyDictionary.at(aProperty);
    } else {
        OOFEM_ERROR("Undefined property ID %d", aProperty);
    }

    return 0.0;
}

double
CrossSection :: give(CrossSectionProperty aProperty, const FloatArray &coords, Element *elem, bool local) const
// Returns the value of the property aProperty of the receiver.
{
    if ( propertyDictionary.includes(aProperty) ) {
        return propertyDictionary.at(aProperty);
    } else {
        OOFEM_ERROR("Undefined property ID %d", aProperty);
    }

    return 0.0;
}


double
CrossSection :: predictRelativeComputationalCost(GaussPoint *gp)
{
    return this->giveRelativeSelfComputationalCost() * this->giveMaterial(gp)->predictRelativeComputationalCost(gp);
}


void CrossSection :: saveContext(DataStream &stream, ContextMode mode)
{
    FEMComponent :: saveContext(stream, mode);

    if ( ( mode & CM_Definition ) ) {
        propertyDictionary.saveContext (stream) ;

        if ( !stream.write(setNumber) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}


void CrossSection :: restoreContext(DataStream &stream, ContextMode mode)
{
    FEMComponent :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        propertyDictionary.restoreContext(stream);

        if ( !stream.read(setNumber) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}

} // end namespace oofem
