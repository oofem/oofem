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

#include "../sm/CrossSections/variablecrosssection.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "scalarfunction.h"
#include "function.h"

#include <string>
#include <sstream>

namespace oofem {
REGISTER_CrossSection(VariableCrossSection);


IRResultType
VariableCrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    if ( ir->hasField(_IFT_SimpleCrossSection_thick) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, thicknessExpr, _IFT_SimpleCrossSection_thick);
    }

    if ( ir->hasField(_IFT_SimpleCrossSection_width) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, widthExpr, _IFT_SimpleCrossSection_width);
    }

    if ( ir->hasField(_IFT_SimpleCrossSection_area) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, areaExpr, _IFT_SimpleCrossSection_area);
    }

    if ( ir->hasField(_IFT_SimpleCrossSection_iy) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, iyExpr, _IFT_SimpleCrossSection_iy);
    }

    if ( ir->hasField(_IFT_SimpleCrossSection_iy) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, izExpr, _IFT_SimpleCrossSection_iz);
    }

    if ( ir->hasField(_IFT_SimpleCrossSection_ik) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, ixExpr, _IFT_SimpleCrossSection_ik);
    }

    if ( ir->hasField(_IFT_SimpleCrossSection_shearareay) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, shearAreayExpr, _IFT_SimpleCrossSection_shearareay);
    }

    if ( ir->hasField(_IFT_SimpleCrossSection_shearareaz) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, shearAreazExpr, _IFT_SimpleCrossSection_shearareaz);
    }

    if ( ir->hasField(_IFT_SimpleCrossSection_drillStiffness) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, drillingStiffnessExpr, _IFT_SimpleCrossSection_drillStiffness);
    }

    IR_GIVE_OPTIONAL_FIELD(ir, this->materialNumber, _IFT_SimpleCrossSection_MaterialNumber);

    directorxExpr.setValue(0.0);
    if ( ir->hasField(_IFT_SimpleCrossSection_directorx) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, directorxExpr, _IFT_SimpleCrossSection_directorx);
    }


    directoryExpr.setValue(0.0);
    if ( ir->hasField(_IFT_SimpleCrossSection_directory) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, directoryExpr, _IFT_SimpleCrossSection_directory);
    }

    directorzExpr.setValue(1.0);
    if ( ir->hasField(_IFT_SimpleCrossSection_directorz) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, directorzExpr, _IFT_SimpleCrossSection_directorz);
    }

    // will read density and other inheritted parameters as constants
    // NOTE: do not call SimpleCrossSection here (as the parameter names are same, but different type is used here!!!!)
    return CrossSection :: initializeFrom(ir);

}


void VariableCrossSection :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralCrossSection :: giveInputRecord(input);

    input.setField(this->thicknessExpr, _IFT_SimpleCrossSection_thick);
    input.setField(this->widthExpr, _IFT_SimpleCrossSection_width);
    input.setField(this->areaExpr, _IFT_SimpleCrossSection_area);
    input.setField(this->ixExpr, _IFT_SimpleCrossSection_ik);
    input.setField(this->iyExpr, _IFT_SimpleCrossSection_iy);
    input.setField(this->izExpr, _IFT_SimpleCrossSection_iz);
    input.setField(this->shearAreayExpr, _IFT_SimpleCrossSection_shearareay);
    input.setField(this->shearAreazExpr, _IFT_SimpleCrossSection_shearareaz);
    input.setField(this->drillingStiffnessExpr, _IFT_SimpleCrossSection_drillStiffness);
    input.setField(this->materialNumber, _IFT_SimpleCrossSection_MaterialNumber);
    input.setField(this->directorxExpr, _IFT_SimpleCrossSection_directorx);
    input.setField(this->directoryExpr, _IFT_SimpleCrossSection_directory);
    input.setField(this->directorzExpr, _IFT_SimpleCrossSection_directorz);
}

void
VariableCrossSection :: giveExpression(const ScalarFunction **expr, CrossSectionProperty aProperty) const
{
    if ( aProperty == CS_Thickness ) {
        * expr = & thicknessExpr;
    } else if ( aProperty ==  CS_Width ) {
        * expr = & widthExpr;
    } else if ( aProperty == CS_Area ) {
        * expr = & areaExpr;
    } else if ( aProperty == CS_TorsionMomentX ) {
        * expr = & ixExpr;
    } else if ( aProperty == CS_InertiaMomentY ) {
        * expr = & iyExpr;
    } else if ( aProperty == CS_InertiaMomentZ ) {
        * expr = & izExpr;
    } else if ( aProperty == CS_ShearAreaY ) {
        * expr = & shearAreayExpr;
    } else if ( aProperty == CS_ShearAreaZ ) {
        * expr = & shearAreazExpr;
    } else if ( aProperty == CS_DrillingStiffness ) {
        * expr = & drillingStiffnessExpr;
    } else if ( aProperty == CS_DirectorVectorX ) {
        * expr = & directorxExpr;
    } else if ( aProperty == CS_DirectorVectorY ) {
        * expr = & directoryExpr;
    } else if ( aProperty == CS_DirectorVectorZ ) {
        * expr = & directorzExpr;
    } else {
        OOFEM_ERROR("called with unknown ID %d", this->giveNumber(), aProperty);
    }
}



double
VariableCrossSection :: give(CrossSectionProperty aProperty, GaussPoint *gpx)
{
    return this->give(aProperty, gpx->giveNaturalCoordinates(), gpx->giveElement(), true);
}


double
VariableCrossSection :: give(CrossSectionProperty aProperty, const FloatArray &coords, Element *elem, bool local)
{
    double value = 0.0;
    const ScalarFunction *expr;

    if ( propertyDictionary.includes(aProperty) ) {
        value = propertyDictionary.at(aProperty);
    } else {
        this->giveExpression(& expr, aProperty);

        FloatArray c;
        if ( this->localFormulationFlag ) {
            if ( local ) {
                c = coords;
            } else {
                // convert given coords into local cs
                if ( !elem->computeLocalCoordinates(c, coords) ) {
                    OOFEM_ERROR( "computeLocalCoordinates failed (element %d)", elem->giveNumber() );
                }
            }
        } else { // global coordinates needed
            if ( local ) {
                // convert given coords into global cs
                if ( !elem->computeGlobalCoordinates(c, coords) ) {
                    OOFEM_ERROR( "computeGlobalCoordinates failed (element %d)", elem->giveNumber() );
                }
            } else {
                c = coords;
            }
        }
        // evaluate the expression
        value = expr->eval( { { "x", c } }, this->giveDomain() );
    }

    return value;
}
} // end namespace oofem
