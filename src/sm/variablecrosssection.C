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

#include "variablecrosssection.h"
#include "gausspoint.h"
#include "structuralmaterial.h"
#include "floatarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "structuralms.h"
#include "parser.h"
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
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // will read density and other inheritted parameters as constants
    // NOTE: do not call SimpleCrossSection here (as the parameter names are same, but different type is used here!!!!)
    this->CrossSection :: initializeFrom(ir);

    this->thicknessExpr = "0.0";
    if ( ir->hasField(_IFT_SimpleCrossSection_thick) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, thicknessExpr, _IFT_SimpleCrossSection_thick);
    }

    this->widthExpr = "0.0";
    if ( ir->hasField(_IFT_SimpleCrossSection_width) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, widthExpr, _IFT_SimpleCrossSection_width);
    }

    this->areaExpr = "0.0";
    if ( ir->hasField(_IFT_SimpleCrossSection_area) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, areaExpr, _IFT_SimpleCrossSection_area);
    }

    this->iyExpr = "0.0";
    if ( ir->hasField(_IFT_SimpleCrossSection_iy) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, iyExpr, _IFT_SimpleCrossSection_iy);
    }
    this->izExpr = "0.0";
    if ( ir->hasField(_IFT_SimpleCrossSection_iy) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, izExpr, _IFT_SimpleCrossSection_iz);
    }
    this->ixExpr = "0.0";
    if ( ir->hasField(_IFT_SimpleCrossSection_ik) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, ixExpr, _IFT_SimpleCrossSection_ik);
    }
    this->shearAreayExpr = "0.0";
    if ( ir->hasField(_IFT_SimpleCrossSection_shearareay) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, shearAreayExpr, _IFT_SimpleCrossSection_shearareay);
    }
    this->shearAreazExpr = "0.0";
    if ( ir->hasField(_IFT_SimpleCrossSection_shearareaz) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, shearAreazExpr, _IFT_SimpleCrossSection_shearareaz);
    }
    this->drillingStiffnessExpr = "0.0";
    if ( ir->hasField(_IFT_SimpleCrossSection_drillStiffness) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, drillingStiffnessExpr, _IFT_SimpleCrossSection_drillStiffness);
    }

    IR_GIVE_OPTIONAL_FIELD(ir, this->materialNumber, _IFT_SimpleCrossSection_MaterialNumber);

    return IRRT_OK;
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
}

std :: string
VariableCrossSection :: giveExpression(CrossSectionProperty aProperty)
{
    if ( aProperty == CS_Thickness ) {
        return this->thicknessExpr;
    } else if ( aProperty ==  CS_Width ) {
        return this->widthExpr;
    } else if ( aProperty == CS_Area ) {
        return this->areaExpr;
    } else if ( aProperty == CS_TorsionMomentX ) {
        return this->ixExpr;
    } else if ( aProperty == CS_InertiaMomentY ) {
        return this->iyExpr;
    } else if ( aProperty == CS_InertiaMomentZ ) {
        return this->izExpr;
    } else if ( aProperty == CS_ShearAreaY ) {
        return this->shearAreayExpr;
    } else if ( aProperty == CS_ShearAreaZ ) {
        return this->shearAreazExpr;
    } else if ( aProperty == CS_DrillingStiffness ) {
        return this->drillingStiffnessExpr;
    } else {
        OOFEM_ERROR3("VariableCrossSection(%d)::give called with unknown ID %d", this->giveNumber(), aProperty);
        return std :: string();
    }
}



double
VariableCrossSection :: give(CrossSectionProperty aProperty, GaussPoint *gpx)
{
    return this->give(aProperty, gpx->giveCoordinates(), gpx->giveElement(), true);
}


double
VariableCrossSection :: give(CrossSectionProperty aProperty, const FloatArray *coords, Element *elem, bool local)
{
    double value = 0.0;
    std :: string expr;

    if ( propertyDictionary->includes(aProperty) ) {
        value = propertyDictionary->at(aProperty);
    } else {
        expr = this->giveExpression(aProperty);

        FloatArray c;
        if ( this->localFormulationFlag ) {
            if ( local ) {
                c = * coords;
            } else {
                // convert given coords into local cs
                if ( !elem->computeLocalCoordinates(c, * coords) ) {
                    OOFEM_ERROR2( "VariableCrossSection::give: computeLocalCoordinates failed (element %d)", elem->giveNumber() );
                }
            }
        } else { // global coordinates needed
            if ( local ) {
                // convert given coords into global cs
                if ( !elem->computeGlobalCoordinates(c, * coords) ) {
                    OOFEM_ERROR2( "VariableCrossSection::give: computeGlobalCoordinates failed (element %d)", elem->giveNumber() );
                }
            } else {
                c = * coords;
            }
        }

        // construct parser expression
        std :: ostringstream buff;
        int err, nsd = c.giveSize();
        const char code[] = "xyz";
        for ( int i = 0; i < nsd; i++ ) {
            buff << code [ i ] << "=" << c(i) << ";";
        }
        buff << expr;
        // evaluate the expression
        value = this->exprParser.eval(buff.str().c_str(), err);
        if ( err ) {
            OOFEM_ERROR2( "VariableCrossSection::give: parser syntax error (expr=\"%s\")", buff.str().c_str() );
        }
    }

    return value;
}
} // end namespace oofem
