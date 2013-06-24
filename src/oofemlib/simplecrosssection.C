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

#include "simplecrosssection.h"
#include "gausspoint.h"
#include "structuralmaterial.h"
#include "floatarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {

REGISTER_CrossSection( SimpleCrossSection );

#if 0 //@todo don't see any difference from base class /JB
void
SimpleCrossSection :: giveRealStresses(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                       const FloatArray &totalStrain, TimeStep *tStep)
//
// this function returns a real stresses corresponding to
// given totalStrain according to stressStrain mode stored
// in each gp.
// IMPORTANT:
//
{
    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );

    if ( mat->hasMaterialModeCapability(mode) ) {
        StructuralCrossSection :: giveRealStresses(answer, form, gp, totalStrain, tStep);
        return;
    } else {
        _error("giveRealStresses : unsupported mode");
    }
}
#endif

void
SimpleCrossSection :: giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                        MatResponseForm form,
                                                        MatResponseMode rMode,
                                                        GaussPoint *gp,
                                                        StructuralMaterial *mat,
                                                        TimeStep *tStep)
//
// only interface to material class, forcing returned matrix to be in reduced form.
//
{
    this->giveMaterialStiffnessMatrixOf(answer, form, rMode, gp, mat, tStep);
}


void
SimpleCrossSection :: giveMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    MatResponseMode rMode,
                                                    GaussPoint *gp,
                                                    StructuralMaterial *mat,
                                                    TimeStep *tStep)
//
// may be only simple interface to material class, forcing returned matrix to be in reduced form.
// otherwise special methods called to obtain required stiffness from 3d case.
//
{
    if ( mat->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
        mat->giveCharacteristicMatrix(answer, form, rMode, gp, tStep);
        return;
    } else {
        OOFEM_ERROR3("GiveMaterialStiffnessMatrixOf: unsupported StressStrainMode %s on Element number %d", __MaterialModeToString(gp->giveMaterialMode()), gp->giveElement()->giveGlobalNumber() );
    }
}


IRResultType
SimpleCrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    double value;

    this->CrossSection :: initializeFrom(ir);

    double thick = 0.0;
    if ( ir->hasField(_IFT_SimpleCrossSection_thick)) {
        IR_GIVE_OPTIONAL_FIELD(ir, thick, _IFT_SimpleCrossSection_thick);
        propertyDictionary->add(CS_Thickness, thick);
    }

    double width = 0.0;
    if ( ir->hasField(_IFT_SimpleCrossSection_width)) {
        IR_GIVE_OPTIONAL_FIELD(ir, width, _IFT_SimpleCrossSection_width);
        propertyDictionary->add(CS_Width, width);
    }

    double area = 0.0;
    if ( ir->hasField(_IFT_SimpleCrossSection_area) ) {
        IR_GIVE_FIELD(ir, area, _IFT_SimpleCrossSection_area);
    } else {
        area = thick*width;
    }
    propertyDictionary->add(CS_Area, area);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_iy);
    propertyDictionary->add(CS_InertiaMomentY, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_iz);
    propertyDictionary->add(CS_InertiaMomentZ, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_ik);
    propertyDictionary->add(CS_TorsionMomentX, value);

    double beamshearcoeff=0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, beamshearcoeff, _IFT_SimpleCrossSection_shearcoeff);
    propertyDictionary->add(CS_BeamShearCoeff, beamshearcoeff);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_shearareay);
    if (value == 0.0) value=beamshearcoeff * area;
    propertyDictionary->add(CS_SHEAR_AREA_Y, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_shearareaz);
    if (value == 0.0) value=beamshearcoeff * area;
    propertyDictionary->add(CS_SHEAR_AREA_Z, value);

    return IRRT_OK;
}


void SimpleCrossSection::giveInputRecord(DynamicInputRecord &input)
{
    StructuralCrossSection :: giveInputRecord(input);
    input.setField(this->give(CS_Thickness), _IFT_SimpleCrossSection_thick);
    input.setField(this->give(CS_Width), _IFT_SimpleCrossSection_width);
    input.setField(this->give(CS_Area), _IFT_SimpleCrossSection_area);
    input.setField(this->give(CS_TorsionMomentX), _IFT_SimpleCrossSection_ik);
    input.setField(this->give(CS_InertiaMomentY), _IFT_SimpleCrossSection_iy);
    input.setField(this->give(CS_InertiaMomentZ), _IFT_SimpleCrossSection_iz);
    input.setField(this->give(CS_SHEAR_AREA_Y), _IFT_SimpleCrossSection_shearareay);
    input.setField(this->give(CS_SHEAR_AREA_Y), _IFT_SimpleCrossSection_shearareaz);
    input.setField(this->give(CS_BeamShearCoeff), _IFT_SimpleCrossSection_shearcoeff);
}


double
SimpleCrossSection :: give(CrossSectionProperty aProperty)
{
    double value = 0.0;

    if ( propertyDictionary->includes(aProperty) ) {
        value = propertyDictionary->at(aProperty);
    } else {
        OOFEM_ERROR3("Simple cross-section Number %d has undefined property ID %d", this->giveNumber(), aProperty);
    }

    return value;
}


void
SimpleCrossSection :: computeStressIndependentStrainVector(FloatArray &answer,
                                                           GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
//
// returns initial strain vector induced by stress independent effects
// like temperatue or shrinkage.
// takes into account form of load vector assumed by engngModel (Incremental or Total Load form).
//
{
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );
    ///@todo  Deprecated or not? If so, remove it! / Mikael
#if 0
    MaterialMode matmode = gp-> giveMaterialMode ();
    FloatArray et, e0, fullAnswer;
    double thick, width;

    if ((matmode == _2dBeam) || (matmode == _3dBeam) || (matmode == _3dShell) || (matmode == _2dPlate)) {

        StructuralElement *elem = (StructuralElement*)gp->giveElement();
        elem -> computeResultingIPTemperatureAt (et, stepN, gp, mode);
        FloatArray redAnswer;

        if (et.giveSize() == 0) {answer.resize(0); return ;}
        if (et.giveSize() < 1) {
            _error ("computeStressIndependentStrainVector - Bad format of TemperatureLoad");
            exit (1);
        }
        mat->giveThermalDilatationVector (e0, gp,stepN);

        if (matmode == _2dBeam) {
            answer.resize (3);
            answer.zero();
            answer.at(1) = e0.at(1) * (et.at(1)- mat->giveReferenceTemperature());
            if (et.giveSize() > 1) {
                thick = this->give(THICKNESS);
                answer.at(2) = e0.at(1) * et.at(2)/ thick;   // kappa_x
            }
        } else if (matmode == _3dBeam) {
            answer.resize (6);
            answer.zero();

            answer.at(1) = e0.at(1) * (et.at(1)- mat->giveReferenceTemperature());
            if (et.giveSize() > 1) {
                thick = this->give(THICKNESS);
                width = this->give(WIDTH);
                answer.at(5) = e0.at(1) * et.at(2)/ thick;   // kappa_y
                if (et.giveSize() > 2)
                    answer.at(6) = e0.at(1) * et.at(3)/ width;   // kappa_z
            }
        } else if (matmode == _2dPlate) {

            if (et.giveSize() > 1) {
                answer.resize (5);
                answer.zero();

                thick = this->give(THICKNESS);
                if (et.giveSize() > 1) {
                    answer.at(1) = e0.at(1) * et.at(2)/ thick;   // kappa_x
                    answer.at(2) = e0.at(2) * et.at(2)/ thick;   // kappa_y
                }
            }
        } else if (matmode == _3dShell) {
            answer.resize (8);
            answer.zero();

            answer.at(1) = e0.at(1) * (et.at(1)- mat->giveReferenceTemperature());
            answer.at(2) = e0.at(2) * (et.at(1)- mat->giveReferenceTemperature());
            if (et.giveSize() > 1) {
                thick = this->give(THICKNESS);
                answer.at(4) = e0.at(1) * et.at(2)/ thick;   // kappa_x
                answer.at(5) = e0.at(2) * et.at(2)/ thick;   // kappa_y
            }
        } else _error ("Unsupported material mode");
    } else {
        mat->computeStressIndependentStrainVector (answer, gp, stepN, mode);
    }
#endif
    mat->computeStressIndependentStrainVector(answer, gp, stepN, mode);
}
} // end namespace oofem
