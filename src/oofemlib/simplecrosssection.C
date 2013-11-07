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

#include "simplecrosssection.h"
#include "gausspoint.h"
#include "structuralmaterial.h"
#include "floatarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "structuralms.h"

namespace oofem {

REGISTER_CrossSection( SimpleCrossSection );


void
SimpleCrossSection :: giveRealStress_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->giveRealStressVector_3d(answer, gp, strain, tStep);
}


void
SimpleCrossSection :: giveRealStress_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
   StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
   mat->giveRealStressVector_PlaneStrain(answer, gp, strain, tStep);
}


void
SimpleCrossSection :: giveRealStress_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->giveRealStressVector_PlaneStress(answer, gp, strain, tStep);
}


void
SimpleCrossSection :: giveRealStress_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->giveRealStressVector_1d(answer, gp, strain, tStep);
}


void
SimpleCrossSection :: giveRealStress_Beam2d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    FloatMatrix tangent;
    this->give2dBeamStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, strain);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus* >( this->giveMaterial(gp)->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
SimpleCrossSection :: giveRealStress_Beam3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    FloatMatrix tangent;
    this->give3dBeamStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, strain);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus* >( this->giveMaterial(gp)->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
SimpleCrossSection :: giveRealStress_Plate(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    FloatMatrix tangent;
    this->give2dPlateStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, strain);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus* >( this->giveMaterial(gp)->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
SimpleCrossSection :: giveRealStress_Shell(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    FloatMatrix tangent;
    this->give3dShellStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, strain);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus* >( this->giveMaterial(gp)->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
SimpleCrossSection :: giveRealStress_MembraneRot(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    FloatMatrix tangent;
    this->giveMembraneRotStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, strain);

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus* >( this->giveMaterial(gp)->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
    

    ///@todo We should support nonlinear behavior for the membrane part. In fact, should be even bundle the rotation part with the membrane?
    /// We gain nothing from this design anyway as the rotation field is always separate. Separate manual integration by the element would be an option.
}


void
SimpleCrossSection :: giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _2dBeam ) {
        this->give2dBeamStiffMtrx(answer, rMode, gp, tStep);
    } else if ( mode == _3dBeam ) {
        this->give3dBeamStiffMtrx(answer, rMode, gp, tStep);
    } else if ( mode == _2dPlate ) {
        this->give2dPlateStiffMtrx(answer, rMode, gp, tStep);
    } else if ( mode == _3dShell ) {
        this->give3dShellStiffMtrx(answer, rMode, gp, tStep);
    } else {

        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
        
        if ( mode == _3dMat ) {
            mat->give3dMaterialStiffnessMatrix(answer, rMode, gp, tStep);
        } else if ( mode == _PlaneStress ) {
            mat->givePlaneStressStiffMtrx(answer, rMode, gp, tStep);
        } else if ( mode == _PlaneStrain ) {
            mat->givePlaneStrainStiffMtrx(answer, rMode, gp, tStep);
        } else if ( mode == _1dMat ) {
            mat->give1dStressStiffMtrx(answer, rMode, gp, tStep);
        } else {
            mat->giveStiffnessMatrix(answer, rMode, gp, tStep);
        }
    }
}

void
SimpleCrossSection :: give2dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    
    FloatMatrix mat3d;
    double area, Iy, shearAreaz;

    mat->give1dStressStiffMtrx(mat3d, rMode, gp, tStep);
    area = this->give(CS_Area);
    Iy   = this->give(CS_InertiaMomentY);
    shearAreaz = this->give(CS_ShearAreaZ);

    answer.resize(3, 3);
    answer.zero();

    answer.at(1, 1) = mat3d.at(1, 1) * area;
    answer.at(2, 2) = mat3d.at(1, 1) * Iy;
    answer.at(3, 3) = shearAreaz * mat->give('G', gp);
}


void
SimpleCrossSection :: give3dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    FloatMatrix mat3d;
    double area, E, G, Iy, Iz, Ik;
    double shearAreay, shearAreaz;

    mat->give1dStressStiffMtrx(mat3d, rMode, gp, tStep);
    E    = mat3d.at(1, 1);
    G    = mat->give('G', gp);
    area = this->give(CS_Area);
    Iy   = this->give(CS_InertiaMomentY);
    Iz   = this->give(CS_InertiaMomentZ);
    Ik   = this->give(CS_TorsionMomentX);

    //shearCoeff = this->give(CS_BeamShearCoeff);
    shearAreay = this->give(CS_ShearAreaY);
    shearAreaz = this->give(CS_ShearAreaZ);

    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) = E * area;
    ///@todo Do this by using the general 3d tangent matrix instead somehow!!!
    answer.at(2, 2) = shearAreay * G;
    answer.at(3, 3) = shearAreaz * G;
    //answer.at(2, 2) = shearCoeff * G * area;
    //answer.at(3, 3) = shearCoeff * G * area;
    answer.at(4, 4) = G * Ik;
    answer.at(5, 5) = E * Iy;
    answer.at(6, 6) = E * Iz;
}


void
SimpleCrossSection :: give2dPlateStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    FloatMatrix mat3d;
    double thickness3, thickness;

    mat->givePlaneStressStiffMtrx(mat3d, rMode, gp, tStep);
    thickness = this->give(CS_Thickness);
    thickness3 = thickness * thickness * thickness;

    answer.resize(5, 5);
    answer.zero();

    for ( int i = 1; i <= 2; i++ ) {
        for ( int j = 1; j <= 2; j++ ) {
            answer.at(i, j) = mat3d.at(i, j) * thickness3 / 12.;
        }
    }

    answer.at(3, 3) = mat3d.at(3, 3) * thickness3 / 12.;
    answer.at(4, 4) = mat3d.at(3, 3) * thickness * ( 5. / 6. );
    answer.at(5, 5) = answer.at(4, 4);
}


void
SimpleCrossSection :: give3dShellStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    FloatMatrix mat3d;

    double thickness = this->give(CS_Thickness);
    double thickness3 = thickness * thickness * thickness;

    mat->givePlaneStressStiffMtrx(mat3d, rMode, gp, tStep);

    answer.resize(8, 8);
    answer.zero();

    for ( int i = 1; i <= 2; i++ ) {
        for ( int j = 1; j <= 2; j++ ) {
            answer.at(i, j) = mat3d.at(i, j) * thickness;
            answer.at(i + 3, j + 3) = mat3d.at(i, j) * thickness3 / 12.0;
        }
    }

    answer.at(3, 1) = mat3d.at(3, 1) * thickness;
    answer.at(3, 2) = mat3d.at(3, 2) * thickness;
    answer.at(3, 3) = mat3d.at(3, 3) * thickness;
    answer.at(6, 4) = mat3d.at(3, 1) * thickness3 / 12.0;
    answer.at(6, 5) = mat3d.at(3, 2) * thickness3 / 12.0;
    answer.at(6, 6) = mat3d.at(3, 3) * thickness3 / 12.0;

    answer.at(8, 8) = answer.at(7, 7) = mat3d.at(3, 3) * thickness * ( 5. / 6. );
}


void
SimpleCrossSection :: giveMembraneRotStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat;
    mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->givePlaneStressStiffMtrx(answer, ElasticStiffness, gp, tStep);
    answer.resizeWithData(4, 4);
    answer.at(4,4) = this->give(CS_DrillingStiffness);
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
    propertyDictionary->add(CS_ShearAreaY, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_shearareaz);
    if (value == 0.0) value=beamshearcoeff * area;
    propertyDictionary->add(CS_ShearAreaZ, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_drillStiffness);
    propertyDictionary->add(CS_DrillingStiffness, value);

    return IRRT_OK;
}


void SimpleCrossSection :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralCrossSection :: giveInputRecord(input);

    if( this->propertyDictionary->includes(CS_Thickness) ) {
    	input.setField(this->give(CS_Thickness), _IFT_SimpleCrossSection_thick);
    }

    if( this->propertyDictionary->includes(CS_Width) ) {
    	input.setField(this->give(CS_Width), _IFT_SimpleCrossSection_width);
    }

    if( this->propertyDictionary->includes(CS_Area) ) {
    	input.setField(this->give(CS_Area), _IFT_SimpleCrossSection_area);
    }

    if( this->propertyDictionary->includes(CS_TorsionMomentX) ) {
    	input.setField(this->give(CS_TorsionMomentX), _IFT_SimpleCrossSection_ik);
    }

    if( this->propertyDictionary->includes(CS_InertiaMomentY) ) {
    	input.setField(this->give(CS_InertiaMomentY), _IFT_SimpleCrossSection_iy);
    }

    if( this->propertyDictionary->includes(CS_InertiaMomentZ) ) {
    	input.setField(this->give(CS_InertiaMomentZ), _IFT_SimpleCrossSection_iz);
    }

    if( this->propertyDictionary->includes(CS_ShearAreaY) ) {
    	input.setField(this->give(CS_ShearAreaY), _IFT_SimpleCrossSection_shearareay);
    }

    if( this->propertyDictionary->includes(CS_ShearAreaY) ) {
		// TODO: Reading shearareaz and setting it to CS_ShearAreaY. Bug or feature ?! // Erik
		input.setField(this->give(CS_ShearAreaY), _IFT_SimpleCrossSection_shearareaz);
    }

    if( this->propertyDictionary->includes(CS_BeamShearCoeff) ) {
    	input.setField(this->give(CS_BeamShearCoeff), _IFT_SimpleCrossSection_shearcoeff);
    }
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


///@todo  Deprecated or not? If so, remove it! / Mikael
#if 0
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
}
#endif


bool SimpleCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode)
{
    return this->domain->giveMaterial(this->giveMaterialNumber())->isCharacteristicMtrxSymmetric(rMode);
}

Material 
*SimpleCrossSection :: giveMaterial(IntegrationPoint *ip) 
{ 
    if ( ip->giveElement()->MAT_GIVEN_BY_CS ) {
        return this->giveDomain()->giveMaterial( this->giveMaterialNumber() ); 
    } else {
        return ip->giveElement()->giveMaterial();
    }
}


// JB
double
SimpleCrossSection :: give(int aProperty, GaussPoint *gp)
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
// atTime allows time dependent behavior to be taken into account
{

    return this->giveMaterial(gp)->give(aProperty, gp);
    
}

} // end namespace oofem
