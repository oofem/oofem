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

#include "../sm/CrossSections/simplecrosssection.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Elements/structuralelement.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "engngm.h"

namespace oofem {
REGISTER_CrossSection(SimpleCrossSection);


void
SimpleCrossSection :: giveRealStress_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->giveRealStressVector_3d(answer, gp, strain, tStep);
}

void
SimpleCrossSection :: giveRealStress_3dDegeneratedShell(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    IntArray strainControl = {
        1, 2, 4, 5, 6
    };
    mat->giveRealStressVector_ShellStressControl(answer, gp, strain, strainControl, tStep);
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
SimpleCrossSection :: giveRealStress_Warping(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->giveRealStressVector_Warping(answer, gp, strain, tStep);
}


void
SimpleCrossSection :: giveStiffnessMatrix_3d(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->give3dMaterialStiffnessMatrix(answer, rMode, gp, tStep);
}


void
SimpleCrossSection :: giveStiffnessMatrix_PlaneStress(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->givePlaneStressStiffMtrx(answer, rMode, gp, tStep);
}


void
SimpleCrossSection :: giveStiffnessMatrix_PlaneStrain(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->givePlaneStrainStiffMtrx(answer, rMode, gp, tStep);
}


void
SimpleCrossSection :: giveStiffnessMatrix_1d(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->give1dStressStiffMtrx(answer, rMode, gp, tStep);
}


void
SimpleCrossSection :: giveGeneralizedStress_Beam2d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    /**Note: (by bp): This assumes that the behaviour is elastic
     * there exist a nuumber of nonlinear integral material models for beams defined directly in terms of integral forces and moments and corresponding
     * deformations and curvatures. This would require to implement support at material model level.
     * Mikael: That would not be a continuum material model, but it would highly depend on the cross-section shape, thus, it should be a special cross-section model instead.
     * This cross-section assumes you can split the response into inertia moments and pure material response. This is only possible for a constant, elastic response (i.e. elastic).
     * Therefore, this cross-section may only be allowed to give the elastic response.
     */
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    FloatArray elasticStrain, et, e0;
    FloatMatrix tangent;
    elasticStrain = strain;
    this->giveTemperatureVector(et, gp, tStep);
    if ( et.giveSize() > 0 ) {
        mat->giveThermalDilatationVector(e0, gp, tStep);
        double thick = this->give(CS_Thickness, gp);
        elasticStrain.at(1) -= e0.at(1) * ( et.at(1) - mat->giveReferenceTemperature() );
        if ( et.giveSize() > 1 ) {
            elasticStrain.at(2) -= e0.at(1) * et.at(2) / thick;     // kappa_x
        }
    }
    this->give2dBeamStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, elasticStrain);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
SimpleCrossSection :: giveGeneralizedStress_Beam3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    /**Note: (by bp): This assumes that the behaviour is elastic
     * there exist a nuumber of nonlinear integral material models for beams defined directly in terms of integral forces and moments and corresponding
     * deformations and curvatures. This would require to implement support at material model level.
     * Mikael: See earlier response to comment
     */
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    FloatArray elasticStrain, et, e0;
    FloatMatrix tangent;
    elasticStrain = strain;
    this->giveTemperatureVector(et, gp, tStep);
    if ( et.giveSize() > 0 ) {
        double thick = this->give(CS_Thickness, gp);
        double width = this->give(CS_Width, gp);
        mat->giveThermalDilatationVector(e0, gp, tStep);
        elasticStrain.at(1) -= e0.at(1) * ( et.at(1) - mat->giveReferenceTemperature() );
        if ( et.giveSize() > 1 ) {
            elasticStrain.at(5) -= e0.at(1) * et.at(2) / thick;     // kappa_y
            if ( et.giveSize() > 2 ) {
                elasticStrain.at(6) -= e0.at(1) * et.at(3) / width;     // kappa_z
            }
        }
    }
    this->give3dBeamStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, elasticStrain);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
SimpleCrossSection :: giveGeneralizedStress_Plate(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    /**Note: (by bp): This assumes that the behaviour is elastic
     * there exist a nuumber of nonlinear integral material models for beams/plates/shells
     * defined directly in terms of integral forces and moments and corresponding
     * deformations and curvatures. This would require to implement support at material model level.
     * Mikael: See earlier response to comment
     */
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    FloatArray elasticStrain, et, e0;
    FloatMatrix tangent;
    elasticStrain = strain;
    this->giveTemperatureVector(et, gp, tStep);
    if ( et.giveSize() > 0 ) {
        double thick = this->give(CS_Thickness, gp);
        mat->giveThermalDilatationVector(e0, gp, tStep);
        if ( et.giveSize() > 1 ) {
            elasticStrain.at(1) -= e0.at(1) * et.at(2) / thick;     // kappa_x
            elasticStrain.at(2) -= e0.at(2) * et.at(2) / thick;     // kappa_y
        }
    }
    this->give2dPlateStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, elasticStrain);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
SimpleCrossSection :: giveGeneralizedStress_Shell(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    /**Note: (by bp): This assumes that the behaviour is elastic
     * there exist a nuumber of nonlinear integral material models for beams/plates/shells
     * defined directly in terms of integral forces and moments and corresponding
     * deformations and curvatures. This would require to implement support at material model level.
     * Mikael: See earlier response to comment
     */
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    FloatArray elasticStrain, et, e0;
    FloatMatrix tangent;
    elasticStrain = strain;
    this->giveTemperatureVector(et, gp, tStep);
    if ( et.giveSize() ) {
        double thick = this->give(CS_Thickness, gp);
        mat->giveThermalDilatationVector(e0, gp, tStep);
        elasticStrain.at(1) -= e0.at(1) * ( et.at(1) - mat->giveReferenceTemperature() );
        elasticStrain.at(2) -= e0.at(2) * ( et.at(1) - mat->giveReferenceTemperature() );
        if ( et.giveSize() > 1 ) {
            elasticStrain.at(4) -= e0.at(1) * et.at(2) / thick;     // kappa_x
            elasticStrain.at(5) -= e0.at(2) * et.at(2) / thick;     // kappa_y
        }
    }
    this->give3dShellStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, elasticStrain);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
SimpleCrossSection :: giveGeneralizedStress_MembraneRot(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    FloatMatrix tangent;
    this->giveMembraneRotStiffMtrx(tangent, ElasticStiffness, gp, tStep);
    answer.beProductOf(tangent, strain);

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveMaterial(gp)->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);


    ///@todo We should support nonlinear behavior for the membrane part. In fact, should be even bundle the rotation part with the membrane?
    /// We gain nothing from this design anyway as the rotation field is always separate. Separate manual integration by the element would be an option.
}

void
SimpleCrossSection :: giveGeneralizedStress_PlateSubSoil(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep)
{
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->giveRealStressVector_2dPlateSubSoil(answer, gp, generalizedStrain, tStep);
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
    } else if ( mode == _3dDegeneratedShell ) {
        this->give3dDegeneratedShellStiffMtrx(answer, rMode, gp, tStep);
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
    area = this->give(CS_Area, gp);
    Iy   = this->give(CS_InertiaMomentY, gp);
    shearAreaz = this->give(CS_ShearAreaZ, gp);

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
    area = this->give(CS_Area, gp);
    Iy   = this->give(CS_InertiaMomentY, gp);
    Iz   = this->give(CS_InertiaMomentZ, gp);
    Ik   = this->give(CS_TorsionMomentX, gp);

    //shearCoeff = this->give(CS_BeamShearCoeff);
    shearAreay = this->give(CS_ShearAreaY, gp);
    shearAreaz = this->give(CS_ShearAreaZ, gp);

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
    thickness = this->give(CS_Thickness, gp);
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

    double thickness = this->give(CS_Thickness, gp);
    double thickness3 = thickness * thickness * thickness;

    mat->givePlaneStressStiffMtrx(mat3d, rMode, gp, tStep);

    answer.resize(8, 8);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = mat3d.at(i, j) * thickness;
        }
    }
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i + 3, j + 3) = mat3d.at(i, j) * thickness3 / 12.0;
        }
    }

    answer.at(8, 8) = answer.at(7, 7) = mat3d.at(3, 3) * thickness * ( 5. / 6. );
}


void
SimpleCrossSection :: give3dDegeneratedShellStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->give3dMaterialStiffnessMatrix(answer, rMode, gp, tStep);

    answer.at(1, 1) -= answer.at(1, 3) * answer.at(3, 1) / answer.at(3, 3);
    answer.at(2, 1) -= answer.at(2, 3) * answer.at(3, 1) / answer.at(3, 3);
    answer.at(1, 2) -= answer.at(1, 3) * answer.at(3, 2) / answer.at(3, 3);
    answer.at(2, 2) -= answer.at(2, 3) * answer.at(3, 2) / answer.at(3, 3);

    answer.at(3, 1) = 0.0;
    answer.at(3, 2) = 0.0;
    answer.at(3, 3) = 0.0;
    answer.at(2, 3) = 0.0;
    answer.at(1, 3) = 0.0;
}

void
SimpleCrossSection :: giveMembraneRotStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat;
    mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->givePlaneStressStiffMtrx(answer, ElasticStiffness, gp, tStep);
    answer.resizeWithData(4, 4);
    answer.at(4, 4) = this->give(CS_DrillingStiffness, gp);
}

void
SimpleCrossSection :: give2dPlateSubSoilStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterial *mat;
    mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    mat->give2dPlateSubSoilStiffMtrx(answer, ElasticStiffness, gp, tStep);
}


IRResultType
SimpleCrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    double value;

    double thick = 0.0;
    if ( ir->hasField(_IFT_SimpleCrossSection_thick) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, thick, _IFT_SimpleCrossSection_thick);
        propertyDictionary.add(CS_Thickness, thick);
    }

    double width = 0.0;
    if ( ir->hasField(_IFT_SimpleCrossSection_width) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, width, _IFT_SimpleCrossSection_width);
        propertyDictionary.add(CS_Width, width);
    }

    double area = 0.0;
    if ( ir->hasField(_IFT_SimpleCrossSection_area) ) {
        IR_GIVE_FIELD(ir, area, _IFT_SimpleCrossSection_area);
    } else {
        area = thick * width;
    }
    propertyDictionary.add(CS_Area, area);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_iy);
    propertyDictionary.add(CS_InertiaMomentY, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_iz);
    propertyDictionary.add(CS_InertiaMomentZ, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_ik);
    propertyDictionary.add(CS_TorsionMomentX, value);

    double beamshearcoeff = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, beamshearcoeff, _IFT_SimpleCrossSection_shearcoeff);
    propertyDictionary.add(CS_BeamShearCoeff, beamshearcoeff);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_shearareay);
    if ( value == 0.0 ) {
        value = beamshearcoeff * area;
    }
    propertyDictionary.add(CS_ShearAreaY, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_shearareaz);
    if ( value == 0.0 ) {
        value = beamshearcoeff * area;
    }
    propertyDictionary.add(CS_ShearAreaZ, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_drillStiffness);
    propertyDictionary.add(CS_DrillingStiffness, value);

    this->materialNumber = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->materialNumber, _IFT_SimpleCrossSection_MaterialNumber);

    if ( ir->hasField(_IFT_SimpleCrossSection_directorx) ) {
        value = 0.0;
        IR_GIVE_FIELD(ir, value, _IFT_SimpleCrossSection_directorx);
        propertyDictionary.add(CS_DirectorVectorX, value);

        value = 0.0;
        IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_directory);
        propertyDictionary.add(CS_DirectorVectorY, value);

        value = 1.0;
        IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_directorz);
        propertyDictionary.add(CS_DirectorVectorZ, value);
    }

    return CrossSection :: initializeFrom(ir);
}


void
SimpleCrossSection :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralCrossSection :: giveInputRecord(input);

    if ( this->propertyDictionary.includes(CS_Thickness) ) {
        input.setField(this->propertyDictionary.at(CS_Thickness), _IFT_SimpleCrossSection_thick);
    }

    if ( this->propertyDictionary.includes(CS_Width) ) {
        input.setField(this->propertyDictionary.at(CS_Width), _IFT_SimpleCrossSection_width);
    }

    if ( this->propertyDictionary.includes(CS_Area) ) {
        input.setField(this->propertyDictionary.at(CS_Area), _IFT_SimpleCrossSection_area);
    }

    if ( this->propertyDictionary.includes(CS_TorsionMomentX) ) {
        input.setField(this->propertyDictionary.at(CS_TorsionMomentX), _IFT_SimpleCrossSection_ik);
    }

    if ( this->propertyDictionary.includes(CS_InertiaMomentY) ) {
        input.setField(this->propertyDictionary.at(CS_InertiaMomentY), _IFT_SimpleCrossSection_iy);
    }

    if ( this->propertyDictionary.includes(CS_InertiaMomentZ) ) {
        input.setField(this->propertyDictionary.at(CS_InertiaMomentZ), _IFT_SimpleCrossSection_iz);
    }

    if ( this->propertyDictionary.includes(CS_ShearAreaY) ) {
        input.setField(this->propertyDictionary.at(CS_ShearAreaY), _IFT_SimpleCrossSection_shearareay);
    }

    if ( this->propertyDictionary.includes(CS_ShearAreaZ) ) {
        input.setField(this->propertyDictionary.at(CS_ShearAreaZ), _IFT_SimpleCrossSection_shearareaz);
    }

    if ( this->propertyDictionary.includes(CS_BeamShearCoeff) ) {
        input.setField(this->propertyDictionary.at(CS_BeamShearCoeff), _IFT_SimpleCrossSection_shearcoeff);
    }

    input.setField(this->materialNumber, _IFT_SimpleCrossSection_MaterialNumber);

    if ( this->propertyDictionary.includes(CS_DirectorVectorX) ) {
        input.setField(this->propertyDictionary.at(CS_DirectorVectorX), _IFT_SimpleCrossSection_directorx);
    }

    if ( this->propertyDictionary.includes(CS_DirectorVectorY) ) {
        input.setField(this->propertyDictionary.at(CS_DirectorVectorY), _IFT_SimpleCrossSection_directory);
    }

    if ( this->propertyDictionary.includes(CS_DirectorVectorZ) ) {
        input.setField(this->propertyDictionary.at(CS_DirectorVectorZ), _IFT_SimpleCrossSection_directorz);
    }
}

void
SimpleCrossSection :: createMaterialStatus(GaussPoint &iGP)
{
    Material *mat = domain->giveMaterial(materialNumber);
    MaterialStatus *matStat = mat->CreateStatus(& iGP);
    iGP.setMaterialStatus(matStat);
}


bool
SimpleCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode)
{
    if ( this->giveMaterialNumber() ) {
        return this->domain->giveMaterial( this->giveMaterialNumber() )->isCharacteristicMtrxSymmetric(rMode);
    } else {
        return false; // Bet false...
    }
}


Material
*SimpleCrossSection :: giveMaterial(IntegrationPoint *ip)
{
    if ( this->giveMaterialNumber() ) {
        return this->giveDomain()->giveMaterial( this->giveMaterialNumber() );
    } else {
        return ip->giveElement()->giveMaterial();
    }
}


double
SimpleCrossSection :: give(int aProperty, GaussPoint *gp)
{
    return this->giveMaterial(gp)->give(aProperty, gp);
}


int
SimpleCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_CrossSectionNumber ) {
        answer.resize(1);
        answer.at(1) = this->giveNumber();
        return 1;
    }
    return this->giveMaterial(ip)->giveIPValue(answer, ip, type, tStep);
}



int
SimpleCrossSection :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    Material *mat = this->giveDomain()->giveMaterial(this->materialNumber);
    if ( !dynamic_cast< StructuralMaterial * >(mat) ) {
        OOFEM_WARNING( "material %s without structural support", mat->giveClassName() );
        result = 0;
    }

    return result;
}


Interface
*SimpleCrossSection :: giveMaterialInterface(InterfaceType t, IntegrationPoint *ip)
{
    return this->giveMaterial(ip)->giveInterface(t);
}




void
SimpleCrossSection :: giveFirstPKStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep)
{
    // This function returns the first Piola-Kirchoff stress in vector format
    // corresponding to a given deformation gradient according to the stress-deformation
    // mode stored in the each gp.

    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    if ( mode == _3dMat ) {
        mat->giveFirstPKStressVector_3d(answer, gp, reducedvF, tStep);
    } else if ( mode == _PlaneStrain ) {
        mat->giveFirstPKStressVector_PlaneStrain(answer, gp, reducedvF, tStep);
    } else if ( mode == _PlaneStress ) {
        mat->giveFirstPKStressVector_PlaneStress(answer, gp, reducedvF, tStep);
    } else if ( mode == _1dMat ) {
        mat->giveFirstPKStressVector_1d(answer, gp, reducedvF, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}


void
SimpleCrossSection :: giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep)
{
    // This function returns the Cauchy stress in vector format
    // corresponding to a given deformation gradient according to the stress-deformation
    // mode stored in the each gp.

    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    if ( mode == _3dMat ) {
        mat->giveCauchyStressVector_3d(answer, gp, reducedvF, tStep);
    } else if ( mode == _PlaneStrain ) {
        mat->giveCauchyStressVector_PlaneStrain(answer, gp, reducedvF, tStep);
    } else if ( mode == _PlaneStress ) {
        mat->giveCauchyStressVector_PlaneStress(answer, gp, reducedvF, tStep);
    } else if ( mode == _1dMat ) {
        mat->giveCauchyStressVector_1d(answer, gp, reducedvF, tStep);
    }
}

void
SimpleCrossSection :: giveEshelbyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    if ( mode == _3dMat ) {
        //        mat->giveCauchyStressVector_3d(answer, gp, reducedvF, tStep);
    } else if ( mode == _PlaneStrain ) {
        mat->giveEshelbyStressVector_PlaneStrain(answer, gp, reducedvF, tStep);
    } else if ( mode == _PlaneStress ) {
        //        mat->giveCauchyStressVector_PlaneStress(answer, gp, reducedvF, tStep);
    } else if ( mode == _1dMat ) {
        //        mat->giveCauchyStressVector_1d(answer, gp, reducedvF, tStep);
    }
}


void
SimpleCrossSection :: giveStiffnessMatrix_dPdF(FloatMatrix &answer,
                                               MatResponseMode rMode, GaussPoint *gp,
                                               TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        mat->give3dMaterialStiffnessMatrix_dPdF(answer, rMode, gp, tStep);
    } else if ( mode == _PlaneStress ) {
        mat->givePlaneStressStiffMtrx_dPdF(answer, rMode, gp, tStep);
    } else if ( mode == _PlaneStrain ) {
        mat->givePlaneStrainStiffMtrx_dPdF(answer, rMode, gp, tStep);
    } else if ( mode == _1dMat ) {
        mat->give1dStressStiffMtrx_dPdF(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}


void
SimpleCrossSection :: giveStiffnessMatrix_dCde(FloatMatrix &answer,
                                               MatResponseMode rMode, GaussPoint *gp,
                                               TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        mat->give3dMaterialStiffnessMatrix_dCde(answer, rMode, gp, tStep);
    } else if ( mode == _PlaneStress ) {
        mat->givePlaneStressStiffMtrx_dCde(answer, rMode, gp, tStep);
    } else if ( mode == _PlaneStrain ) {
        mat->givePlaneStrainStiffMtrx_dCde(answer, rMode, gp, tStep);
    } else if ( mode == _1dMat ) {
        mat->give1dStressStiffMtrx_dCde(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}


void
SimpleCrossSection :: giveTemperatureVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    Element *elem = gp->giveElement();
    answer.clear();
    //sum up all prescribed temperatures over an element
    StructuralElement *selem = dynamic_cast< StructuralElement * >(elem);
    selem->computeResultingIPTemperatureAt(answer, tStep, gp, VM_Total);

    /* add external source, if provided */
    FieldManager *fm = this->domain->giveEngngModel()->giveContext()->giveFieldManager();
    FieldPtr tf;

    if ( ( tf = fm->giveField(FT_Temperature) ) ) {
        // temperature field registered
        FloatArray gcoords, et2;
        int err;
        elem->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates() );
        if ( ( err = tf->evaluateAt(et2, gcoords, VM_Total, tStep) ) ) {
            OOFEM_ERROR("tf->evaluateAt failed, element %d, error code %d", elem->giveNumber(), err);
        }
        if ( et2.isNotEmpty() ) {
            if ( answer.isEmpty() ) {
                answer = et2;
            } else {
                answer.at(1) += et2.at(1);
            }
        }
    }
}

int
SimpleCrossSection :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveMaterial(gp)->packUnknowns(buff, tStep, gp);
}

int
SimpleCrossSection :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveMaterial(gp)->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
SimpleCrossSection :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->giveMaterial(gp)->estimatePackSize(buff, gp);
}
} // end namespace oofem
