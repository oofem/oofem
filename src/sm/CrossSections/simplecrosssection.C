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

#include "sm/CrossSections/simplecrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Elements/structuralelement.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "datastream.h"
#include "contextioerr.h"
#include "engngm.h"

namespace oofem {
REGISTER_CrossSection(SimpleCrossSection);


FloatArrayF<6>
SimpleCrossSection :: giveRealStress_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->giveRealStressVector_3d(strain, gp, tStep);
}

FloatArrayF<6>
SimpleCrossSection :: giveRealStress_3dDegeneratedShell(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    IntArray strainControl = {
        1, 2, 4, 5, 6
    };
    return mat->giveRealStressVector_ShellStressControl(strain, strainControl, gp, tStep);
}


FloatArrayF<4>
SimpleCrossSection :: giveRealStress_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->giveRealStressVector_PlaneStrain(strain, gp, tStep);
}


FloatArrayF<3>
SimpleCrossSection :: giveRealStress_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->giveRealStressVector_PlaneStress(strain, gp, tStep);
}


FloatArrayF<1>
SimpleCrossSection :: giveRealStress_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->giveRealStressVector_1d(strain, gp, tStep);
}


FloatArrayF<2>
SimpleCrossSection :: giveRealStress_Warping(const FloatArrayF<2> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->giveRealStressVector_Warping(strain, gp, tStep);
}


FloatMatrixF<6,6>
SimpleCrossSection :: giveStiffnessMatrix_3d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->give3dMaterialStiffnessMatrix(rMode, gp, tStep);
}


FloatMatrixF<3,3>
SimpleCrossSection :: giveStiffnessMatrix_PlaneStress(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->givePlaneStressStiffMtrx(rMode, gp, tStep);
}


FloatMatrixF<4,4>
SimpleCrossSection :: giveStiffnessMatrix_PlaneStrain(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->givePlaneStrainStiffMtrx(rMode, gp, tStep);
}


FloatMatrixF<1,1>
SimpleCrossSection :: giveStiffnessMatrix_1d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->give1dStressStiffMtrx(rMode, gp, tStep);
}


FloatArrayF<3>
SimpleCrossSection :: giveGeneralizedStress_Beam2d(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    /**Note: (by bp): This assumes that the behaviour is elastic
     * there exist a nuumber of nonlinear integral material models for beams defined directly in terms of integral forces and moments and corresponding
     * deformations and curvatures. This would require to implement support at material model level.
     * Mikael: That would not be a continuum material model, but it would highly depend on the cross-section shape, thus, it should be a special cross-section model instead.
     * This cross-section assumes you can split the response into inertia moments and pure material response. This is only possible for a constant, elastic response (i.e. elastic).
     * Therefore, this cross-section may only be allowed to give the elastic response.
     */
    auto mat = static_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    auto elasticStrain = strain;
    FloatArray et;
    this->giveTemperatureVector(et, gp, tStep);
    if ( et.giveSize() > 0 ) {
        auto e0 = mat->giveThermalDilatationVector(gp, tStep);
        double thick = this->give(CS_Thickness, gp);
        elasticStrain.at(1) -= e0.at(1) * ( et.at(1) - mat->giveReferenceTemperature() );
        if ( et.giveSize() > 1 ) {
            elasticStrain.at(2) -= e0.at(1) * et.at(2) / thick;     // kappa_x
        }
    }
    auto tangent = this->give2dBeamStiffMtrx(ElasticStiffness, gp, tStep);
    auto answer = dot(tangent, elasticStrain);

    auto status = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);

    return answer;
}


FloatArrayF<6>
SimpleCrossSection :: giveGeneralizedStress_Beam3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    /**Note: (by bp): This assumes that the behaviour is elastic
     * there exist a nuumber of nonlinear integral material models for beams defined directly in terms of integral forces and moments and corresponding
     * deformations and curvatures. This would require to implement support at material model level.
     * Mikael: See earlier response to comment
     */
    auto mat = static_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    auto elasticStrain = strain;
    FloatArray et;
    this->giveTemperatureVector(et, gp, tStep);
    if ( et.giveSize() > 0 ) {
        double thick = this->give(CS_Thickness, gp);
        double width = this->give(CS_Width, gp);
        auto e0 = mat->giveThermalDilatationVector(gp, tStep);
        elasticStrain.at(1) -= e0.at(1) * ( et.at(1) - mat->giveReferenceTemperature() );
        if ( et.giveSize() > 1 ) {
            elasticStrain.at(5) -= e0.at(1) * et.at(2) / thick;     // kappa_y
            if ( et.giveSize() > 2 ) {
                elasticStrain.at(6) -= e0.at(1) * et.at(3) / width;     // kappa_z
            }
        }
    }
    auto tangent = this->give3dBeamStiffMtrx(ElasticStiffness, gp, tStep);
    auto answer = dot(tangent, elasticStrain);

    auto status = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);

    return answer;
}


FloatArrayF<5>
SimpleCrossSection :: giveGeneralizedStress_Plate(const FloatArrayF<5> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    /**Note: (by bp): This assumes that the behaviour is elastic
     * there exist a number of nonlinear integral material models for beams/plates/shells
     * defined directly in terms of integral forces and moments and corresponding
     * deformations and curvatures. This would require to implement support at material model level.
     * Mikael: See earlier response to comment
     */
    auto mat = static_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    auto elasticStrain = strain;
    FloatArray et;
    this->giveTemperatureVector(et, gp, tStep); // FIXME use return channel, maybe fixed size/std::pair?
    if ( et.giveSize() > 0 ) {
        double thick = this->give(CS_Thickness, gp);
        auto e0 = mat->giveThermalDilatationVector(gp, tStep);
        if ( et.giveSize() > 1 ) {
            elasticStrain.at(1) -= e0.at(1) * et.at(2) / thick;     // kappa_x
            elasticStrain.at(2) -= e0.at(2) * et.at(2) / thick;     // kappa_y
        }
    }
    auto tangent = this->give2dPlateStiffMtrx(ElasticStiffness, gp, tStep);
    auto answer = dot(tangent, elasticStrain);

    auto status = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
    
    return answer;
}


FloatArrayF<8>
SimpleCrossSection :: giveGeneralizedStress_Shell(const FloatArrayF<8> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    /**Note: (by bp): This assumes that the behaviour is elastic
     * there exist a nuumber of nonlinear integral material models for beams/plates/shells
     * defined directly in terms of integral forces and moments and corresponding
     * deformations and curvatures. This would require to implement support at material model level.
     * Mikael: See earlier response to comment
     */
    auto mat = static_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    auto elasticStrain = strain;
    FloatArray et;
    this->giveTemperatureVector(et, gp, tStep);
    if ( et.giveSize() ) {
        double thick = this->give(CS_Thickness, gp);
        auto e0 = mat->giveThermalDilatationVector(gp, tStep);
        elasticStrain.at(1) -= e0.at(1) * ( et.at(1) - mat->giveReferenceTemperature() );
        elasticStrain.at(2) -= e0.at(2) * ( et.at(1) - mat->giveReferenceTemperature() );
        if ( et.giveSize() > 1 ) {
            elasticStrain.at(4) -= e0.at(1) * et.at(2) / thick;     // kappa_x
            elasticStrain.at(5) -= e0.at(2) * et.at(2) / thick;     // kappa_y
        }
    }
    auto tangent = this->give3dShellStiffMtrx(ElasticStiffness, gp, tStep);
    auto answer = dot(tangent, elasticStrain);

    auto status = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
    
    return answer;
}


FloatArrayF<4>
SimpleCrossSection :: giveGeneralizedStress_MembraneRot(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto tangent = this->giveMembraneRotStiffMtrx(ElasticStiffness, gp, tStep);
    auto answer = dot(tangent, strain);

    auto status = static_cast< StructuralMaterialStatus * >( this->giveMaterial(gp)->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);

    ///@todo We should support nonlinear behavior for the membrane part. In fact, should be even bundle the rotation part with the membrane?
    /// We gain nothing from this design anyway as the rotation field is always separate. Separate manual integration by the element would be an option.
    return answer;
}

FloatArrayF<3>
SimpleCrossSection :: giveGeneralizedStress_PlateSubSoil(const FloatArrayF<3> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = static_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->giveRealStressVector_2dPlateSubSoil(generalizedStrain, gp, tStep);
}


void
SimpleCrossSection :: giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _2dBeam ) {
        answer = this->give2dBeamStiffMtrx(rMode, gp, tStep);
    } else if ( mode == _3dBeam ) {
        answer = this->give3dBeamStiffMtrx(rMode, gp, tStep);
    } else if ( mode == _2dPlate ) {
        answer = this->give2dPlateStiffMtrx(rMode, gp, tStep);
    } else if ( mode == _3dShell ) {
        answer = this->give3dShellStiffMtrx(rMode, gp, tStep);
    } else if ( mode == _3dDegeneratedShell ) {
        answer = this->give3dDegeneratedShellStiffMtrx(rMode, gp, tStep);
    } else {
        auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

        if ( mode == _3dMat ) {
            answer = mat->give3dMaterialStiffnessMatrix(rMode, gp, tStep);
        } else if ( mode == _PlaneStress ) {
            answer = mat->givePlaneStressStiffMtrx(rMode, gp, tStep);
        } else if ( mode == _PlaneStrain ) {
            answer = mat->givePlaneStrainStiffMtrx(rMode, gp, tStep);
        } else if ( mode == _1dMat ) {
            answer = mat->give1dStressStiffMtrx(rMode, gp, tStep);
        } else {
            mat->giveStiffnessMatrix(answer, rMode, gp, tStep);
        }
    }
}


FloatMatrixF<3,3>
SimpleCrossSection :: give2dBeamStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    auto mat1d = mat->give1dStressStiffMtrx(rMode, gp, tStep);
    auto area = this->give(CS_Area, gp);
    auto Iy   = this->give(CS_InertiaMomentY, gp);
    auto shearAreaz = this->give(CS_ShearAreaZ, gp);

    FloatMatrixF<3,3> answer;
    answer.at(1, 1) = mat1d.at(1, 1) * area;
    answer.at(2, 2) = mat1d.at(1, 1) * Iy;
    answer.at(3, 3) = shearAreaz * mat->give('G', gp);
    return answer;
}


FloatMatrixF<6,6>
SimpleCrossSection :: give3dBeamStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    auto mat1d = mat->give1dStressStiffMtrx(rMode, gp, tStep);
    double E    = mat1d.at(1, 1);
    double G    = mat->give('G', gp);
    double area = this->give(CS_Area, gp);
    double Iy   = this->give(CS_InertiaMomentY, gp);
    double Iz   = this->give(CS_InertiaMomentZ, gp);
    double Ik   = this->give(CS_TorsionConstantX, gp); // St. Venant torsional constant

    //shearCoeff = this->give(CS_BeamShearCoeff);
    double shearAreay = this->give(CS_ShearAreaY, gp);
    double shearAreaz = this->give(CS_ShearAreaZ, gp);

    FloatMatrixF<6,6> answer;

    answer.at(1, 1) = E * area;
    ///@todo Do this by using the general 3d tangent matrix instead somehow!!!
    answer.at(2, 2) = shearAreay * G;
    answer.at(3, 3) = shearAreaz * G;
    //answer.at(2, 2) = shearCoeff * G * area;
    //answer.at(3, 3) = shearCoeff * G * area;
    answer.at(4, 4) = G * Ik;
    answer.at(5, 5) = E * Iy;
    answer.at(6, 6) = E * Iz;
    return answer;
}


FloatMatrixF<5,5>
SimpleCrossSection :: give2dPlateStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    auto mat2d = mat->givePlaneStressStiffMtrx(rMode, gp, tStep);
    double thickness = this->give(CS_Thickness, gp);
    double thickness3 = thickness * thickness * thickness;

    FloatMatrixF<5, 5> answer;

    for ( int i = 1; i <= 2; i++ ) {
        for ( int j = 1; j <= 2; j++ ) {
            answer.at(i, j) = mat2d.at(i, j) * thickness3 / 12.;
        }
    }

    answer.at(3, 3) = mat2d.at(3, 3) * thickness3 / 12.;
    answer.at(4, 4) = mat2d.at(3, 3) * thickness * ( 5. / 6. );
    answer.at(5, 5) = answer.at(4, 4);
    return answer;
}


FloatMatrixF<8,8>
SimpleCrossSection :: give3dShellStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );

    double thickness = this->give(CS_Thickness, gp);
    double thickness3 = thickness * thickness * thickness;

    auto mat2d = mat->givePlaneStressStiffMtrx(rMode, gp, tStep);

    FloatMatrixF<8,8> answer;

    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = mat2d.at(i, j) * thickness;
        }
    }
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i + 3, j + 3) = mat2d.at(i, j) * thickness3 / 12.0;
        }
    }

    answer.at(8, 8) = answer.at(7, 7) = mat2d.at(3, 3) * thickness * ( 5. / 6. );
    return answer;
}


FloatMatrixF<6,6>
SimpleCrossSection :: give3dDegeneratedShellStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    auto answer = mat->give3dMaterialStiffnessMatrix(rMode, gp, tStep);

    answer.at(1, 1) -= answer.at(1, 3) * answer.at(3, 1) / answer.at(3, 3);
    answer.at(2, 1) -= answer.at(2, 3) * answer.at(3, 1) / answer.at(3, 3);
    answer.at(1, 2) -= answer.at(1, 3) * answer.at(3, 2) / answer.at(3, 3);
    answer.at(2, 2) -= answer.at(2, 3) * answer.at(3, 2) / answer.at(3, 3);

    answer.at(3, 1) = 0.0;
    answer.at(3, 2) = 0.0;
    answer.at(3, 3) = 0.0;
    answer.at(2, 3) = 0.0;
    answer.at(1, 3) = 0.0;
    return answer;
}

FloatMatrixF<4,4>
SimpleCrossSection :: giveMembraneRotStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    auto d = mat->givePlaneStressStiffMtrx(ElasticStiffness, gp, tStep);
    auto ds = assemble<4,4>(d, {0, 1, 2}, {0, 1, 2});
    ds.at(4, 4) = this->give(CS_DrillingStiffness, gp);
    return ds;
}

FloatMatrixF<3,3>
SimpleCrossSection :: give2dPlateSubSoilStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    return mat->give2dPlateSubSoilStiffMtrx(ElasticStiffness, gp, tStep);
}


void
SimpleCrossSection :: initializeFrom(InputRecord &ir)
{
    CrossSection :: initializeFrom(ir);

    double value;

    double thick = 0.0;
    if ( ir.hasField(_IFT_SimpleCrossSection_thick) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, thick, _IFT_SimpleCrossSection_thick);
        propertyDictionary.add(CS_Thickness, thick);
    }

    double width = 0.0;
    if ( ir.hasField(_IFT_SimpleCrossSection_width) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, width, _IFT_SimpleCrossSection_width);
        propertyDictionary.add(CS_Width, width);
    }

    double area = 0.0;
    if ( ir.hasField(_IFT_SimpleCrossSection_area) ) {
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
    propertyDictionary.add(CS_TorsionConstantX, value);

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

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_relDrillStiffness);
    propertyDictionary.add(CS_RelDrillingStiffness, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_drillType);
    propertyDictionary.add(CS_DrillingType, value);

    this->materialNumber = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->materialNumber, _IFT_SimpleCrossSection_MaterialNumber);

    if ( ir.hasField(_IFT_SimpleCrossSection_directorx) ) {
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

    if ( this->propertyDictionary.includes(CS_TorsionConstantX) ) {
        input.setField(this->propertyDictionary.at(CS_TorsionConstantX), _IFT_SimpleCrossSection_ik);
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
SimpleCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode) const
{
    if ( this->giveMaterialNumber() ) {
        return this->domain->giveMaterial( this->giveMaterialNumber() )->isCharacteristicMtrxSymmetric(rMode);
    } else {
        return false; // Bet false...
    }
}


Material *
SimpleCrossSection :: giveMaterial(IntegrationPoint *ip) const
{
    if ( this->giveMaterialNumber() ) {
        return this->giveDomain()->giveMaterial( this->giveMaterialNumber() );
    } else {
        return ip->giveElement()->giveMaterial();
    }
}


double
SimpleCrossSection :: give(int aProperty, GaussPoint *gp) const
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
        answer = mat->giveFirstPKStressVector_3d(reducedvF, gp, tStep);
    } else if ( mode == _PlaneStrain ) {
        answer = mat->giveFirstPKStressVector_PlaneStrain(reducedvF, gp, tStep);
    } else if ( mode == _PlaneStress ) {
        answer = mat->giveFirstPKStressVector_PlaneStress(reducedvF, gp, tStep);
    } else if ( mode == _1dMat ) {
        answer = mat->giveFirstPKStressVector_1d(reducedvF, gp, tStep);
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
        answer = mat->give3dMaterialStiffnessMatrix_dPdF(rMode, gp, tStep);
    } else if ( mode == _PlaneStress ) {
        answer = mat->givePlaneStressStiffMtrx_dPdF(rMode, gp, tStep);
    } else if ( mode == _PlaneStrain ) {
        answer = mat->givePlaneStrainStiffMtrx_dPdF(rMode, gp, tStep);
    } else if ( mode == _1dMat ) {
        answer = mat->give1dStressStiffMtrx_dPdF(rMode, gp, tStep);
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
SimpleCrossSection :: giveTemperatureVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const
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

void
SimpleCrossSection :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralCrossSection :: saveContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.write(materialNumber) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(czMaterialNumber) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}

void
SimpleCrossSection :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralCrossSection :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.read(materialNumber) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(czMaterialNumber) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}

} // end namespace oofem
