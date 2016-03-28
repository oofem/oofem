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

#include "../sm/CrossSections/layeredcrosssection.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "material.h"
#include "floatarray.h"
#include "contextioerr.h"
#include "gaussintegrationrule.h"
#include "mathfem.h"
#include "classfactory.h"
#include "lobattoir.h"
#include "dynamicinputrecord.h"
#include "cltypes.h"

namespace oofem {
REGISTER_CrossSection(LayeredCrossSection);


void
LayeredCrossSection :: giveRealStress_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    if ( gp->giveIntegrationRule()->giveIntegrationDomain() == _Cube || gp->giveIntegrationRule()->giveIntegrationDomain() == _Wedge ) {
        // Determine which layer the gp belongs to. This code assumes that the gauss point are created consistently (through CrossSection::setupIntegrationPoints)
        int ngps = gp->giveIntegrationRule()->giveNumberOfIntegrationPoints();
        int gpnum = gp->giveNumber();
        int gpsperlayer = ngps / this->numberOfLayers;
        int layer = ( gpnum - 1 ) / gpsperlayer + 1;
        Material *layerMat = this->domain->giveMaterial( this->giveLayerMaterial(layer) );
        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);

            FloatArray rotStress;
            FloatArray rotStrain(6);
            rotStrain.at(1) = c * c * strain.at(1) - c *s *strain.at(6) + s *s *strain.at(2);
            rotStrain.at(2) = c * c * strain.at(2) + c *s *strain.at(6) + s *s *strain.at(1);
            rotStrain.at(3) = strain.at(3);
            rotStrain.at(4) = c * strain.at(4) + s *strain.at(5);
            rotStrain.at(5) = c * strain.at(5) - s *strain.at(4);
            rotStrain.at(6) = ( c * c - s * s ) * strain.at(6) + 2 * c * s * ( strain.at(1) - strain.at(2) );

            static_cast< StructuralMaterial * >(layerMat)->giveRealStressVector_3d(rotStress, gp, rotStrain, tStep);

            answer = {
                c *c * rotStress.at(1) + 2 * c * s * rotStress.at(6) + s * s * rotStress.at(2),
                c * c * rotStress.at(2) - 2 * c * s * rotStress.at(6) + s * s * rotStress.at(1),
                rotStress.at(3),
                c * rotStress.at(4) - s * rotStress.at(5),
                c * rotStress.at(5) + s * rotStress.at(4),
                ( c * c - s * s ) * rotStress.at(6) - c * s * ( rotStress.at(1) - rotStress.at(2) ),
            };
        } else {
            static_cast< StructuralMaterial * >(layerMat)->giveRealStressVector_3d(answer, gp, strain, tStep);
        }
    } else {
        OOFEM_ERROR("Only cubes and wedges are meaningful for layered cross-sections");
    }
}


void
LayeredCrossSection :: giveRealStress_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
LayeredCrossSection :: giveRealStress_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
LayeredCrossSection :: giveRealStress_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
LayeredCrossSection :: giveRealStress_Warping(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
LayeredCrossSection :: giveStiffnessMatrix_3d(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    if ( gp->giveIntegrationRule()->giveIntegrationDomain() == _Cube || gp->giveIntegrationRule()->giveIntegrationDomain() == _Wedge ) {
        // Determine which layer the gp belongs to. This code assumes that the gauss point are created consistently (through CrossSection::setupIntegrationPoints)
        int ngps = gp->giveIntegrationRule()->giveNumberOfIntegrationPoints();
        int gpnum = gp->giveNumber();
        int gpsperlayer = ngps / this->numberOfLayers;
        int layer = ( gpnum - 1 ) / gpsperlayer + 1;
        Material *layerMat = this->domain->giveMaterial( this->giveLayerMaterial(layer) );
        static_cast< StructuralMaterial * >(layerMat)->give3dMaterialStiffnessMatrix(answer, rMode, gp, tStep);

        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);

            FloatMatrix rotTangent = {
                {  c *c,    s *s, 0,  0,  0,    -c *s },
                {  s *s,    c *c, 0,  0,  0,     c *s },
                {    0,      0, 1,  0,  0,       0 },
                {    0,      0, 0,  c,  s,       0 },
                {    0,      0, 0, -s,  c,       0 },
                { 2 * c * s, -2 * c * s, 0,  0,  0, c * c - s * s }
            };
            answer.rotatedWith(rotTangent, 't');
        }
    } else {
        OOFEM_ERROR("Only cubes and wedges are meaningful for layered cross-sections");
    }
}


void
LayeredCrossSection :: giveStiffnessMatrix_PlaneStress(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
LayeredCrossSection :: giveStiffnessMatrix_PlaneStrain(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
LayeredCrossSection :: giveStiffnessMatrix_1d(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
LayeredCrossSection :: giveGeneralizedStress_Beam2d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    FloatArray layerStrain, reducedLayerStress;
    StructuralElement *element = static_cast< StructuralElement * >( gp->giveElement() );
    LayeredCrossSectionInterface *interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );

    answer.resize(3);
    answer.zero();

    // perform integration over layers
    bottom = this->give(CS_BottomZCoord, gp);
    top = this->give(CS_TopZCoord, gp);

    if ( interface == NULL ) {
        OOFEM_ERROR("element with no layer support encountered");
    }

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        GaussPoint *layerGp = this->giveSlaveGaussPoint(gp, layer - 1);
        StructuralMaterial *layerMat = static_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(layer) ) );

        // resolve current layer z-coordinate
        layerThick = this->layerThicks.at(layer);
        layerWidth = this->layerWidths.at(layer);
        layerZeta = layerGp->giveNaturalCoordinate(3);
        layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, gp, layerGp, tStep);

        if ( this->layerRots.at(layer) != 0. ) {
            OOFEM_ERROR("Rotation not supported for beams");
        } else {
            layerMat->giveRealStressVector_2dBeamLayer(reducedLayerStress, layerGp, layerStrain, tStep);
        }

        answer.at(1) += reducedLayerStress.at(1) * layerWidth * layerThick;
        answer.at(2) += reducedLayerStress.at(1) * layerWidth * layerThick * layerZCoord;
        answer.at(3) += reducedLayerStress.at(2) * layerWidth * layerThick;
    }

    // Create material status according to the first layer material
    ///@todo This should be replaced with a general "CrossSectionStatus"
    //CrossSectionStatus *status = new CrossSectionStatus(gp);
    //gp->setMaterialStatus(status);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >
                                       ( domain->giveMaterial( layerMaterials.at(1) )->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
LayeredCrossSection :: giveGeneralizedStress_Beam3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
LayeredCrossSection :: giveGeneralizedStress_Plate(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    FloatArray layerStrain, reducedLayerStress;
    StructuralElement *element = static_cast< StructuralElement * >( gp->giveElement() );
    LayeredCrossSectionInterface *interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );

    answer.resize(5);
    answer.zero();

    // perform integration over layers
    bottom = this->give(CS_BottomZCoord, gp);
    top = this->give(CS_TopZCoord, gp);

    if ( interface == NULL ) {
        OOFEM_ERROR("element with no layer support encountered");
    }

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        GaussPoint *layerGp = this->giveSlaveGaussPoint(gp, layer - 1);
        StructuralMaterial *layerMat = static_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(layer) ) );

        // resolve current layer z-coordinate
        layerThick = this->layerThicks.at(layer);
        layerWidth = this->layerWidths.at(layer);
        layerZeta = layerGp->giveNaturalCoordinate(3);
        layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, gp, layerGp, tStep);

        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);

            FloatArray rotStress;
            FloatArray rotStrain(5);
            rotStrain.at(1) = c * c * strain.at(1) - c *s *strain.at(5) + s *s *strain.at(2);
            rotStrain.at(2) = c * c * strain.at(2) + c *s *strain.at(5) + s *s *strain.at(1);
            rotStrain.at(3) = c * strain.at(3) + s *strain.at(4);
            rotStrain.at(4) = c * strain.at(4) - s *strain.at(3);
            rotStrain.at(5) = ( c * c - s * s ) * strain.at(5) + c * s * ( strain.at(1) - strain.at(2) );

            layerMat->giveRealStressVector_PlateLayer(rotStress, gp, rotStrain, tStep);

            answer = {
                c *c * rotStress.at(1) + 2 * c * s * rotStress.at(5) + s * s * rotStress.at(2),
                c * c * rotStress.at(2) - 2 * c * s * rotStress.at(5) + s * s * rotStress.at(1),
                c * rotStress.at(3) - s * rotStress.at(4),
                c * rotStress.at(4) + s * rotStress.at(3),
                ( c * c - s * s ) * rotStress.at(5) - c * s * ( rotStress.at(1) - rotStress.at(2) ),
            };
        } else {
            layerMat->giveRealStressVector_PlateLayer(reducedLayerStress, layerGp, layerStrain, tStep);
        }

        answer.at(1) += reducedLayerStress.at(1) * layerWidth * layerThick * layerZCoord;
        answer.at(2) += reducedLayerStress.at(2) * layerWidth * layerThick * layerZCoord;
        answer.at(3) += reducedLayerStress.at(5) * layerWidth * layerThick * layerZCoord;
        answer.at(4) += reducedLayerStress.at(4) * layerWidth * layerThick;
        answer.at(5) += reducedLayerStress.at(3) * layerWidth * layerThick;
    }

    // now we must update master gp
    // Create material status according to the first layer material
    ///@todo This should be replaced with a general "CrossSectionStatus"
    //CrossSectionStatus *status = new CrossSectionStatus(gp);
    //gp->setMaterialStatus(status);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >
                                       ( domain->giveMaterial( layerMaterials.at(1) )->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
LayeredCrossSection :: giveGeneralizedStress_Shell(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    FloatArray layerStrain, reducedLayerStress;
    StructuralElement *element = static_cast< StructuralElement * >( gp->giveElement() );
    LayeredCrossSectionInterface *interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );

    answer.resize(8);
    answer.zero();

    // perform integration over layers
    bottom = this->give(CS_BottomZCoord, gp);
    top = this->give(CS_TopZCoord, gp);

    if ( interface == NULL ) {
        OOFEM_ERROR("element with no layer support encountered");
    }

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        GaussPoint *layerGp = this->giveSlaveGaussPoint(gp, layer - 1);
        StructuralMaterial *layerMat = static_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(layer) ) );

        // resolve current layer z-coordinate
        layerThick = this->layerThicks.at(layer);
        layerWidth = this->layerWidths.at(layer);
        layerZeta = layerGp->giveNaturalCoordinate(3);
        layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, gp, layerGp, tStep);

        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);

            FloatArray rotStress;
            FloatArray rotStrain = {
                c *c * strain.at(1) - c * s * strain.at(5) + s * s * strain.at(2),
                c * c * strain.at(2) + c * s * strain.at(5) + s * s * strain.at(1),
                c * strain.at(3) + s * strain.at(4),
                c * strain.at(4) - s * strain.at(3),
                ( c * c - s * s ) * strain.at(5) + c * s * ( strain.at(1) - strain.at(2) ),
            };

            layerMat->giveRealStressVector_PlateLayer(rotStress, gp, rotStrain, tStep);

            answer = {
                c *c * rotStress.at(1) + 2 * c * s * rotStress.at(5) + s * s * rotStress.at(2),
                c * c * rotStress.at(2) - 2 * c * s * rotStress.at(5) + s * s * rotStress.at(1),
                c * rotStress.at(3) - s * rotStress.at(4),
                c * rotStress.at(4) + s * rotStress.at(3),
                ( c * c - s * s ) * rotStress.at(5) - c * s * ( rotStress.at(1) - rotStress.at(2) ),
            };
        } else {
            layerMat->giveRealStressVector_PlateLayer(reducedLayerStress, layerGp, layerStrain, tStep);
        }

        // 1) membrane terms sx, sy, sxy
        answer.at(1) += reducedLayerStress.at(1) * layerWidth * layerThick;
        answer.at(2) += reducedLayerStress.at(2) * layerWidth * layerThick;
        answer.at(3) += reducedLayerStress.at(5) * layerWidth * layerThick;
        // 2) bending terms mx, my, mxy
        answer.at(4) += reducedLayerStress.at(1) * layerWidth * layerThick * layerZCoord;
        answer.at(5) += reducedLayerStress.at(2) * layerWidth * layerThick * layerZCoord;
        answer.at(6) += reducedLayerStress.at(5) * layerWidth * layerThick * layerZCoord;
        // 3) shear terms qx, qy
        answer.at(7) += reducedLayerStress.at(4) * layerWidth * layerThick;
        answer.at(8) += reducedLayerStress.at(3) * layerWidth * layerThick;
    }

    // now we must update master gp
    ///@todo This should be replaced with a general "CrossSectionStatus"
    //CrossSectionStatus *status = new CrossSectionStatus(gp);
    //gp->setMaterialStatus(status);
    // Create material status according to the first layer material
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >
                                       ( domain->giveMaterial( layerMaterials.at(1) )->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
LayeredCrossSection :: giveGeneralizedStress_MembraneRot(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported in given cross-section (yet).");
}

void 
LayeredCrossSection :: giveGeneralizedStress_PlateSubSoil(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported in given cross-section (yet).");
}

void
LayeredCrossSection :: giveCharMaterialStiffnessMatrix(FloatMatrix &answer,
                                                       MatResponseMode rMode,
                                                       GaussPoint *gp,
                                                       TimeStep *tStep)
//
// only interface to material class, forcing returned matrix to be in reduced form.
//
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
        int ngps = gp->giveIntegrationRule()->giveNumberOfIntegrationPoints();
        int gpnum = gp->giveNumber();
        int gpsperlayer = ngps / this->numberOfLayers;
        int layer = ( gpnum - 1 ) / gpsperlayer + 1;
        StructuralMaterial *mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->giveLayerMaterial(layer) ) );
        if ( mat->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
            mat->giveStiffnessMatrix(answer, rMode, gp, tStep);
        } else {
            OOFEM_ERROR("unsupported StressStrainMode");
        }
    }
}


void
LayeredCrossSection :: give2dPlateStiffMtrx(FloatMatrix &answer,
                                            MatResponseMode rMode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)

//
// assumption sigma_z = 0.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    FloatMatrix layerMatrix;
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    double layerZCoord2;

    answer.resize(5, 5);
    answer.zero();

    // perform integration over layers
    bottom = this->give(CS_BottomZCoord, gp);
    top = this->give(CS_TopZCoord, gp);

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        GaussPoint *layerGp = giveSlaveGaussPoint(gp, layer - 1);

        ///@todo Just using the gp number doesn't nicely support more than 1 gp per layer. Must rethink.
        StructuralMaterial *mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->giveLayerMaterial(layer) ) );
        mat->givePlateLayerStiffMtrx(layerMatrix, rMode, layerGp, tStep);
        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);
            FloatMatrix rotTangent = {
                {  c *c,    s *s,  0,  0,    -c *s },
                {  s *s,    c *c,  0,  0,     c *s },
                {    0,      0,  c,  s,       0 },
                {    0,      0, -s,  c,       0 },
                { 2 * c * s, -2 * c * s,  0,  0, c * c - s * s }
            };
            layerMatrix.rotatedWith(rotTangent, 't');
        }

        //
        // resolve current layer z-coordinate
        //
        layerThick = this->layerThicks.at(layer);
        layerWidth  = this->layerWidths.at(layer);
        layerZeta   = layerGp->giveNaturalCoordinate(3);
        layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );
        layerZCoord2 = layerZCoord * layerZCoord;
        //
        // perform integration
        //
        // 1) bending terms mx, my, mxy

        answer.at(1, 1) += layerMatrix.at(1, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(1, 2) += layerMatrix.at(1, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(1, 3) += layerMatrix.at(1, 5) * layerWidth * layerThick * layerZCoord2;

        answer.at(2, 1) += layerMatrix.at(2, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(2, 2) += layerMatrix.at(2, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(2, 3) += layerMatrix.at(2, 5) * layerWidth * layerThick * layerZCoord2;

        answer.at(3, 1) += layerMatrix.at(5, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(3, 2) += layerMatrix.at(5, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(3, 3) += layerMatrix.at(5, 5) * layerWidth * layerThick * layerZCoord2;

        // 2) shear terms qx = qxz, qy = qyz
        answer.at(4, 4) += layerMatrix.at(4, 4) * layerWidth * layerThick;
        answer.at(4, 5) += layerMatrix.at(4, 3) * layerWidth * layerThick;
        answer.at(5, 4) += layerMatrix.at(3, 4) * layerWidth * layerThick;
        answer.at(5, 5) += layerMatrix.at(3, 3) * layerWidth * layerThick;
    }
}


void
LayeredCrossSection :: give3dShellStiffMtrx(FloatMatrix &answer,
                                            MatResponseMode rMode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
//
// assumption sigma_z = 0.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    FloatMatrix layerMatrix;
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    // double zi, zi1;
    double layerZCoord2;

    answer.resize(8, 8);
    answer.zero();
    // perform integration over layers
    bottom = this->give(CS_BottomZCoord, gp);
    top = this->give(CS_TopZCoord, gp);

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        GaussPoint *layerGp = giveSlaveGaussPoint(gp, layer - 1);

        ///@todo The logic in this whole class is pretty messy to support both slave-gp's and normal gps. Rethinking the approach is necessary.
        /// Just using the gp number doesn't nicely support more than 1 gp per layer. Must rethink.
        StructuralMaterial *mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->giveLayerMaterial(layer) ) );
        mat->givePlateLayerStiffMtrx(layerMatrix, rMode, layerGp, tStep);
        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot);
            double s = sin(rot);
            FloatMatrix rotTangent = {
                {  c *c,    s *s,  0,  0,    -c *s },
                {  s *s,    c *c,  0,  0,     c *s },
                {    0,      0,  c,  s,       0 },
                {    0,      0, -s,  c,       0 },
                { 2 * c * s, -2 * c * s,  0,  0, c * c - s * s }
            };
            layerMatrix.rotatedWith(rotTangent, 't');
        }

        //
        // resolve current layer z-coordinate
        //
        layerThick = this->layerThicks.at(layer);
        layerWidth  = this->layerWidths.at(layer);
        layerZeta   = layerGp->giveNaturalCoordinate(3);
        layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );
        layerZCoord2 = layerZCoord * layerZCoord;
        //
        // perform integration
        //
        // 1) membrane terms sx, sy, sxy
        answer.at(1, 1) += layerMatrix.at(1, 1) * layerWidth * layerThick;
        answer.at(1, 2) += layerMatrix.at(1, 2) * layerWidth * layerThick;
        answer.at(1, 3) += layerMatrix.at(1, 5) * layerWidth * layerThick;

        answer.at(2, 1) += layerMatrix.at(2, 1) * layerWidth * layerThick;
        answer.at(2, 2) += layerMatrix.at(2, 2) * layerWidth * layerThick;
        answer.at(2, 3) += layerMatrix.at(2, 5) * layerWidth * layerThick;

        answer.at(3, 1) += layerMatrix.at(5, 1) * layerWidth * layerThick;
        answer.at(3, 2) += layerMatrix.at(5, 2) * layerWidth * layerThick;
        answer.at(3, 3) += layerMatrix.at(5, 5) * layerWidth * layerThick;

        // 2) bending terms mx, my, mxy

        answer.at(4, 4) += layerMatrix.at(1, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(4, 5) += layerMatrix.at(1, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(4, 6) += layerMatrix.at(1, 5) * layerWidth * layerThick * layerZCoord2;

        answer.at(5, 4) += layerMatrix.at(2, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(5, 5) += layerMatrix.at(2, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(5, 6) += layerMatrix.at(2, 5) * layerWidth * layerThick * layerZCoord2;

        answer.at(6, 4) += layerMatrix.at(5, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(6, 5) += layerMatrix.at(5, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(6, 6) += layerMatrix.at(5, 5) * layerWidth * layerThick * layerZCoord2;

        // 3) shear terms qx, qy
        answer.at(7, 7) += layerMatrix.at(4, 4) * layerWidth * layerThick;
        answer.at(7, 8) += layerMatrix.at(4, 3) * layerWidth * layerThick;
        answer.at(8, 7) += layerMatrix.at(3, 4) * layerWidth * layerThick;
        answer.at(8, 8) += layerMatrix.at(3, 3) * layerWidth * layerThick;
    }
}


void
LayeredCrossSection :: give2dBeamStiffMtrx(FloatMatrix &answer,
                                           MatResponseMode rMode,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
//
// assumption sigma_z = 0.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    FloatMatrix layerMatrix;
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    double layerZCoord2;

    // perform integration over layers
    bottom = this->give(CS_BottomZCoord, gp);
    top = this->give(CS_TopZCoord, gp);

    answer.resize(3, 3);
    answer.zero();

    for ( int i = 1; i <= numberOfLayers; i++ ) {
        GaussPoint *layerGp = giveSlaveGaussPoint(gp, i - 1);

        ///@todo The logic in this whole class is pretty messy to support both slave-gp's and normal gps. Rethinking the approach is necessary.
        /// Just using the gp number doesn't nicely support more than 1 gp per layer. Must rethink.
        StructuralMaterial *mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->giveLayerMaterial(i) ) );
        mat->give2dBeamLayerStiffMtrx(layerMatrix, rMode, layerGp, tStep);
        if ( this->layerRots.at(i) != 0. ) {
            OOFEM_ERROR("Doesn't support layer rotations.");
        }

        //
        // resolve current layer z-coordinate
        //
        layerThick = this->layerThicks.at(i);
        layerWidth  = this->layerWidths.at(i);
        layerZeta   = layerGp->giveNaturalCoordinate(3);
        layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );
        layerZCoord2 = layerZCoord * layerZCoord;
        //
        // perform integration
        //
        // 1) membrane terms sx
        answer.at(1, 1) += layerMatrix.at(1, 1) * layerWidth * layerThick;
        answer.at(1, 3) += layerMatrix.at(1, 2) * layerWidth * layerThick;
        // 2) bending terms my
        answer.at(2, 2) += layerMatrix.at(1, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(2, 3) += layerMatrix.at(1, 2) * layerWidth * layerThick * layerZCoord2;
        // 3) shear terms qx
        answer.at(3, 1) += layerMatrix.at(2, 1) * layerWidth * layerThick;
        answer.at(3, 3) += layerMatrix.at(2, 2) * layerWidth * layerThick;
    }
}


void
LayeredCrossSection :: give3dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
}


void
LayeredCrossSection :: giveMembraneRotStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
}

void
LayeredCrossSection :: give2dPlateSubSoilStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
}



FloatArray *
LayeredCrossSection :: imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStressVector3d)
//
// returns modified gradient of stress vector, which is used to
// bring stresses back to yield surface.
//
// imposes zeros on places, where zero stress occurs. if energetically connected
// strain is zero, we do not impose zero there, because stress exist and
// must be taken into account when computing yeld function. In such case
// a problem is assumed to be full 3d with some explicit strain equal to 0.
//
// On the other hand, if some stress is imposed to be zero, we understand
// such case as subspace of 3d case (like a classical plane stess problem, with no
// tracing of ez, sigma_z)
//
{
    MaterialMode mode = gp->giveMaterialMode();
    int size = gradientStressVector3d->giveSize();
    if ( size != 6 ) {
        OOFEM_ERROR("size mismatch");
    }

    switch ( mode ) {
    case _PlateLayer:
        gradientStressVector3d->at(3) = 0.;
        break;
    case _2dBeamLayer:
        for ( int i = 2; i <= 5; i++ ) {
            gradientStressVector3d->at(i) = 0.;
        }

        break;
    default:
        StructuralCrossSection :: imposeStressConstrainsOnGradient(gp, gradientStressVector3d);
    }

    return gradientStressVector3d;
}


FloatArray *
LayeredCrossSection :: imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStrainVector3d)
//
// returns modified gradient of strain vector, which is used to
// compute plastic strain increment.
//
// imposes zeros on places, where zero strain occurs or energetically connected stress
// is prescribed to be zero.
//
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( gradientStrainVector3d->giveSize() != 6 ) {
        OOFEM_ERROR("size mismatch");
    }

    switch ( mode ) {
    case _PlateLayer:
        gradientStrainVector3d->at(3) = 0.;
        break;
    case _2dBeamLayer:
        for ( int i = 2; i <= 5; i++ ) {
            gradientStrainVector3d->at(i) = 0.;
        }

        break;
    default:
        StructuralCrossSection :: imposeStrainConstrainsOnGradient(gp, gradientStrainVector3d);
    }

    return gradientStrainVector3d;
}


IRResultType
LayeredCrossSection :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = CrossSection :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    IR_GIVE_FIELD(ir, numberOfLayers, _IFT_LayeredCrossSection_nlayers);
    IR_GIVE_FIELD(ir, layerMaterials, _IFT_LayeredCrossSection_layermaterials);
    IR_GIVE_FIELD(ir, layerThicks, _IFT_LayeredCrossSection_thicks);
    layerWidths.resize(numberOfLayers);
    layerWidths.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, layerWidths, _IFT_LayeredCrossSection_widths);
    layerRots.resize(numberOfLayers);
    layerRots.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, layerRots, _IFT_LayeredCrossSection_layerRotations);

    if ( numberOfLayers != layerMaterials.giveSize() ||
        numberOfLayers != layerThicks.giveSize()  ||
        numberOfLayers != layerRots.giveSize() ) {  //|| ( numberOfLayers != layerWidths.giveSize() ) )
        OOFEM_WARNING("numberOfLayers does not equal given number of thicknesses. ");
        return IRRT_BAD_FORMAT;
    }

    if ( numberOfLayers <= 0 ) {
        OOFEM_WARNING("numberOfLayers<= 0 is not allowed");
        return IRRT_BAD_FORMAT;
    }

    // Interface materials // add check if correct numbers
    interfacerMaterials.resize(numberOfLayers - 1);
    interfacerMaterials.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, interfacerMaterials, _IFT_LayeredCrossSection_interfacematerials);

    numberOfIntegrationPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfIntegrationPoints, _IFT_LayeredCrossSection_nintegrationpoints);

    // read z-coordinate of mid-surface measured from bottom layer
    midSurfaceZcoordFromBottom = 0.5 * this->computeIntegralThick();  // Default: geometric midplane
    midSurfaceXiCoordFromBottom = 1.0; // add to IR
    IR_GIVE_OPTIONAL_FIELD(ir, midSurfaceZcoordFromBottom, _IFT_LayeredCrossSection_midsurf);

    this->setupLayerMidPlanes();

    return IRRT_OK;
}

void LayeredCrossSection :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralCrossSection :: giveInputRecord(input);

    input.setField(this->numberOfLayers, _IFT_LayeredCrossSection_nlayers);
    input.setField(this->layerMaterials, _IFT_LayeredCrossSection_layermaterials);
    input.setField(this->layerThicks, _IFT_LayeredCrossSection_thicks);
    input.setField(this->layerWidths, _IFT_LayeredCrossSection_widths);
    input.setField(this->layerRots, _IFT_LayeredCrossSection_layerRotations);
    input.setField(this->interfacerMaterials, _IFT_LayeredCrossSection_interfacematerials);
    input.setField(this->numberOfIntegrationPoints, _IFT_LayeredCrossSection_nintegrationpoints);
    input.setField(this->midSurfaceZcoordFromBottom, _IFT_LayeredCrossSection_midsurf);
}

void LayeredCrossSection :: createMaterialStatus(GaussPoint &iGP)
{
    for ( int i = 1; i <= numberOfLayers; i++ ) {
        GaussPoint *layerGp = giveSlaveGaussPoint(& iGP, i - 1);
        StructuralMaterial *mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->giveLayerMaterial(i) ) );
        MaterialStatus *matStat = mat->CreateStatus(layerGp);
        layerGp->setMaterialStatus(matStat);
    }
}


void
LayeredCrossSection :: setupLayerMidPlanes()
{
    // z-coord of each layer midplane measured from the global cross section z-coord
    this->layerMidZ.resize(this->numberOfLayers);
    double layerBottomZ = -midSurfaceZcoordFromBottom; // initialize to the bottom coord
    for ( int j = 1; j <= numberOfLayers; j++ ) {
        double thickness = this->layerThicks.at(j);
        this->layerMidZ.at(j) = layerBottomZ + thickness * 0.5;
        layerBottomZ += thickness;
    }
}


Material *
LayeredCrossSection :: giveMaterial(IntegrationPoint *ip)
{
    return this->domain->giveMaterial( this->giveLayerMaterial(ip->giveNumber()) );

}


int
LayeredCrossSection :: setupIntegrationPoints(IntegrationRule &irule, int npoints, Element *element)
{
    ///@todo We must send arrays for integration points instead of just a single scalar.
    if ( element->giveIntegrationDomain() == _Cube ) {
#if 0
        ///@todo "npoints" should be an intarray
        return irule.SetUpPointsOnCubeLayers(npoints.at(1), npoints.at(2), this->numberOfIntegrationPoints,
                                             element->giveMaterialMode(), this->layerThicks);

#else
        int points1 = ( int ) floor(cbrt( double ( npoints ) ) + 0.5);
        // If numberOfIntegrationPoints > 0 then use that, otherwise use the element's default.
        return irule.SetUpPointsOnCubeLayers(points1, points1, this->numberOfIntegrationPoints ? numberOfIntegrationPoints : points1,
                                             element->giveMaterialMode(), this->layerThicks);

#endif
    } else if ( element->giveIntegrationDomain() == _Wedge ) {
#if 0
        ///@todo "npoints" should be an intarray
        return irule.SetUpPointsOnWedgeLayers(npoints.at(1), this->numberOfIntegrationPoints,
                                              element->giveMaterialMode(), this->layerThicks);

#else
        if ( npoints == 2 ) {
            return irule.SetUpPointsOnWedgeLayers(1, this->numberOfIntegrationPoints,
                                                  element->giveMaterialMode(), this->layerThicks);
        } else {
            return irule.SetUpPointsOnWedgeLayers(3, this->numberOfIntegrationPoints,
                                                  element->giveMaterialMode(), this->layerThicks);
        }
#endif
    } else {
        return irule.setUpIntegrationPoints( element->giveIntegrationDomain(), npoints, element->giveMaterialMode() );
    }
}

GaussPoint *
LayeredCrossSection :: giveSlaveGaussPoint(GaussPoint *masterGp, int i)
//
// return the i-th slave gauss point of master gp
// if slave gp don't exists - create them
//
{
    GaussPoint *slave = masterGp->giveSlaveGaussPoint(i);
    if ( slave == NULL ) {
        // check for proper dimensions - slave can be NULL if index too high or if not
        // slaves previously defined
        if ( i > this->numberOfLayers ) {
            OOFEM_ERROR("no such layer defined");
        }

        // create new slave record in masterGp
        // (requires that this is friend of gp)
        double currentZTopCoord, currentZCoord,  bottom, top;
        const FloatArray &masterCoords = masterGp->giveNaturalCoordinates();
        // resolve slave material mode
        MaterialMode slaveMode, masterMode = masterGp->giveMaterialMode();
        slaveMode = this->giveCorrespondingSlaveMaterialMode(masterMode);

        bottom = this->give(CS_BottomZCoord, masterGp);
        top = this->give(CS_TopZCoord, masterGp);

        ///@todo Generalize to multiple integration points per layer
        masterGp->gaussPoints.resize( numberOfLayers );
        currentZTopCoord = -midSurfaceZcoordFromBottom;
        for ( int j = 0; j < numberOfLayers; j++ ) {
            FloatArray zCoord(3);
            currentZTopCoord += this->layerThicks.at(j + 1);
            currentZCoord = currentZTopCoord - this->layerThicks.at(j + 1) / 2.0; // z-coord of layer mid surface
            if ( masterCoords.giveSize() > 0 ) {
                zCoord.at(1) = masterCoords.at(1); // gp x-coord of mid surface
            }

            if ( masterCoords.giveSize() > 1 ) {
                zCoord.at(2) = masterCoords.at(2); // gp y-coord of mid surface
            }

            zCoord.at(3) = ( 2.0 * currentZCoord - top - bottom ) / ( top - bottom );
            // in gp - is stored isoparametric coordinate (-1,1) of z-coordinate
            //masterGp->gaussPoints [ j ] = new GaussPoint(masterGp->giveIntegrationRule(), j + 1, zCoord, 0., slaveMode);

            // test - remove!
            masterGp->gaussPoints [ j ] = new GaussPoint(masterGp->giveIntegrationRule(), j + 1, zCoord, 1.0, slaveMode);
        }

        slave = masterGp->gaussPoints [ i ];
    }

    return slave;
}






double
LayeredCrossSection :: computeIntegralThick()
//
// computes total thickness of receiver
//
{
    if ( totalThick == 0 ) {
        totalThick = layerThicks.sum();
    }

    return totalThick;
}


void
LayeredCrossSection :: printYourself()
// Prints the receiver on screen.
{
    printf("Cross Section with properties: \n");
    propertyDictionary.printYourself();
    printf("Layer Materials: \n");
    layerMaterials.printYourself();
    printf("Thickness of each layer: \n");
    layerThicks.printYourself();
    if ( layerWidths.giveSize() ) {
        printf("Width of each layer: \n");
        layerWidths.printYourself();
    }
    printf("Number of integration points per layer: %i \n", this->numberOfIntegrationPoints);
    printf("MidSurfaceZCoordinate from bottom: %f \n", midSurfaceZcoordFromBottom);
}


contextIOResultType
LayeredCrossSection :: saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *masterGp)
//
// saves full material context (saves state variables, that completely describe
// current state)
// stores also slaves records of master gp
//
{
    contextIOResultType iores;

    if ( ( iores = CrossSection :: saveIPContext(stream, mode, masterGp) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // saved master gp record;

    // and now save slave gp of master:
    for ( int i = 1; i <= numberOfLayers; i++ ) {
        GaussPoint *slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1);
        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(i) ) );
        if ( ( iores = mat->saveIPContext(stream, mode, slaveGP) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType
LayeredCrossSection :: restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *masterGp)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
// restores also slaves of master gp
//
{
    contextIOResultType iores;

    if ( ( iores = CrossSection :: restoreIPContext(stream, mode, masterGp) ) != CIO_OK ) {
        THROW_CIOERR(iores);                                                                    // saved masterGp
    }

    // and now save slave gp of master:
    for ( int i = 1; i <= numberOfLayers; i++ ) {
        // creates also slaves if they don't exists
        GaussPoint *slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1);
        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(i) ) );
        if ( ( iores = mat->restoreIPContext(stream, mode, slaveGP) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


MaterialMode
LayeredCrossSection :: giveCorrespondingSlaveMaterialMode(MaterialMode masterMode)
//
// returns corresponding slave material mode to master mode
//
{
    if ( masterMode == _2dPlate ) {
        return _PlateLayer;
    } else if ( masterMode == _2dBeam ) {
        return _2dBeamLayer;
    } else if ( masterMode == _3dShell ) {
        return _PlateLayer;
    } else if ( masterMode == _3dMat ) {
        return _3dMat;
    } else {
        OOFEM_ERROR("unsupported mode");
    }

    return _Unknown;
}


double
LayeredCrossSection :: give(CrossSectionProperty aProperty, GaussPoint *gp)
{
    if ( aProperty == CS_Thickness ) {
        return this->computeIntegralThick();
    } else if ( aProperty == CS_TopZCoord ) {
        this->computeIntegralThick();
        return totalThick - midSurfaceZcoordFromBottom;
    } else if ( aProperty == CS_BottomZCoord ) {
        return -midSurfaceZcoordFromBottom;
    } else if ( aProperty == CS_Area ) {
        return this->giveArea();
    } else if ( aProperty == CS_NumLayers ) {
        return this->numberOfLayers;
    }

    return CrossSection :: give(aProperty, gp);
}
double
LayeredCrossSection :: give(CrossSectionProperty aProperty, const FloatArray &coords, Element *elem, bool local)
{
    if ( aProperty == CS_Thickness ) {
        return this->computeIntegralThick();
    } else if ( aProperty == CS_TopZCoord ) {
        this->computeIntegralThick();
        return totalThick - midSurfaceZcoordFromBottom;
    } else if ( aProperty == CS_BottomZCoord ) {
        return -midSurfaceZcoordFromBottom;
    } else if ( aProperty == CS_Area ) {
        return this->giveArea();
    } else if ( aProperty == CS_NumLayers ) {
        return this->numberOfLayers;
    }

    return CrossSection :: give(aProperty, coords, elem, local);
}


int
LayeredCrossSection :: giveNumberOfLayers()
{
    return this->numberOfLayers;
}

double
LayeredCrossSection :: giveArea()
{
    if ( this->area <= 0.0 ) {
        this->area = this->layerThicks.dotProduct(this->layerWidths);
    }

    return area;
}


bool LayeredCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode)
{
    for ( int i = 1; i <= this->numberOfLayers; i++ ) {
        if ( !this->domain->giveMaterial( this->giveLayerMaterial(i) )->isCharacteristicMtrxSymmetric(rMode) ) {
            return false;
        }
    }
    return true;
}


void
LayeredCrossSection :: giveInterfaceXiCoords(FloatArray &answer)
{
    // returns an array with the xi-coords corresponding to the boundaries where
    // the layers meet (size = number of layers -1)

    int numInterfaces = this->giveNumberOfLayers() - 1;
    answer.resize(numInterfaces);
    double totalThickness = this->computeIntegralThick();
    for ( int i = 1; i <= numInterfaces; i++  ) {
        double midZ = this->giveLayerMidZ(i);
        double interfaceZ  = midZ + this->giveLayerThickness(i) * 0.5;
        answer.at(i) = interfaceZ * ( 2.0 / totalThickness );
    }
}

void
LayeredCrossSection :: setupLayeredIntegrationRule(std :: vector< std :: unique_ptr< IntegrationRule > > &integrationRulesArray, Element *el, int numInPlanePoints)
{
    // Loop over each layer and set up an integration rule as if each layer was an independent element
    // @todo - only works for wedge integration at the moment
    int numberOfLayers     = this->giveNumberOfLayers();
    int numPointsThickness = this->giveNumIntegrationPointsInLayer();

    integrationRulesArray.clear();
    integrationRulesArray.reserve( numberOfLayers );
    for ( int i = 0; i < numberOfLayers; i++ ) {
        integrationRulesArray.emplace_back( new LayeredIntegrationRule(i + 1, el) );
        integrationRulesArray.back()->SetUpPointsOnWedge(numInPlanePoints, numPointsThickness, _3dMat);
    }
    this->mapLayerGpCoordsToShellCoords(integrationRulesArray);
}


void
LayeredCrossSection :: mapLayerGpCoordsToShellCoords(std :: vector< std :: unique_ptr< IntegrationRule > > &layerIntegrationRulesArray)
/*
 * Maps the local xi-coord (z-coord) in each layer [-1,1] to the corresponding
 * xi-coord in the cross section coordinate system.
 * Also renames the gp numbering from layerwise to global (1,2,1,2 -> 1,2,3,4)
 *  xi
 * --------  1               --------  1
 |           |                         |
 |           |                         |
 | -------- -1       =>      --------  x
 | --------  1               --------  x
 |           |                         |
 | -------- -1               -------- -1
 */
{
    double scaleFactor = 0.999; // Will be numerically unstable with xfem if the endpoints lie at +-1
    double totalThickness = this->computeIntegralThick();
    int number = 1;
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        for ( GaussPoint *gp: *layerIntegrationRulesArray [ layer - 1 ] ) {

            // Map local layer cs to local shell cs
            double zMid_i = this->giveLayerMidZ(layer); // global z-coord
            double xiMid_i = 1.0 - 2.0 * ( totalThickness - this->midSurfaceZcoordFromBottom - zMid_i ) / totalThickness; // local z-coord
            double deltaxi = gp->giveNaturalCoordinates().at(3) * this->giveLayerThickness(layer) / totalThickness; // distance from layer mid
            double xinew = xiMid_i + deltaxi * scaleFactor;
            FloatArray lcoords = gp->giveNaturalCoordinates();
            lcoords.at(3) = xinew;
            gp->setNaturalCoordinates(lcoords);
            gp->number = number;   // fix gp ordering
            number++;
        }
    }
}


LayeredIntegrationRule :: LayeredIntegrationRule(int n, Element *e,
                                                 int startIndx, int endIndx, bool dynamic) :
    IntegrationRule(n, e, startIndx, endIndx, dynamic) { }

LayeredIntegrationRule :: LayeredIntegrationRule(int n, Element *e) :
    IntegrationRule(n, e) { }

LayeredIntegrationRule :: ~LayeredIntegrationRule()
{ }

int
LayeredIntegrationRule :: SetUpPointsOnWedge(int nPointsTri, int nPointsThickness, MaterialMode mode)
{
    // Set up integration rule for a specific layer

    int nPoints = nPointsTri * nPointsThickness;
    //@todo - is not really a Gauss point but rather a hybrid.
    this->gaussPoints.resize( nPoints );

    // uses Gauss integration in the plane and Lobatto in the thickness
    FloatArray coords_xi1, coords_xi2, coords_xi, weights_tri, weights_thickness;
    GaussIntegrationRule   :: giveTriCoordsAndWeights(nPointsTri, coords_xi1, coords_xi2, weights_tri);
    //LobattoIntegrationRule :: giveLineCoordsAndWeights(nPointsThickness, coords_xi, weights_thickness );
    GaussIntegrationRule :: giveLineCoordsAndWeights(nPointsThickness, coords_xi, weights_thickness);

    // Assumes that the integration rules of the layers are the same such that the ordering of the ip's are also
    // the same =>  upperInterfacePoints.at(i) of one layer is paired with lowerInterfacePoints.at(i) of the next.
    // This will be used to estimate interlaminar stresses, sice values in the two ip will generally be different
    // due to beam/plate/shell theory assumptions.
    if ( nPointsThickness != 1 ) { // otherwise there are no points on the interface
        this->lowerInterfacePoints.resize(nPointsTri);
        this->upperInterfacePoints.resize(nPointsTri);
    }
    for ( int i = 1, ind = 0; i <= nPointsThickness; i++ ) {
        for ( int j = 1; j <= nPointsTri; j++ ) {
            this->gaussPoints [ ind ] =
                new GaussPoint(this, 1, {coords_xi1.at(j), coords_xi2.at(j), coords_xi.at(i)},
                               weights_tri.at ( j ) *weights_thickness.at ( i ), mode);

            // store interface points
            if ( i == 1 && nPointsThickness > 1 ) { //then lower surface
                this->lowerInterfacePoints.at(j) = ind;
            } else if ( i == nPointsThickness && nPointsThickness > 1 ) {  //then upper surface
                this->upperInterfacePoints.at(j) = ind;
            }
            ind++;
        }
    }
    return this->giveNumberOfIntegrationPoints();
}


int
LayeredCrossSection :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    for ( int i = 1; this->giveNumberOfLayers(); i++ ) {
        Material *mat = this->giveDomain()->giveMaterial( this->giveLayerMaterial(i) );
        if ( !dynamic_cast< StructuralMaterial * >(mat) ) {
            OOFEM_WARNING("material %s without structural support", mat->giveClassName() );
            result = 0;
            continue;
        }
    }
    return result;
}


int
LayeredCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( gp->giveIntegrationRule()->giveIntegrationDomain() == _Cube || gp->giveIntegrationRule()->giveIntegrationDomain() == _Wedge ) {
        // Determine which layer the gp belongs to. This code assumes that the gauss point are created consistently (through CrossSection::setupIntegrationPoints)
        int ngps = gp->giveIntegrationRule()->giveNumberOfIntegrationPoints();
        int gpnum = gp->giveNumber();
        int gpsperlayer = ngps / this->numberOfLayers;
        int layer = ( gpnum - 1 ) / gpsperlayer + 1;
        Material *layerMat = this->domain->giveMaterial( this->giveLayerMaterial(layer) );
        if ( this->layerRots.at(layer) != 0. ) {
            FloatArray rotVal; // the requested value in the material c.s.
            InternalStateValueType valType = giveInternalStateValueType(type);
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);

            int ret = layerMat->giveIPValue(rotVal, gp, type, tStep);
            if ( ret == 0 ) {
                return 0;
            }

            // Determine how to rotate it according to the value type
            if ( valType == ISVT_TENSOR_S3 ) {
                answer = {
                    c *c * rotVal.at(1) + 2 * c * s * rotVal.at(6) + s * s * rotVal.at(2),
                    c * c * rotVal.at(2) - 2 * c * s * rotVal.at(6) + s * s * rotVal.at(1),
                    rotVal.at(3),
                    c * rotVal.at(4) - s * rotVal.at(5),
                    c * rotVal.at(5) + s * rotVal.at(4),
                    ( c * c - s * s ) * rotVal.at(6) - c * s * ( rotVal.at(1) - rotVal.at(2) ),
                };
            } else if ( valType == ISVT_TENSOR_S3E ) {
                answer = {
                    c *c * rotVal.at(1) + c * s * rotVal.at(6) + s * s * rotVal.at(2),
                    c * c * rotVal.at(2) - c * s * rotVal.at(6) + s * s * rotVal.at(1),
                    rotVal.at(3),
                    c * rotVal.at(4) - s * rotVal.at(5),
                    c * rotVal.at(5) + s * rotVal.at(4),
                    ( c * c - s * s ) * rotVal.at(6) - 2 * c * s * ( rotVal.at(1) - rotVal.at(2) ),
                };
            } else if ( valType == ISVT_VECTOR ) {
                answer = {
                    c *rotVal.at(1) - s * rotVal.at(2), s * rotVal.at(1), +c * rotVal.at(2), rotVal.at(3)
                };
            } else if ( valType == ISVT_SCALAR ) {
                answer = rotVal;
            } else {
                return 0;
            }
            return 1;
        } else {
            return layerMat->giveIPValue(answer, gp, type, tStep);
        }
    } else {
        //return CrossSection :: giveIPValue(answer, gp, type, tStep);

        ///@todo so far this only works for el where each layer has its own integration rule
        int layer = gp->giveIntegrationRule()->giveNumber();
        return this->giveDomain()->giveMaterial( this->giveLayerMaterial(layer) )->giveIPValue(answer, gp, type, tStep);
    }
}


double
LayeredCrossSection :: give(int aProperty, GaussPoint* gp)
{
    double average = 0.;
    for ( int layer = 1; layer <= numberOfLayers; ++layer ) {
        Material *mat = this->giveDomain()->giveMaterial( giveLayerMaterial(layer) );
        average += mat->give(aProperty, gp) * giveLayerThickness(layer);
    }
    return average / this->totalThick;
}

} // end namespace oofem
