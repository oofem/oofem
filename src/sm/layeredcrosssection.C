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

#include "layeredcrosssection.h"
#include "structuralelement.h"
#include "gausspoint.h"
#include "material.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "floatarray.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {

REGISTER_CrossSection( LayeredCrossSection );


void
LayeredCrossSection :: giveRealStress_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    if ( gp->giveIntegrationRule()->giveIntegrationDomain() == _Cube || gp->giveIntegrationRule()->giveIntegrationDomain() == _Wedge ) {
        // Determine which layer the gp belongs to. This code assumes that the gauss point are created consistently (through CrossSection::setupIntegrationPoints)
        int ngps = gp->giveIntegrationRule()->giveNumberOfIntegrationPoints();
        int gpnum = gp->giveNumber();
        int gpsperlayer = ngps / this->numberOfLayers;
        int layer = (gpnum - 1) / gpsperlayer + 1;
        Material *layerMat = this->domain->giveMaterial( this->giveLayerMaterial(layer) );
        static_cast< StructuralMaterial * >( layerMat )->giveRealStressVector_3d(answer, gp, strain, tStep);
    } else {
        OOFEM_ERROR("LayeredCrossSection :: giveRealStress_3d - Only cubes and wedges are meaningful for layered cross-sections");
    }
}


void
LayeredCrossSection :: giveRealStress_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("LayeredCrossSection :: giveRealStress_PlanteStrain - Not supported");
}


void
LayeredCrossSection :: giveRealStress_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("LayeredCrossSection :: giveRealStress_PlaneStress - Not supported");
}


void
LayeredCrossSection :: giveRealStress_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("LayeredCrossSection :: giveRealStress_1d - Not supported");
}


void
LayeredCrossSection :: giveRealStress_Beam2d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    FloatArray layerStrain, reducedLayerStress;
    StructuralElement *element = static_cast< StructuralElement * >( gp->giveElement() );
    LayeredCrossSectionInterface *interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );

    answer.resize(3);
    answer.zero();

    // perform integration over layers
    bottom = this->give(CS_BottomZCoord);
    top = this->give(CS_TopZCoord);

    if ( interface == NULL ) {
        _error("giveRealStresses - element with no layer support encountered");
    }

    for ( int i = 1; i <= numberOfLayers; i++ ) {
        GaussPoint *layerGp = this->giveSlaveGaussPoint(gp, i - 1);
        StructuralMaterial *layerMat = static_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(i) ) );

        // resolve current layer z-coordinate
        layerThick = this->layerThicks.at(i);
        layerWidth = this->layerWidths.at(i);
        layerZeta = layerGp->giveCoordinate(3);
        layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, layerGp, tStep);

        layerMat->giveRealStressVector_2dBeamLayer(reducedLayerStress, layerGp, layerStrain, tStep);

        answer.at(1) += reducedLayerStress.at(1) * layerWidth * layerThick;
        answer.at(2) += reducedLayerStress.at(1) * layerWidth * layerThick * layerZCoord;
        answer.at(3) += reducedLayerStress.at(2) * layerWidth * layerThick;
    }

    // now we must update master gp
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( gp->giveMaterial()->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
LayeredCrossSection :: giveRealStress_Beam3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("LayeredCrossSection :: giveRealStress_Beam3d - Not supported");
}


void
LayeredCrossSection :: giveRealStress_Plate(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    FloatArray layerStrain, reducedLayerStress;
    StructuralElement *element = static_cast< StructuralElement * >( gp->giveElement() );
    LayeredCrossSectionInterface *interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );

    answer.resize(5);
    answer.zero();

    // perform integration over layers
    bottom = this->give(CS_BottomZCoord);
    top = this->give(CS_TopZCoord);

    if ( interface == NULL ) {
        _error("giveRealStresses - element with no layer support encountered");
    }

    for ( int i = 1; i <= numberOfLayers; i++ ) {
        GaussPoint *layerGp = this->giveSlaveGaussPoint(gp, i - 1);
        StructuralMaterial *layerMat = static_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(i) ) );

        // resolve current layer z-coordinate
        layerThick = this->layerThicks.at(i);
        layerWidth = this->layerWidths.at(i);
        layerZeta = layerGp->giveCoordinate(3);
        layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, layerGp, tStep);

        layerMat->giveRealStressVector_PlateLayer(reducedLayerStress, layerGp, layerStrain, tStep);

        answer.at(1) += reducedLayerStress.at(1) * layerWidth * layerThick * layerZCoord;
        answer.at(2) += reducedLayerStress.at(2) * layerWidth * layerThick * layerZCoord;
        answer.at(3) += reducedLayerStress.at(5) * layerWidth * layerThick * layerZCoord;
        answer.at(4) += reducedLayerStress.at(4) * layerWidth * layerThick;
        answer.at(5) += reducedLayerStress.at(3) * layerWidth * layerThick;
    }

    // now we must update master gp
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( gp->giveMaterial()->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
LayeredCrossSection :: giveRealStress_Shell(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    FloatArray layerStrain, reducedLayerStress;
    StructuralElement *element = static_cast< StructuralElement * >( gp->giveElement() );
    LayeredCrossSectionInterface *interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );

    answer.resize(8);
    answer.zero();

    // perform integration over layers
    bottom = this->give(CS_BottomZCoord);
    top = this->give(CS_TopZCoord);

    if ( interface == NULL ) {
        _error("giveRealStresses - element with no layer support encountered");
    }

    for ( int i = 1; i <= numberOfLayers; i++ ) {
        GaussPoint *layerGp = this->giveSlaveGaussPoint(gp, i - 1);
        StructuralMaterial *layerMat = static_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(i) ) );

        // resolve current layer z-coordinate
        layerThick = this->layerThicks.at(i);
        layerWidth = this->layerWidths.at(i);
        layerZeta = layerGp->giveCoordinate(3);
        layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, layerGp, tStep);

        layerMat->giveRealStressVector_PlateLayer(reducedLayerStress, layerGp, layerStrain, tStep);

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
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( gp->giveMaterial()->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
LayeredCrossSection :: giveRealStress_MembraneRot(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("LayeredCrossSection :: giveRealStress_MembraneRot - Not supported in given cross-section (yet).");
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
        int layer = (gpnum - 1) / gpsperlayer + 1;
        StructuralMaterial *mat = static_cast< StructuralMaterial* >( domain->giveMaterial( this->giveLayerMaterial(layer) ) );
        if ( mat->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
            mat->giveStiffnessMatrix(answer, rMode, gp, tStep);
        } else {
            _error("giveMaterialStiffnessMatrixOf: unsupported StressStrainMode");
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
    bottom = this->give(CS_BottomZCoord);
    top = this->give(CS_TopZCoord);

    for ( int i = 1; i <= numberOfLayers; i++ ) {
        GaussPoint *layerGp = giveSlaveGaussPoint(gp, i - 1);

        ///@todo Just using the gp number doesn't nicely support more than 1 gp per layer. Must rethink.
        StructuralMaterial *mat = static_cast< StructuralMaterial* >( domain->giveMaterial( this->giveLayerMaterial(layerGp->giveNumber()) ) );
        mat->givePlateLayerStiffMtrx(layerMatrix, rMode, layerGp, tStep);

        //
        // resolve current layer z-coordinate
        //
        layerThick = this->layerThicks.at(i);
        layerWidth  = this->layerWidths.at(i);
        layerZeta   = layerGp->giveCoordinate(3);
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
    bottom = this->give(CS_BottomZCoord);
    top = this->give(CS_TopZCoord);

    for ( int i = 1; i <= numberOfLayers; i++ ) {
        GaussPoint *layerGp = giveSlaveGaussPoint(gp, i - 1);

        ///@todo The logic in this whole class is pretty messy to support both slave-gp's and normal gps. Rethinking the approach is necessary.
        /// Just using the gp number doesn't nicely support more than 1 gp per layer. Must rethink.
        StructuralMaterial *mat = static_cast< StructuralMaterial* >( domain->giveMaterial( this->giveLayerMaterial(layerGp->giveNumber()) ) );
        mat->givePlateLayerStiffMtrx(layerMatrix, rMode, layerGp, tStep);

        //
        // resolve current layer z-coordinate
        //
        layerThick = this->layerThicks.at(i);
        layerWidth  = this->layerWidths.at(i);
        layerZeta   = layerGp->giveCoordinate(3);
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
    bottom = this->give(CS_BottomZCoord);
    top = this->give(CS_TopZCoord);

    answer.resize(3, 3);
    answer.zero();

    for ( int i = 1; i <= numberOfLayers; i++ ) {
        GaussPoint *layerGp = giveSlaveGaussPoint(gp, i - 1);

        ///@todo The logic in this whole class is pretty messy to support both slave-gp's and normal gps. Rethinking the approach is necessary.
        /// Just using the gp number doesn't nicely support more than 1 gp per layer. Must rethink.
        StructuralMaterial *mat = static_cast< StructuralMaterial* >( domain->giveMaterial( this->giveLayerMaterial(layerGp->giveNumber()) ) );
        mat->give2dBeamLayerStiffMtrx(layerMatrix, rMode, layerGp, tStep);

        //
        // resolve current layer z-coordinate
        //
        layerThick = this->layerThicks.at(i);
        layerWidth  = this->layerWidths.at(i);
        layerZeta   = layerGp->giveCoordinate(3);
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
    OOFEM_ERROR("LayeredCrossSection :: give3dBeamStiffMtrx - Not implemented\n");
}


void
LayeredCrossSection :: giveMembraneRotStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("LayeredCrossSection :: giveMembraneRotStiffMtrx - Not implemented\n");
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
        _error("ImposeStressConstrainsOnGradient: gradientStressVector3d size mismatch");
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
        _error("ImposeStrainConstrainsOnGradient: gradientStrainVector3d size mismatch");
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
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, numberOfLayers, _IFT_LayeredCrossSection_nlayers);
    IR_GIVE_FIELD(ir, layerMaterials, _IFT_LayeredCrossSection_layermaterials);
    IR_GIVE_FIELD(ir, layerThicks, _IFT_LayeredCrossSection_thicks);
    IR_GIVE_OPTIONAL_FIELD(ir, layerWidths, _IFT_LayeredCrossSection_widths);

    if ( ( numberOfLayers != layerMaterials.giveSize() ) ||
        ( numberOfLayers != layerThicks.giveSize() ) )   //|| ( numberOfLayers != layerWidths.giveSize() ) ) 
    {
        _error("initializeFrom : numberOfLayers does not equal given number of thicknesses. ");
    }

    if ( numberOfLayers <= 0 ) {
        _error("instanciateFrom : numberOfLayers<= 0 is not allowed");
    }


    numberOfIntegrationPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfIntegrationPoints, _IFT_LayeredCrossSection_nintegrationpoints);

    // read z-coordinate of mid-surface measured from bottom layer
    midSurfaceZcoordFromBottom = 0.5*this->computeIntegralThick();  // Default: geometric midplane
    midSurfaceXiCoordFromBottom = 1.0; // add to IR
    IR_GIVE_OPTIONAL_FIELD(ir, midSurfaceZcoordFromBottom, _IFT_LayeredCrossSection_midsurf);

    this->setupLayerMidPlanes();

    return IRRT_OK;
}


void
LayeredCrossSection :: setupLayerMidPlanes() 
{
    // z-coord of each layer midplane measured from the global cross section z-coord
    this->layerMidZ.resize(this->numberOfLayers);
    double layerTopZ = -midSurfaceZcoordFromBottom;
    for ( int j = 1; j <= numberOfLayers; j++ ) {
        double thickness = this->layerThicks.at(j);
        layerTopZ += thickness;
        this->layerMidZ.at(j) = layerTopZ - thickness * 0.5; 
    }
}


int
LayeredCrossSection :: setupIntegrationPoints(IntegrationRule &irule, int npoints, Element *element)
{
    ///@todo We must send arrays for integration points instead of just a single scalar.
    if ( element->giveIntegrationDomain() == _Cube ) {
        int points1 = floor(cbrt( double( npoints ) ) + 0.5);
        // If numberOfIntegrationPoints > 0 then use that, otherwise use the element's default.
        return irule.SetUpPointsOnCubeLayers(points1, points1, this->numberOfIntegrationPoints ? numberOfIntegrationPoints : points1,
                                             element->giveMaterialMode(), this->layerThicks);
    } else if ( element->giveIntegrationDomain() == _Wedge ) {
        if ( npoints == 2 ) {
            return irule.SetUpPointsOnWedgeLayers(1, this->numberOfIntegrationPoints ? numberOfIntegrationPoints : 2,
                                                  element->giveMaterialMode(), this->layerThicks);
        } else {
            return irule.SetUpPointsOnWedgeLayers(3, this->numberOfIntegrationPoints ? numberOfIntegrationPoints : 3,
                                                  element->giveMaterialMode(), this->layerThicks);
        }
    } else {
        return irule.setUpIntegrationPoints(element->giveIntegrationDomain(), npoints, element->giveMaterialMode());
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
            _error("giveSlaveGaussPoint: no such layer defined");
        }

        // create new slave record in masterGp
        // (requires that this is friend of gp)
        double currentZTopCoord, currentZCoord,  bottom, top;
        FloatArray *zCoord, *masterCoords = masterGp->giveCoordinates();
        // resolve slave material mode
        MaterialMode slaveMode, masterMode = masterGp->giveMaterialMode();
        slaveMode = this->giveCorrespondingSlaveMaterialMode(masterMode);

        bottom = this->give(CS_BottomZCoord);
        top = this->give(CS_TopZCoord);

        masterGp->numberOfGp = this->numberOfLayers;                        // Generalize to multiple integration points per layer
        masterGp->gaussPointArray = new GaussPoint * [ numberOfLayers ];
        currentZTopCoord = -midSurfaceZcoordFromBottom;
        for ( int j = 0; j < numberOfLayers; j++ ) {
            currentZTopCoord += this->layerThicks.at(j + 1);
            currentZCoord = currentZTopCoord - this->layerThicks.at(j + 1) / 2.0; // z-coord of layer mid surface
            zCoord = new FloatArray(3);
            zCoord->zero();
            if ( masterCoords->giveSize() > 0 ) {
                zCoord->at(1) = masterCoords->at(1); // gp x-coord of mid surface
            }

            if ( masterCoords->giveSize() > 1 ) {
                zCoord->at(2) = masterCoords->at(2); // gp y-coord of mid surface
            }

            zCoord->at(3) = ( 2.0 * currentZCoord - top - bottom ) / ( top - bottom );
            // in gp - is stored isoparametric coordinate (-1,1) of z-coordinate
            //masterGp->gaussPointArray [ j ] = new GaussPoint(masterGp->giveIntegrationRule(), j + 1, zCoord, 0., slaveMode);
            
            // test - remove!
            masterGp->gaussPointArray [ j ] = new GaussPoint(masterGp->giveIntegrationRule(), j + 1, zCoord, 1.0, slaveMode);
        }

        slave = masterGp->gaussPointArray [ i ];
    }

    return slave;
}




void
LayeredCrossSection :: mapLayerGpCoordsToShellCoords(LayeredCrossSection *layeredCS, IntegrationRule **layerIntegrationRulesArray)
/*
  Maps the local xi-coord (z-coord) in each layer [-1,1] to the corresponding 
  xi-coord in the cross section coordinate system.
  Also renames the gp numbering from layerwise to glabal

        --------  1               --------  1
               |                         |  
               |                         |
        -------- -1       =>      --------  x
        --------  1               --------  x
               |                         |  
        -------- -1               -------- -1
*/
{
    double totalThickness = this->computeIntegralThick();
    int number = 1;
    for( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = layerIntegrationRulesArray [layer-1]; 

        for( int j = 1; j <= iRule->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(j-1);

            // Map local layer cs to local shell cs
            double zMid_i = layeredCS->giveLayerMidZ(layer);
            double xiMid_i = 1.0 - 2.0*(totalThickness - this->midSurfaceZcoordFromBottom - zMid_i)/totalThickness; 
            double xi = xiMid_i; 
            double xi2 = gp->coordinates->at(3)*layeredCS->giveLayerThickness(layer)/totalThickness;
            double xinew = xi+xi2;
            iRule->getIntegrationPoint(j-1)->coordinates->at(3) = xinew;
            iRule->getIntegrationPoint(j-1)->number = number;
            number++;
        }
    }

}



double
LayeredCrossSection :: computeIntegralThick()
//
// computes total thickness of receiver
//
{
    if ( totalThick == 0 ) {
        for ( int i = 1; i <= numberOfLayers; i++ ) {
            totalThick += layerThicks.at(i);
        }
    }

    return totalThick;
}


void
LayeredCrossSection :: printYourself()
// Prints the receiver on screen.
{
    printf("Cross Section with properties: \n");
    propertyDictionary->printYourself();
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
LayeredCrossSection :: saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *masterGp)
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
    StructuralMaterial *mat;
    GaussPoint *slaveGP;
    for ( int i = 1; i <= numberOfLayers; i++ ) {
        slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1);
        mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(i) ) );
        if ( ( iores = mat->saveIPContext(stream, mode, slaveGP) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType
LayeredCrossSection :: restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *masterGp)
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
    StructuralMaterial *mat;
    GaussPoint *slaveGP = NULL;
    for ( int i = 1; i <= numberOfLayers; i++ ) {
        // creates also slaves if they don't exists
        slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1);
        mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(i) ) );
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
        _error("giveCorrespondingSlaveMaterialMode : unsupported mode");
    }

    return _Unknown;
}


double
LayeredCrossSection :: give(CrossSectionProperty aProperty)
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

    return CrossSection :: give(aProperty);
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
        this->area = 0.0;
        for ( int i = 1; i <= numberOfLayers; i++ ) {
            this->area += this->layerThicks.at(i) * this->layerWidths.at(i);
        }
    }

    return area;
}

void
LayeredCrossSection :: computeStressIndependentStrainVector(FloatArray &answer,
                                                            GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
{
    StructuralElement *elem = static_cast< StructuralElement * >( gp->giveElement() );
    FloatArray et;

    elem->computeResultingIPTemperatureAt(et, stepN, gp, mode);

    if ( et.isNotEmpty() ) {
        _error("computeStressIndependentStrainVector: temperature loading not supported");
    }
}

bool LayeredCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode, int mat)
{
    for ( int i = 1; i <= this->numberOfLayers; i++ ) {
        if ( !this->domain->giveMaterial(this->giveLayerMaterial(i))->isCharacteristicMtrxSymmetric(rMode) ) {
            return false;
        }
    }
    return true;
}

} // end namespace oofem
