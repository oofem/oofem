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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "layeredcrosssection.h"
#include "structuralelement.h"
#include "gausspnt.h"
#include "material.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "flotarry.h"
#include "contextioerr.h"

namespace oofem {
void
LayeredCrossSection ::  giveRealStresses(FloatArray &answer, MatResponseForm form,
                                         GaussPoint *gp,
                                         const FloatArray &totalStrain, TimeStep *tStep)
//
// this function returns a real stresses corresponding to
// given strainIncrement according to stressStrain mode stored
// in each gp.
// IMPORTANT:
//
{
    FloatArray stressVector3d;
    FloatArray layerStrain, fullLayerStrain, fullStressVect, stressVect;
    StructuralElement *element = ( StructuralElement * ) gp->giveElement();
    Material *layerMat;
    StructuralMaterialStatus *status;
    LayeredCrossSectionInterface *interface;

    if ( ( interface = ( LayeredCrossSectionInterface * ) element->giveInterface(LayeredCrossSectionInterfaceType) ) == NULL ) {
        _error("giveRealStresses - element with no layer support encountered");
    }

    GaussPoint *layerGp;

    for ( int i = 1; i <= numberOfLayers; i++ ) {
        // the question is whether this function should exist ?
        // if yes the element details will be hidden.
        // good idea also should be existence of element::GiveBmatrixOfLayer
        // and computing strains here - but first idea looks better
        layerGp = this->giveSlaveGaussPoint(gp, i - 1);
        layerMat = domain->giveMaterial( layerMaterials.at(i) );
        // but treating of geometric non-linearities may become more complicated
        // another aproach - use several functions with assumed
        // kinematic constraints
        //prevLayerStrain = (((StructuralMaterialStatus*) layerMat->giveStatus(layerGp))
        //  ->giveStrainVector());
        interface->computeStrainVectorInLayer(fullLayerStrain, gp, layerGp, tStep);
        this->giveReducedCharacteristicVector(layerStrain, layerGp, fullLayerStrain);

        /*
         * if (prevLayerStrain.giveSize()) {
         * layerStrainIncrement = layerStrain;
         * layerStrainIncrement.subtract(prevLayerStrain);
         * //layerStrainIncrement = prevLayerStrain->GiveCopy()->negated();
         * //layerStrainIncrement -> add (layerStrain);
         * } else {
         * layerStrainIncrement = layerStrain;
         * }
         */

        ( ( StructuralMaterial * ) layerMat )
        ->giveRealStressVector(stressVector3d, FullForm, layerGp, layerStrain, tStep);
        // reducedStressIncrement = this -> GiveReducedStressVector (gp, stressIncrement3d);
    }

    this->giveIntegrated3dShellStress(fullStressVect, gp);
    this->giveReducedCharacteristicVector(stressVect, gp, fullStressVect);
    answer = stressVect;
    status = ( StructuralMaterialStatus * ) ( gp->giveMaterial()->giveStatus(gp) );

    // now we must update master gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stressVect);
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
    this->giveMaterialStiffnessMatrixOf(answer, ReducedForm, rMode, gp,
                                        dynamic_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() ),
                                        tStep);
}


void
LayeredCrossSection :: giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
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
LayeredCrossSection :: giveMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                     MatResponseForm form,
                                                     MatResponseMode rMode,
                                                     GaussPoint *gp,
                                                     StructuralMaterial *mat,
                                                     TimeStep *tStep)

//
// only interface to material class, forcing returned matrix to be in form form.
//
{
    //Material *mat = gp->giveElement()->giveMaterial();
    MaterialMode mode = gp->giveMaterialMode();
    if ( ( mode == _2dPlate ) || ( mode == _3dShell ) || ( mode == _2dBeam ) ) {
        this->giveDerivedMaterialStiffnessMatrix(answer, form, rMode, gp, mat, tStep);
    } else if ( mat->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
        mat->giveCharacteristicMatrix(answer, form, rMode, gp, tStep);
    } else {
        _error("giveMaterialStiffnessMatrixOf: unsupported StressStrainMode");
    }
}


void
LayeredCrossSection :: giveDerivedMaterialStiffnessMatrix(FloatMatrix &answer,
                                                          MatResponseForm form,
                                                          MatResponseMode rMode,
                                                          GaussPoint *gp,
                                                          StructuralMaterial *mat,
                                                          TimeStep *tStep)
//
// return material stiffness matrix for derived types of stressStreinState
//
{
    MaterialMode mode = gp->giveMaterialMode();
    //
    // integral layer modes
    //
    if ( mode == _2dPlate ) {
        this->give2dPlateMaterialStiffnessMatrix(answer, form, rMode, gp, mat, tStep);
    } else if ( mode == _2dBeam ) {
        this->give2dBeamMaterialStiffnessMatrix(answer, form, rMode, gp, mat, tStep);
    } else if ( mode == _3dShell ) {
        this->give3dShellMaterialStiffness(answer, form, rMode, gp, mat, tStep);
    } else {
        _error("giveDerivedMaterialStiffnessMatrix: unsupported StressStrainMode");
    }
}


void
LayeredCrossSection :: give2dPlateMaterialStiffnessMatrix(FloatMatrix &answer,
                                                          MatResponseForm form,
                                                          MatResponseMode rMode,
                                                          GaussPoint *gp,
                                                          StructuralMaterial *mat,
                                                          TimeStep *tStep)

//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    FloatMatrix layerMatrix;
    GaussPoint *layerGp;
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    double layerZCoord2;
    int i, jj;

    if ( mode != _2dPlate ) {
        _error("give2dPlateMaterialStiffnessMatrix : unsupported mode");
    }

    // if (form != ReducedForm) error ("give2dPlateMaterialStiffnessMatrix : full form unsupported");

    if ( form == ReducedForm ) {
        answer.resize(5, 5);
        answer.zero();
    } else {
        answer.resize(8, 8);
        answer.zero();
    }

    // perform integration over layers
    this->computeIntegralThick(); // ensure that total thick has been computed
    bottom = -midSurfaceZcoordFromBottom;
    top    = totalThick - midSurfaceZcoordFromBottom;

    for ( i = 1; i <= numberOfLayers; i++ ) {
        layerGp = giveSlaveGaussPoint(gp, i - 1);
        this->giveLayerMaterialStiffnessMatrix(layerMatrix, FullForm, rMode, layerGp, tStep);

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
        jj = 3; // FullForm assumed
        if ( form == ReducedForm ) {
            jj = 0;
        }

        answer.at(1 + jj, 1 + jj) += layerMatrix.at(1, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(1 + jj, 2 + jj) += layerMatrix.at(1, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(1 + jj, 3 + jj) += layerMatrix.at(1, 6) * layerWidth * layerThick * layerZCoord2;

        answer.at(2 + jj, 1 + jj) += layerMatrix.at(2, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(2 + jj, 2 + jj) += layerMatrix.at(2, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(2 + jj, 3 + jj) += layerMatrix.at(2, 6) * layerWidth * layerThick * layerZCoord2;

        answer.at(3 + jj, 1 + jj) += layerMatrix.at(6, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(3 + jj, 2 + jj) += layerMatrix.at(6, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(3 + jj, 3 + jj) += layerMatrix.at(6, 6) * layerWidth * layerThick * layerZCoord2;

        // 2) shear terms qx = qxz, qy = qyz
        answer.at(4 + jj, 4 + jj) += layerMatrix.at(5, 5) * layerWidth * layerThick;
        answer.at(4 + jj, 5 + jj) += layerMatrix.at(5, 4) * layerWidth * layerThick;
        answer.at(5 + jj, 4 + jj) += layerMatrix.at(4, 5) * layerWidth * layerThick;
        answer.at(5 + jj, 5 + jj) += layerMatrix.at(4, 4) * layerWidth * layerThick;
    }
}


void
LayeredCrossSection :: give3dShellMaterialStiffness(FloatMatrix &answer, MatResponseForm form,
                                                    MatResponseMode rMode,
                                                    GaussPoint *gp,
                                                    StructuralMaterial *mat,
                                                    TimeStep *tStep)
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    //Material * mat = gp->giveElement()->giveMaterial();
    FloatMatrix layerMatrix;
    GaussPoint *layerGp;
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    // double zi, zi1;
    double layerZCoord2;
    int i;

    if ( mode != _3dShell ) {
        _error("give3dShellMaterialStiffness : unsupported mode");
    }

    // if (form != ReducedForm) error ("give3dShellMaterialStiffness : full form unsupported");

    answer.resize(8, 8);
    answer.zero();
    // perform integration over layers
    this->computeIntegralThick(); // ensure that total thick has been conputed
    bottom = -midSurfaceZcoordFromBottom;
    top    = totalThick - midSurfaceZcoordFromBottom;

    for ( i = 1; i <= numberOfLayers; i++ ) {
        layerGp = giveSlaveGaussPoint(gp, i - 1);
        this->giveLayerMaterialStiffnessMatrix(layerMatrix, FullForm, rMode, layerGp, tStep);
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
        answer.at(1, 3) += layerMatrix.at(1, 6) * layerWidth * layerThick;

        answer.at(2, 1) += layerMatrix.at(2, 1) * layerWidth * layerThick;
        answer.at(2, 2) += layerMatrix.at(2, 2) * layerWidth * layerThick;
        answer.at(2, 3) += layerMatrix.at(2, 6) * layerWidth * layerThick;

        answer.at(3, 1) += layerMatrix.at(6, 1) * layerWidth * layerThick;
        answer.at(3, 2) += layerMatrix.at(6, 2) * layerWidth * layerThick;
        answer.at(3, 3) += layerMatrix.at(6, 6) * layerWidth * layerThick;

        // 2) bending terms mx, my, mxy

        answer.at(4, 4) += layerMatrix.at(1, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(4, 5) += layerMatrix.at(1, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(4, 6) += layerMatrix.at(1, 6) * layerWidth * layerThick * layerZCoord2;

        answer.at(5, 4) += layerMatrix.at(2, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(5, 5) += layerMatrix.at(2, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(5, 6) += layerMatrix.at(2, 6) * layerWidth * layerThick * layerZCoord2;

        answer.at(6, 4) += layerMatrix.at(6, 1) * layerWidth * layerThick * layerZCoord2;
        answer.at(6, 5) += layerMatrix.at(6, 2) * layerWidth * layerThick * layerZCoord2;
        answer.at(6, 6) += layerMatrix.at(6, 6) * layerWidth * layerThick * layerZCoord2;

        // 3) shear terms qx, qy
        answer.at(7, 7) += layerMatrix.at(5, 5) * layerWidth * layerThick;
        answer.at(7, 8) += layerMatrix.at(5, 4) * layerWidth * layerThick;
        answer.at(8, 7) += layerMatrix.at(4, 5) * layerWidth * layerThick;
        answer.at(8, 8) += layerMatrix.at(4, 4) * layerWidth * layerThick;
    }
}


void
LayeredCrossSection :: give2dBeamMaterialStiffnessMatrix(FloatMatrix &answer,
                                                         MatResponseForm form,
                                                         MatResponseMode rMode,
                                                         GaussPoint *gp,
                                                         StructuralMaterial *mat,
                                                         TimeStep *tStep)
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    FloatMatrix layerMatrix;
    GaussPoint *layerGp;
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    double layerZCoord2;
    int i;

    if ( mode != _2dBeam ) {
        _error("give2dBeamMaterialStiffnessMatrix : unsupported mode");
    }

    // perform integration over layers
    this->computeIntegralThick(); // ensure that total thick has been conputed
    bottom = -midSurfaceZcoordFromBottom;
    top    = totalThick - midSurfaceZcoordFromBottom;

    if ( form == ReducedForm ) {
        answer.resize(3, 3);
    } else {
        answer.resize(8, 8);
    }

    answer.zero();

    for ( i = 1; i <= numberOfLayers; i++ ) {
        layerGp = giveSlaveGaussPoint(gp, i - 1);
        this->giveLayerMaterialStiffnessMatrix(layerMatrix, ReducedForm, rMode,
                                               layerGp, tStep);
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
        if ( form == ReducedForm ) {
            // 1) membrane terms sx
            answer.at(1, 1) += layerMatrix.at(1, 1) * layerWidth * layerThick;
            answer.at(1, 3) += layerMatrix.at(1, 2) * layerWidth * layerThick;
            // 2) bending terms my
            answer.at(2, 2) += layerMatrix.at(1, 1) * layerWidth * layerThick * layerZCoord2;
            answer.at(2, 3) += layerMatrix.at(1, 2) * layerWidth * layerThick * layerZCoord2;
            // 3) shear terms qx
            answer.at(3, 1) += layerMatrix.at(2, 1) * layerWidth * layerThick;
            answer.at(3, 3) += layerMatrix.at(2, 2) * layerWidth * layerThick;
        } else {
            // 1) membrane terms sx
            answer.at(1, 1) += layerMatrix.at(1, 1) * layerWidth * layerThick;
            answer.at(1, 7) += layerMatrix.at(1, 2) * layerWidth * layerThick;
            // 2) bending terms my
            answer.at(5, 5) += layerMatrix.at(1, 1) * layerWidth * layerThick * layerZCoord2;
            answer.at(5, 7) += layerMatrix.at(1, 2) * layerWidth * layerThick * layerZCoord2;
            // 3) shear terms qx
            answer.at(7, 1) += layerMatrix.at(2, 1) * layerWidth * layerThick;
            answer.at(7, 7) += layerMatrix.at(2, 2) * layerWidth * layerThick;
        }
    }
}


void
LayeredCrossSection :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                       const FloatArray &charVector3d)
//
// returns reduced charVector3d from full 3d  vector reduced
// to  vector required by gp->giveStressStrainMode()
//
// enhanced method in order to support cases with integral bending (2dplate, 3dshell..)
// in such cases full strain vector has the form:
// strainVector {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    int size = charVector3d.giveSize();

    if ( mode == _2dPlateLayer ) {
        if ( size != 6 ) {
            _error("giveReducedCharacteristicVector: stressVector3d size mismatch");
        }

        answer.resize(5);

        answer.at(1) = charVector3d.at(1);
        answer.at(2) = charVector3d.at(2);
        answer.at(3) = charVector3d.at(4);
        answer.at(4) = charVector3d.at(5);
        answer.at(5) = charVector3d.at(6);
    } else if ( mode == _2dBeamLayer ) {
        if ( size != 6 ) {
            _error("giveReducedCharacteristicVector: stressVector3d size mismatch");
        }

        answer.resize(2);

        answer.at(1) = charVector3d.at(1);
        answer.at(2) = charVector3d.at(5);
    } else if ( mode == _3dShellLayer ) {
        if ( size != 6 ) {
            _error("giveReducedCharacteristicVector: stressVector3d size mismatch");
        }

        answer.resize(5);

        answer.at(1) = charVector3d.at(1);
        answer.at(2) = charVector3d.at(2);
        answer.at(3) = charVector3d.at(4);
        answer.at(4) = charVector3d.at(5);
        answer.at(5) = charVector3d.at(6);
    } else if ( mode == _2dPlate ) {
        if ( size != 8 ) {
            _error("giveReducedCharacteristicVector: stressVector3d size mismatch");
        }

        answer.resize(5);

        answer.at(1) = charVector3d.at(4);
        answer.at(2) = charVector3d.at(5);
        answer.at(3) = charVector3d.at(6);
        answer.at(4) = charVector3d.at(7);
        answer.at(5) = charVector3d.at(8);
    } else if ( mode == _2dBeam ) {
        if ( size != 8 ) {
            _error("giveReducedCharacteristicVector: stressVector3d size mismatch");
        }

        answer.resize(3);

        answer.at(1) = charVector3d.at(1);
        answer.at(2) = charVector3d.at(4);
        answer.at(3) = charVector3d.at(7);
    } else if ( mode == _3dShell ) {
        if ( size != 8 ) {
            _error("giveReducedCharacteristicVector: stressVector3d size mismatch");
        }

        answer = charVector3d;
    } else {
        this->StructuralCrossSection :: giveReducedCharacteristicVector(answer, gp, charVector3d);
    }
}

void
LayeredCrossSection :: giveFullCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                    const FloatArray &charVector)
//
// returns reduced charVector3d from full 3d  vector reduced
// to  vector required by gp->giveStressStrainMode()
//
// enhaced method in order to support cases with integral bending (2dplate, 3dshell..)
// in such cases full strain vector has the form:
// strainVector {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    int size = charVector.giveSize();

    if ( mode == _2dPlateLayer ) {
        if ( size != 5 ) {
            _error("giveFullCharacteristicVector: stressVector size mismatch");
        }

        answer.resize(6);
        answer.zero();

        answer.at(1) = charVector.at(1);
        answer.at(2) = charVector.at(2);
        answer.at(4) = charVector.at(3);
        answer.at(5) = charVector.at(4);
        answer.at(6) = charVector.at(5);
    } else if ( mode == _2dBeamLayer ) {
        if ( size != 2 ) {
            _error("giveFullCharacteristicVector: stressVector size mismatch");
        }

        answer.resize(6);
        answer.zero();

        answer.at(1) = charVector.at(1);
        answer.at(5) = charVector.at(2);
    } else if ( mode == _3dShellLayer ) {
        if ( size != 5 ) {
            _error("giveFullCharacteristicVector: stressVector size mismatch");
        }

        answer.resize(6);
        answer.zero();

        answer.at(1) = charVector.at(1);
        answer.at(2) = charVector.at(2);
        answer.at(4) = charVector.at(3);
        answer.at(5) = charVector.at(4);
        answer.at(6) = charVector.at(5);
    } else if ( mode == _2dPlate ) {
        if ( size != 5 ) {
            _error("giveFullCharacteristicVector: stressVector size mismatch");
        }

        answer.resize(8);
        answer.zero();

        answer.at(4) = charVector.at(1);
        answer.at(5) = charVector.at(2);
        answer.at(6) = charVector.at(3);
        answer.at(7) = charVector.at(4);
        answer.at(8) = charVector.at(5);
    } else if ( mode == _2dBeam ) {
        if ( size != 3 ) {
            _error("giveFullCharacteristicVector: stressVector size mismatch");
        }

        answer.resize(8);
        answer.zero();

        answer.at(1) = charVector.at(1);
        answer.at(4) = charVector.at(2);
        answer.at(7) = charVector.at(3);
    } else if ( mode == _3dShell ) {
        if ( size != 8 ) {
            _error("giveFullCharacteristicVector: stressVector size mismatch");
        }

        answer = charVector;
    } else {
        this->StructuralCrossSection :: giveFullCharacteristicVector(answer, gp, charVector);
    }
}


FloatArray *
LayeredCrossSection :: imposeStressConstrainsOnGradient(GaussPoint *gp,
                                                        FloatArray *gradientStressVector3d)
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
    int i, size = gradientStressVector3d->giveSize();
    if ( size != 6 ) {
        _error("ImposeStressConstrainsOnGradient: gradientStressVector3d size mismatch");
    }

    switch ( mode ) {
    case _3dShellLayer:
        gradientStressVector3d->at(3) = 0.;
        break;
    case _2dPlateLayer:
        gradientStressVector3d->at(3) = 0.;
        break;
    case _2dBeamLayer:
        for ( i = 2; i <= 5; i++ ) {
            gradientStressVector3d->at(i) = 0.;
        }

        break;
    default:
        StructuralCrossSection :: imposeStressConstrainsOnGradient(gp, gradientStressVector3d);
    }

    return gradientStressVector3d;
}



FloatArray *
LayeredCrossSection :: imposeStrainConstrainsOnGradient(GaussPoint *gp,
                                                        FloatArray *gradientStrainVector3d)
//
// returns modified gradient of strain vector, which is used to
// compute plastic strain increment.
//
// imposes zeros on places, where zero strain occurs or energetically connected stress
// is prescribed to be zero.
//
{
    MaterialMode mode = gp->giveMaterialMode();
    int i, size = gradientStrainVector3d->giveSize();
    if ( size != 6 ) {
        _error("ImposeStrainConstrainsOnGradient: gradientStrainVector3d size mismatch");
    }

    switch ( mode ) {
    case _3dShellLayer:
        gradientStrainVector3d->at(3) = 0.;
        break;
    case _2dPlateLayer:
        gradientStrainVector3d->at(3) = 0.;
        break;
    case _2dBeamLayer:
        for ( i = 2; i <= 5; i++ ) {
            gradientStrainVector3d->at(i) = 0.;
        }

        break;
    default:
        StructuralCrossSection :: imposeStrainConstrainsOnGradient(gp, gradientStrainVector3d);
    }

    return gradientStrainVector3d;
}


void
LayeredCrossSection :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                            MaterialMode mmode, StructuralMaterial *mat) const
{
    //
    // this function returns mask of reduced(if form == ReducedForm)
    // or Full(if form==FullForm) stressStrain vector in full or
    // reduced StressStrainVector
    // acording to stressStrain mode of given gp.
    //
    // mask has size of reduced or full StressStrain Vector and  i-th component
    // is index to full or reduced StressStrainVector where corresponding
    // stressStrain resides.
    //
    int i;
    if ( mat->hasMaterialModeCapability(mmode) ) {
        mat->giveStressStrainMask(answer, form, mmode);
        return;
    } else {
        if ( form == ReducedForm ) {
            switch ( mmode ) {
            case _2dPlate:
                answer.resize(5);
                answer.at(1) = 4;
                answer.at(2) = 5;
                answer.at(3) = 6;
                answer.at(4) = 7;
                answer.at(5) = 8;
                break;
            case _2dBeam:
                answer.resize(3);
                answer.at(1) = 1;
                answer.at(2) = 5;
                answer.at(3) = 7;
                break;
            /*   case _3dRotContinuum:
             * indx->at(1) = 1;
             * indx->at(2) = 2;
             * indx->at(3) = 3;
             * indx->at(4) = 6;
             * break;*/
            case _3dShell:
                answer.resize(8);
                for ( i = 1; i <= 8; i++ ) {
                    answer.at(i) = i;
                }

                break;
            default:
                _error2( "giveStressStrainMask : unknown mode (%s)", __MaterialModeToString(mmode) );
            }
        } else if ( form == FullForm ) {
            switch ( mmode ) {
            case _2dPlate:
                answer.resize(8);
                answer.zero();
                answer.at(4) = 1;
                answer.at(5) = 2;
                answer.at(6) = 3;
                answer.at(7) = 4;
                answer.at(8) = 5;
                break;
            case _2dBeam:
                answer.resize(8);
                answer.zero();
                answer.at(1) = 1;
                answer.at(5) = 2;
                answer.at(7) = 3;
                break;
            /*   case _3dRotContinuum:
             * indx->at(1) = 1;
             * indx->at(2) = 2;
             * indx->at(3) = 3;
             * indx->at(6) = 4;
             * break;*/
            case _3dShell:
                answer.resize(8);
                answer.zero();
                for ( i = 1; i <= 8; i++ ) {
                    answer.at(i) = i;
                }

                break;
            default:
                _error2( "giveStressStrainMask : unknown mode (%s)", __MaterialModeToString(mmode) );
            }
        } else {
            _error("giveStressStrainMask : unknown form mode");
        }
    }
}



IRResultType
LayeredCrossSection :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, numberOfLayers, IFT_LayeredCrossSection_nlayers, "nlayers"); // Macro
    IR_GIVE_FIELD(ir, layerMaterials, IFT_LayeredCrossSection_layermaterials, "layermaterials"); // Macro
    IR_GIVE_FIELD(ir, layerThicks, IFT_LayeredCrossSection_thicks, "thicks"); // Macro
    IR_GIVE_FIELD(ir, layerWidths, IFT_LayeredCrossSection_widths, "widths"); // Macro

    // read z - coordinate of mid-surface measured from bottom layer
    IR_GIVE_FIELD(ir, midSurfaceZcoordFromBottom, IFT_LayeredCrossSection_midsurf, "midsurf"); // Macro

    if ( ( numberOfLayers != layerMaterials.giveSize() ) ||
        ( numberOfLayers != layerThicks.giveSize() )    ||
        ( numberOfLayers != layerWidths.giveSize() ) ) {
        _error("instanciateFrom : Array size mismatch ");
    }

    if ( numberOfLayers <= 0 ) {
        _error("instanciateFrom : numberOfLayers<= 0 is not alloved");
    }

    return IRRT_OK;
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
        double currentZTopCoord = 0., currentZCoord = 0.,  bottom, top;
        FloatArray *zCoord, *masterCoords = masterGp->giveCoordinates();
        // resolve slave material mode
        MaterialMode slaveMode, masterMode = masterGp->giveMaterialMode();
        slaveMode = this->giveCorrespondingSlaveMaterialMode(masterMode);

        this->computeIntegralThick(); // ensure that total thic has been conputed
        bottom = -midSurfaceZcoordFromBottom;
        top    = totalThick - midSurfaceZcoordFromBottom;

        masterGp->numberOfGp = this->numberOfLayers;
        masterGp->gaussPointArray = new GaussPoint * [ numberOfLayers ];
        currentZTopCoord = -midSurfaceZcoordFromBottom;
        for ( int j = 0; j < numberOfLayers; j++ ) {
            currentZTopCoord += this->layerThicks.at(j + 1);
            currentZCoord = currentZTopCoord - this->layerThicks.at(j + 1) / 2.0;
            zCoord = new FloatArray(3);
            zCoord->zero();
            if ( masterCoords->giveSize() > 0 ) {
                zCoord->at(1) = masterCoords->at(1);
            }

            if ( masterCoords->giveSize() > 1 ) {
                zCoord->at(2) = masterCoords->at(2);
            }

            zCoord->at(3) = ( 2.0 * ( currentZCoord ) - top - bottom ) / ( top - bottom );
            // in gp - is stored isoparametric coordinate (-1,1) of z-coordinate
            masterGp->gaussPointArray [ j ] = new GaussPoint(masterGp->giveIntegrationRule(), j + 1, zCoord, 0., slaveMode);
        }

        slave = masterGp->gaussPointArray [ i ];
    }

    return slave;
}


double
LayeredCrossSection :: computeIntegralThick()
//
// computes total thick of receiver
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
    printf("Cross Section with properties : \n");
    propertyDictionary->printYourself();
    printf("Layer Materials: \n");
    layerMaterials.printYourself();
    printf("Layer Thicks   : \n");
    layerThicks.printYourself();
    printf("Layer Widths   : \n");
    layerWidths.printYourself();
    printf("MidSurfaceZCoordinate from bottom %f\n", midSurfaceZcoordFromBottom);
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
// returns cooresponding slave material mode to master mode
//
{
    if ( masterMode == _2dPlate ) {
        return _2dPlateLayer;
    } else if ( masterMode == _2dBeam ) {
        return _2dBeamLayer;
    } else if ( masterMode == _3dShell ) {
        return _3dShellLayer;
    } else {
        _error("giveCorrespondingSlaveMaterialMode : unsupported mode");
    }

    return _Unknown;
}


void
LayeredCrossSection :: giveIntegrated3dShellStress(FloatArray &answer, GaussPoint *masterGp)
//
// computes integral internal forces for current mode
// this functions assumes that slave gp's have updatet stresses.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    Material *layerMat;
    StructuralMaterialStatus *layerStatus;
    FloatArray layerStress, reducedLayerStress;
    GaussPoint *layerGp;
    double layerThick, layerWidth, layerZCoord, top, bottom, layerZeta;
    int i;

    answer.resize(8);
    // perform integration over layers
    this->computeIntegralThick(); // ensure that total thick has been conputed
    bottom = -midSurfaceZcoordFromBottom;
    top    = totalThick - midSurfaceZcoordFromBottom;

    for ( i = 1; i <= numberOfLayers; i++ ) {
        layerGp = giveSlaveGaussPoint(masterGp, i - 1);
        layerMat = domain->giveMaterial( layerMaterials.at(i) );
        layerStatus = ( ( StructuralMaterialStatus * ) layerMat->giveStatus(layerGp) );

        if ( layerStatus->giveTempStressVector().giveSize() ) { // there exist total sress in gp
            reducedLayerStress = layerStatus->giveTempStressVector();
            giveFullCharacteristicVector(layerStress, layerGp, reducedLayerStress);
        } else { // no total stress
            continue; // skip gp without stress
        }

        //
        // resolve current layer z-coordinate
        //
        layerThick = this->layerThicks.at(i);
        layerWidth  = this->layerWidths.at(i);
        layerZeta   = layerGp->giveCoordinate(3);
        layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );
        //
        // perform integration
        //
        // 1) membrane terms sx, sy, sxy
        answer.at(1) += layerStress.at(1) * layerWidth * layerThick;
        answer.at(2) += layerStress.at(2) * layerWidth * layerThick;
        answer.at(3) += layerStress.at(6) * layerWidth * layerThick;
        // 2) bending terms mx, my, mxy
        answer.at(4) += layerStress.at(1) * layerWidth * layerThick * layerZCoord;
        answer.at(5) += layerStress.at(2) * layerWidth * layerThick * layerZCoord;
        answer.at(6) += layerStress.at(6) * layerWidth * layerThick * layerZCoord;
        // 3) shear terms qx, qy
        answer.at(7) += layerStress.at(5) * layerWidth * layerThick;
        answer.at(8) += layerStress.at(4) * layerWidth * layerThick;
    }
}


double
LayeredCrossSection :: give(CrossSectionProperty aProperty)
{
    if ( aProperty == CS_Thickness ) {
        return this->computeIntegralThick();
    } else if ( aProperty == CS_TopZCoord )   {
        this->computeIntegralThick();
        return totalThick - midSurfaceZcoordFromBottom;
    } else if ( aProperty == CS_BottomZCoord )   {
        return -midSurfaceZcoordFromBottom;
    } else if ( aProperty == CS_Area )   {
        return this->giveArea();
    }

    return CrossSection :: give(aProperty);
}


double
LayeredCrossSection :: giveArea()
{
    int i;
    if ( this->area <= 0.0 ) {
        this->area = 0.0;
        for ( i = 1; i <= numberOfLayers; i++ ) {
            this->area += this->layerThicks.at(i) * this->layerWidths.at(i);
        }
    }

    return area;
}


void
LayeredCrossSection :: giveLayerMaterialStiffnessMatrix(FloatMatrix &layerMatrix, MatResponseForm form,
                                                        MatResponseMode rMode, GaussPoint *layerGp,
                                                        TimeStep *tStep)
{
    this->giveMaterialStiffnessMatrixOf(layerMatrix, form, rMode, layerGp,
                                        dynamic_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at( layerGp->giveNumber() ) ) ),
                                        tStep);
}


void
LayeredCrossSection :: computeStressIndependentStrainVector(FloatArray &answer,
                                                            GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
{
    StructuralElement *elem = ( StructuralElement * ) gp->giveElement();
    FloatArray et;

    elem->computeResultingIPTemperatureAt(et, stepN, gp, mode);

    if ( et.isNotEmpty() ) {
        _error("computeStressIndependentStrainVector: temperature loading not supported");
    }
}
} // end namespace oofem
