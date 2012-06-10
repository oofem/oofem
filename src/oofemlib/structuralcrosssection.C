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

#include "structuralcrosssection.h"
#include "gausspnt.h"
#include "element.h"
#include "structuralmaterial.h"
#include "flotarry.h"

namespace oofem {
void
StructuralCrossSection ::  giveRealStresses(FloatArray &answer, MatResponseForm form,
                                            GaussPoint *gp,
                                            const FloatArray &strain,
                                            TimeStep *tStep)
//
// this function returns a real stresses corresponding to
// given strain according to stressStrain mode stored
// in each gp.
// IMPORTANT:
//
{
    MaterialMode mode = gp->giveMaterialMode();
    Material *mat = gp->giveElement()->giveMaterial();

    if ( mat->hasMaterialModeCapability(mode) ) {
        ( ( StructuralMaterial * ) gp->giveElement()->giveMaterial() )
        ->giveRealStressVector(answer, form, gp, strain, tStep);
        return;
    } else {
        _error("giveRealStresses : unsupported mode");
    }
}


void
StructuralCrossSection :: giveCharMaterialStiffnessMatrix(FloatMatrix &answer,
                                                          MatResponseMode rMode, GaussPoint *gp,
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
StructuralCrossSection :: giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                            MatResponseForm form,
                                                            MatResponseMode rMode,
                                                            GaussPoint *gp, StructuralMaterial *mat,
                                                            TimeStep *tStep)
//
// only interface to material class, forcing returned matrix to be in reduced form.
//
{
    this->giveMaterialStiffnessMatrixOf(answer, form, rMode, gp, mat, tStep);
}


void
StructuralCrossSection :: giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                                      MatResponseForm form,
                                                      MatResponseMode rMode,
                                                      GaussPoint *gp,
                                                      TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );
    this->giveMaterialStiffnessMatrixOf(answer, form, rMode, gp,
                                        mat, tStep);
}



void
StructuralCrossSection :: giveMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                        MatResponseForm form,
                                                        MatResponseMode rMode,
                                                        GaussPoint *gp,
                                                        StructuralMaterial *mat,
                                                        TimeStep *tStep)
//
// only interface to material class, forcing returned matrix to be in reduced form.
//
{
    // Material *mat = gp->giveElement()->giveMaterial();
    ( ( StructuralMaterial * ) mat )->giveCharacteristicMatrix(answer, form, rMode, gp, tStep);
}



void
StructuralCrossSection :: giveCharMaterialComplianceMatrix(FloatMatrix &answer,
                                                           MatResponseMode rMode, GaussPoint *gp,
                                                           TimeStep *tStep)
{
    /* returns compliance matrix according to stress strain mode in gp */
    FloatMatrix redInvAnswer;

    this->giveCharMaterialStiffnessMatrix(redInvAnswer, rMode, gp, tStep);
    answer.beInverseOf(redInvAnswer);
}



void
StructuralCrossSection :: giveCharMaterialComplianceMatrixOf(FloatMatrix &answer,
                                                             MatResponseForm form,
                                                             MatResponseMode rMode,
                                                             GaussPoint *gp, StructuralMaterial *mat,
                                                             TimeStep *tStep)
{
    /* returns compliance matrix according to stress strain mode in gp */
    FloatMatrix redInvAnswer, redAnswer;
    IntArray mask;

    this->giveCharMaterialStiffnessMatrixOf(redInvAnswer, ReducedForm, rMode, gp, mat,
                                            tStep);
    redAnswer.beInverseOf(redInvAnswer);

    if ( form == FullForm ) {
        this->giveStressStrainMask( mask, ReducedForm, gp->giveMaterialMode(),
                                   static_cast< StructuralMaterial * >( gp->giveMaterial() ) );
        answer.beSubMatrixOfSizeOf(redAnswer, mask, 6);
    } else if ( form == ReducedForm ) {
        answer = redAnswer;
    } else {
        _error("giveCharMaterialComplianceMatrix - unsupported form mode");
    }
}

void
StructuralCrossSection :: giveFullCharacteristicVector(FloatArray &answer,
                                                       GaussPoint *gp,
                                                       const FloatArray &strainVector)
//
// returns full 3d general strain vector from strainVector in reducedMode
// based on StressStrainMode in gp. Included are strains which
// perform nonzero work.
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
// you must assigng your stress strain mode to one of the folloving modes (or add new)
// FullForm of MaterialStiffnessMatrix must have the same form.
//
{
    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );

    if ( ( mode == _3dMat ) || ( mode == _3dMicroplane ) ) {
        answer = strainVector;
        return;
    } else {
        mat->giveFullCharacteristicVector(answer, gp, strainVector);
    }
}



void
StructuralCrossSection :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                          const FloatArray &charVector3d)
//
// returns reduced stressVector or strainVector from full 3d vector reduced
// to vector required by gp->giveStressStrainMode()
//
{
    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );

    if ( ( mode == _3dMat ) || ( mode == _3dMicroplane ) ) {
        answer = charVector3d;
        return;
    } else {
        mat->giveReducedCharacteristicVector(answer, gp, charVector3d);
    }
}




FloatArray *
StructuralCrossSection :: imposeStressConstrainsOnGradient(GaussPoint *gp,
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

    if ( mode == _3dMat ) {
        return gradientStressVector3d;
    }


    switch ( mode ) {
    case _PlaneStress:
        gradientStressVector3d->at(3) = 0.;
        gradientStressVector3d->at(4) = 0.;
        gradientStressVector3d->at(5) = 0.;
        break;
    case _PlaneStrain:
        // gradientStressVector3d ->at(3) = 0.;
        gradientStressVector3d->at(4) = 0.;
        gradientStressVector3d->at(5) = 0.;
        break;
    case _1dMat:
        for ( i = 2; i <= 6; i++ ) {
            gradientStressVector3d->at(i) = 0.;
        }

        break;
    default:
        _error2( "ImposeStressConstrainsOnGradient: unknown mode (%s)", __MaterialModeToString(mode) );
        break;
    }

    return gradientStressVector3d;
}



FloatArray *
StructuralCrossSection :: imposeStrainConstrainsOnGradient(GaussPoint *gp,
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

    if ( mode == _3dMat ) {
        return gradientStrainVector3d;
    }


    switch ( mode ) {
    case _PlaneStress:
        gradientStrainVector3d->at(3) = 0.;
        gradientStrainVector3d->at(4) = 0.;
        gradientStrainVector3d->at(5) = 0.;
        break;
    case _PlaneStrain:
        gradientStrainVector3d->at(3) = 0.;
        gradientStrainVector3d->at(4) = 0.;
        gradientStrainVector3d->at(5) = 0.;
        break;
    case _1dMat:
        for ( i = 2; i <= 6; i++ ) {
            gradientStrainVector3d->at(i) = 0.;
        }

        break;
    default:
        _error2( "ImposeStrainConstrainsOnGradient: unknown mode (%s)", __MaterialModeToString(mode) );
        break;
    }

    return gradientStrainVector3d;
}


void
StructuralCrossSection :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
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
    //MaterialMode mode = gp-> giveMaterialMode ();
    //StructuralMaterial * mat = (StructuralMaterial*) gp->giveElement()->giveMaterial();
    if ( mat->hasMaterialModeCapability(mmode) ) {
        mat->giveStressStrainMask(answer, form, mmode);
        return;
    } else {
        _error("giveStressStrainMask : unsupported mode");
    }
}


void
StructuralCrossSection :: computeStressIndependentStrainVector(FloatArray &answer,
                                                               GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
//
// returns initial strain vector induced by stress independent effects
// like temperatue or shrinkage.
// takes into account form of load vector assumed by engngModel (Incremental or Total Load form).
//
{
    StructuralMaterial *mat = ( StructuralMaterial * ) gp->giveElement()->giveMaterial();
    FloatArray e0, fullAnswer;

    //
    // add parts caused by  material
    //

    mat->computeStressIndependentStrainVector(answer, gp, stepN, mode);
}
} // end namespace oofem
