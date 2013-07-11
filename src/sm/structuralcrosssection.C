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

#include "structuralcrosssection.h"
#include "gausspoint.h"
#include "element.h"
#include "structuralmaterial.h"
#include "floatarray.h"

namespace oofem {

void
StructuralCrossSection ::  giveFirstPKStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &F, TimeStep *tStep)
{
    // This function returns the first Piola-Kirchoff stress in vector format
    // corresponding to a given deformation gradient according to the stress-deformation
    // mode stored in the each gp.

    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );
    if ( mat->hasMaterialModeCapability(mode) ) {
        mat->giveFirstPKStressVector(answer, gp, F, tStep);
        return;
    } else {
        _error("giveFirstPKStresses : unsupported MaterialMode");
    }
}


void
StructuralCrossSection ::  giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &F, TimeStep *tStep)
{
    // This function returns the Cauchy stress in vector format
    // corresponding to a given deformation gradient according to the stress-deformation
    // mode stored in the each gp.

    MaterialMode mode = gp->giveMaterialMode();
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );
    if ( mat->hasMaterialModeCapability(mode) ) {
        mat->giveCauchyStressVector(answer, gp, F, tStep);
        return;
    } else {
        _error("giveCauchyStresses : unsupported MaterialMode");
    }
}


void
StructuralCrossSection :: giveStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                   MatResponseMode rMode, GaussPoint *gp,
                                                   TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );
    mat->giveStiffnessMatrix_dPdF(answer, rMode, gp, tStep);
}


void
StructuralCrossSection :: giveStiffnessMatrix_dCde(FloatMatrix &answer,
                                                   MatResponseMode rMode, GaussPoint *gp,
                                                   TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );
    mat->giveStiffnessMatrix_dCde(answer, rMode, gp, tStep);
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
    if ( gradientStressVector3d->giveSize() != 6 ) {
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
        for ( int i = 2; i <= 6; i++ ) {
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
    if ( gradientStrainVector3d->giveSize() != 6 ) {
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
        for ( int i = 2; i <= 6; i++ ) {
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
StructuralCrossSection :: computeStressIndependentStrainVector(FloatArray &answer,
                                                               GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
//
// returns initial strain vector induced by stress independent effects
// like temperatue or shrinkage.
// takes into account form of load vector assumed by engngModel (Incremental or Total Load form).
//
{
    StructuralMaterial *mat = static_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );
    // add parts caused by  material
    mat->computeStressIndependentStrainVector(answer, gp, stepN, mode);
}
} // end namespace oofem
