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

#include "structuralcrosssection.h"
#include "gausspoint.h"
#include "element.h"
#include "structuralmaterial.h"
#include "floatarray.h"
#include "structuralms.h"

namespace oofem {
void
StructuralCrossSection :: giveRealStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _2dBeam ) {
        this->giveGeneralizedStress_Beam2d(answer, gp, strain, tStep);
    } else if ( mode == _3dBeam ) {
        this->giveGeneralizedStress_Beam3d(answer, gp, strain, tStep);
    } else if ( mode == _2dPlate ) {
        this->giveGeneralizedStress_Plate(answer, gp, strain, tStep);
    } else if ( mode == _3dShell ) {
        this->giveGeneralizedStress_Shell(answer, gp, strain, tStep);
    } else if ( mode == _3dMat ) {
        this->giveRealStress_3d(answer, gp, strain, tStep);
    } else if ( mode == _PlaneStrain ) {
        this->giveRealStress_PlaneStrain(answer, gp, strain, tStep);
    } else if ( mode == _PlaneStress ) {
        this->giveRealStress_PlaneStress(answer, gp, strain, tStep);
    } else if ( mode == _1dMat ) {
        this->giveRealStress_1d(answer, gp, strain, tStep);
    } else {
        // This should never happen ?
        ///@todo this part only works for simple cross section and will be removed soon when new interface elements are done /JB
        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( gp->giveMaterial() );
        if ( mat->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
            mat->giveRealStressVector(answer, gp, strain, tStep);
        } else {
            OOFEM_ERROR("giveRealStresses : unsupported mode");
        }
    }
}

//
//void
//StructuralCrossSection :: giveFirstPKStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep)
//{
//    // This function returns the first Piola-Kirchoff stress in vector format
//    // corresponding to a given deformation gradient according to the stress-deformation
//    // mode stored in the each gp.
//
//    MaterialMode mode = gp->giveMaterialMode();
//    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
//
//    if ( mode == _3dMat ) {
//        this->giveFirstPKStress_3d(answer, gp, reducedvF, tStep);
//    } else if ( mode == _PlaneStrain ) {
//        this->giveFirstPKStress_PlaneStrain(answer, gp, reducedvF, tStep);
//    } else if ( mode == _PlaneStress ) {
//        this->giveFirstPKStress_PlaneStress(answer, gp, reducedvF, tStep);
//    } else if ( mode == _1dMat ) {
//        this->giveFirstPKStress_1d(answer, gp, reducedvF, tStep);
//    } else {
//        OOFEM_ERROR("StructuralCrossSection :: giveStiffnessMatrix_dPdF : unknown mode (%s)", __MaterialModeToString(mode) );
//    }
//}
//
//
//void
//StructuralCrossSection :: giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep)
//{
//    // This function returns the Cauchy stress in vector format
//    // corresponding to a given deformation gradient according to the stress-deformation
//    // mode stored in the each gp.
//
//    MaterialMode mode = gp->giveMaterialMode();
//    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
//
//    if ( mode == _3dMat ) {
//        this->giveCauchyStress_3d(answer, gp, reducedvF, tStep);
//    } else if ( mode == _PlaneStrain ) {
//        this->giveCauchyStress_PlaneStrain(answer, gp, reducedvF, tStep);
//    } else if ( mode == _PlaneStress ) {
//        this->giveCauchyStress_PlaneStress(answer, gp, reducedvF, tStep);
//    } else if ( mode == _1dMat ) {
//        this->giveCauchyStress_1d(answer, gp, reducedvF, tStep);
//    }
//}
//
//
//void
//StructuralCrossSection :: giveStiffnessMatrix_dPdF(FloatMatrix &answer,
//                                                   MatResponseMode rMode, GaussPoint *gp,
//                                                   TimeStep *tStep)
//{
//    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
//
//    MaterialMode mode = gp->giveMaterialMode();
//    if ( mode == _3dMat ) {
//        this->give3dMaterialStiffnessMatrix_dPdF(answer, rMode, gp, tStep);
//    } else if ( mode == _PlaneStress ) {
//        this->givePlaneStressStiffMtrx_dPdF(answer, rMode, gp, tStep);
//    } else if ( mode == _PlaneStrain ) {
//        this->givePlaneStrainStiffMtrx_dPdF(answer, rMode, gp, tStep);
//    } else if ( mode == _1dMat ) {
//        this->give1dStressStiffMtrx_dPdF(answer, rMode, gp, tStep);
//    } else {
//        OOFEM_ERROR("StructuralCrossSection :: giveStiffnessMatrix_dPdF : unknown mode (%s)", __MaterialModeToString(mode) );
//    }
//}
//
//
//void
//StructuralCrossSection :: giveStiffnessMatrix_dCde(FloatMatrix &answer,
//                                                   MatResponseMode rMode, GaussPoint *gp,
//                                                   TimeStep *tStep)
//{
//    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
//
//    MaterialMode mode = gp->giveMaterialMode();
//    if ( mode == _3dMat ) {
//        this->give3dMaterialStiffnessMatrix_dCde(answer, rMode, gp, tStep);
//    } else if ( mode == _PlaneStress ) {
//        this->givePlaneStressStiffMtrx_dCde(answer, rMode, gp, tStep);
//    } else if ( mode == _PlaneStrain ) {
//        this->givePlaneStrainStiffMtrx_dCde(answer, rMode, gp, tStep);
//    } else if ( mode == _1dMat ) {
//        this->give1dStressStiffMtrx_dCde(answer, rMode, gp, tStep);
//    } else {
//        OOFEM_ERROR("StructuralCrossSection :: giveStiffnessMatrix_dCde : unknown mode (%s)", __MaterialModeToString(mode) );
//    }
//}


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
        OOFEM_ERROR("ImposeStressConstrainsOnGradient: gradientStressVector3d size mismatch");
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
        OOFEM_ERROR( "ImposeStressConstrainsOnGradient: unknown mode (%s)", __MaterialModeToString(mode) );
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
        OOFEM_ERROR("ImposeStrainConstrainsOnGradient: gradientStrainVector3d size mismatch");
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
        OOFEM_ERROR( "ImposeStrainConstrainsOnGradient: unknown mode (%s)", __MaterialModeToString(mode) );
        break;
    }

    return gradientStrainVector3d;
}


IRResultType
StructuralCrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    return IRRT_OK;
}

///@todo make static? /JB - No. Remove it. / Mikael
int
StructuralCrossSection :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain  || mode == _1dMat ||
           mode == _PlateLayer || mode == _2dBeamLayer || mode == _Fiber;
}
} // end namespace oofem
