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

namespace oofem {





void
StructuralCrossSection :: giveRealStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _2dBeam ) {
        this->giveRealStress_Beam2d(answer, gp, strain, tStep);
    } else if ( mode == _3dBeam ) {
        this->giveRealStress_Beam3d(answer, gp, strain, tStep);
    } else if ( mode == _2dPlate ) {
        this->giveRealStress_Plate(answer, gp, strain, tStep);
    } else if ( mode == _3dShell ) {
        this->giveRealStress_Shell(answer, gp, strain, tStep);
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
        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );
        if ( mat->hasMaterialModeCapability(gp->giveMaterialMode()) ) {
            mat->giveRealStressVector(answer, gp, strain, tStep);
        } else {
            _error("giveRealStresses : unsupported mode");
        }
    }
}


void
StructuralCrossSection :: giveFirstPKStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep)
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
        OOFEM_ERROR2("StructuralCrossSection :: giveStiffnessMatrix_dPdF : unknown mode (%s)", __MaterialModeToString(mode) );
    }
}


void
StructuralCrossSection :: giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep)
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
StructuralCrossSection :: giveStiffnessMatrix_dPdF(FloatMatrix &answer,
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
        OOFEM_ERROR2("StructuralCrossSection :: giveStiffnessMatrix_dPdF : unknown mode (%s)", __MaterialModeToString(mode) );
    }
}


void
StructuralCrossSection :: giveStiffnessMatrix_dCde(FloatMatrix &answer,
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
        OOFEM_ERROR2("StructuralCrossSection :: giveStiffnessMatrix_dCde : unknown mode (%s)", __MaterialModeToString(mode) );
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
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
    // add parts caused by  material
    mat->computeStressIndependentStrainVector(answer, gp, stepN, mode);
}



IRResultType
StructuralCrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // Read a cohesive zone material
    IR_GIVE_OPTIONAL_FIELD(ir, this->materialNumber, _StructuralCrossSection_MaterialNumber);

    // Read a cohesive zone material
    IR_GIVE_OPTIONAL_FIELD(ir, this->czMaterialNumber, _StructuralCrossSection_czMaterialNumber);

    return IRRT_OK;
}


int
StructuralCrossSection :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain  || mode == _1dMat ||
           mode == _PlateLayer || mode == _2dBeamLayer || mode == _Fiber;
}


int
StructuralCrossSection :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    Material *mat = this->giveDomain()->giveMaterial(this->materialNumber);
    if ( !dynamic_cast< StructuralMaterial * >( mat ) ) {
        _warning2("checkConsistency : material %s without structural support", mat->giveClassName());
        result = 0;
    }

    return result;
}


Material 
*StructuralCrossSection :: giveMaterial(IntegrationPoint *ip) 
{ 
    if ( ip->giveElement()->MAT_GIVEN_BY_CS ) {
        return this->giveDomain()->giveMaterial( this->giveMaterialNumber() ); 
    } else {
        return ip->giveElement()->giveMaterial();
    }
}

} // end namespace oofem
