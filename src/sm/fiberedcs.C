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

#include "fiberedcs.h"
#include "structuralelement.h"
#include "gausspnt.h"
#include "material.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "flotarry.h"
#include "verbose.h"
#include "contextioerr.h"

namespace oofem {
void
FiberedCrossSection ::  giveRealStresses(FloatArray &answer, MatResponseForm form,
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
    FloatArray fullFiberStrain, fiberStrain, *fullStressVect, stressVect;
    StructuralElement *element = ( StructuralElement * ) gp->giveElement();
    Material *fiberMat;
    StructuralMaterialStatus *status;
    FiberedCrossSectionInterface *interface;

    if ( ( interface = ( FiberedCrossSectionInterface * ) element->giveInterface(FiberedCrossSectionInterfaceType) ) == NULL ) {
        _error("giveRealStresses - element with no fiber support encountered");
    }

    GaussPoint *fiberGp;

    for ( int i = 1; i <= numberOfFibers; i++ ) {
        // the question is whether this function should exist ?
        // if yes the element details will be hidden.
        // good idea also should be existence of element::GiveBmatrixOfLayer
        // and computing strains here - but first idea looks better
        fiberGp = this->giveSlaveGaussPoint(gp, i - 1);
        fiberMat = domain->giveMaterial( fiberMaterials.at(i) );
        // but treating of geometric non-linearities may become more complicated
        // another approach - use several functions with assumed
        // kinematic constraints
        interface->FiberedCrossSectionInterface_computeStrainVectorInFiber(fullFiberStrain, gp, fiberGp, tStep);
        this->giveReducedCharacteristicVector(fiberStrain, fiberGp, fullFiberStrain);


        ( ( StructuralMaterial * ) fiberMat )->giveRealStressVector(stressVector3d, FullForm, fiberGp, fiberStrain, tStep);
    }

    fullStressVect = this->GiveIntegrated3dBeamStress(gp);
    giveReducedCharacteristicVector(stressVect, gp, * fullStressVect);
    delete fullStressVect;
    answer = stressVect;
    status = ( StructuralMaterialStatus * ) ( gp->giveMaterial()->giveStatus(gp) );

    // now we must update master gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stressVect);
}


void
FiberedCrossSection :: giveCharMaterialStiffnessMatrix(FloatMatrix &answer,
                                                       MatResponseMode rMode,
                                                       GaussPoint *gp,
                                                       TimeStep *tStep)
//
// only interface to material class, forcing returned matrix to be in reduced form.
//
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( gp->giveElement()->giveMaterial() );
    this->giveMaterialStiffnessMatrixOf(answer, ReducedForm, rMode, gp,
                                        mat, tStep);
}


void
FiberedCrossSection :: giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
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
FiberedCrossSection :: giveMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                     MatResponseForm form,
                                                     MatResponseMode rMode,
                                                     GaussPoint *gp,
                                                     StructuralMaterial *mat,
                                                     TimeStep *tStep)

//
// only interface to material class, forcing returned matrix to be in form form.
//
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dBeam ) {
        this->give3dBeamMaterialStiffnessMatrix(answer, form, rMode, gp, mat, tStep);
    } else if ( mat->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
        mat->giveCharacteristicMatrix(answer, form, rMode, gp, tStep);
    } else {
        _error("giveMaterialStiffnessMatrixOf: unsupported StressStrainMode");
    }
}


void
FiberedCrossSection :: give3dBeamMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                                         MatResponseMode rMode,
                                                         GaussPoint *gp,
                                                         StructuralMaterial *mat,
                                                         TimeStep *tStep)
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain fiber vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x, gamma_xz, gamma_xy, \der{phi_x}{x}, kappa_y, kappa_z}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    //Material * mat = gp->giveElement()->giveMaterial();
    FloatMatrix fiberMatrix;
    GaussPoint *fiberGp;
    double fiberThick, fiberWidth, fiberZCoord, fiberYCoord;
    double fiberZCoord2, fiberYCoord2, Ip = 0.0, A = 0.0, Ik, G = 0.0;
    int i;

    if ( mode != _3dBeam ) {
        _error("give3dBeamMaterialStiffness : unsupported mode");
    }

    // if (form != ReducedForm) error ("give3dShellMaterialStiffness : full form unsupported");

    answer.resize(6, 6);
    answer.zero();
    // perform integration over layers

    for ( i = 1; i <= numberOfFibers; i++ ) {
        fiberGp = giveSlaveGaussPoint(gp, i - 1);
        this->giveFiberMaterialStiffnessMatrix(fiberMatrix, FullForm, rMode, fiberGp, tStep);
        //
        // resolve current layer z-coordinate
        //
        fiberThick  = this->fiberThicks.at(i);
        fiberWidth  = this->fiberWidths.at(i);
        fiberZCoord = fiberZcoords.at(i);
        fiberYCoord = fiberYcoords.at(i);
        fiberYCoord2 = fiberYCoord * fiberYCoord;
        fiberZCoord2 = fiberZCoord * fiberZCoord;
        //
        // perform integration
        //
        // 1) membrane terms N, Qz, Qy
        answer.at(1, 1) += fiberMatrix.at(1, 1) * fiberWidth * fiberThick;

        answer.at(2, 2) += fiberMatrix.at(5, 5) * fiberWidth * fiberThick;

        answer.at(3, 3) += fiberMatrix.at(6, 6) * fiberWidth * fiberThick;

        // 2) bending terms mx, my, mz

        Ip             += fiberWidth * fiberThick * fiberZCoord2 + fiberWidth * fiberThick * fiberYCoord2;
        A              += fiberWidth * fiberThick;
        G               = fiberMatrix.at(5, 5) * fiberWidth * fiberThick;

        answer.at(5, 5) += fiberMatrix.at(1, 1) * fiberWidth * fiberThick * fiberZCoord2;
        answer.at(6, 6) += fiberMatrix.at(1, 1) * fiberWidth * fiberThick * fiberYCoord2;
    }

    G /= A;
    Ik = A * A * A * A / ( 40.0 * Ip );
    answer.at(4, 4) = G * Ik;
}


void
FiberedCrossSection :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                       const FloatArray &charVector3d)
//
// returns reduced charVector3d from full 3d  vector reduced
// to  vector required by gp->giveStressStrainMode()
//
// enhaced method in order to support cases with integral bending (3dbeam)
// in such cases full strain vector has the form:
// 2) strainVectorShell {eps_x, gamma_xz, gamma_xy, \der{phi_x}{x}, kappa_y, kappa_z}

//
{
    MaterialMode mode = gp->giveMaterialMode();
    int size = charVector3d.giveSize();

    if ( mode == _1dFiber ) {
        if ( size != 6 ) {
            _error("giveReducedCharacteristicVector: stressVector3d size mismatch");
        }

        answer.resize(3);

        answer.at(1) = charVector3d.at(1);
        answer.at(2) = charVector3d.at(5);
        answer.at(3) = charVector3d.at(6);
    } else if ( mode == _3dBeam ) {
        if ( size != 6 ) {
            _error("giveReducedCharacteristicVector: stressVector3d size mismatch");
        }

        answer = charVector3d;
    } else {
        this->StructuralCrossSection :: giveReducedCharacteristicVector(answer, gp, charVector3d);
    }
}


void
FiberedCrossSection :: giveFullCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                    const FloatArray &charVector)
//
// returns reduced charVector3d from full 3d  vector reduced
// to  vector required by gp->giveStressStrainMode()
//
// enhaced method in order to support cases with integral bending (3dbeam..)
// in such cases full strain vector has the form:
// 2) strainVectorShell {eps_x, gamma_xz, gamma_xy, \der{phi_x}{x}, kappa_y, kappa_z}

//
{
    MaterialMode mode = gp->giveMaterialMode();
    int size = charVector.giveSize();

    if ( mode == _1dFiber ) {
        if ( size != 3 ) {
            _error("giveFullCharacteristicVector: stressVector size mismatch");
        }

        answer.resize(6);
        answer.zero();

        answer.at(1) = charVector.at(1);
        answer.at(5) = charVector.at(2);
        answer.at(6) = charVector.at(3);
    } else if ( mode == _3dBeam ) {
        if ( size != 6 ) {
            _error("giveFullCharacteristicVector: stressVector size mismatch");
        }

        answer = charVector;
    } else {
        this->StructuralCrossSection :: giveFullCharacteristicVector(answer, gp, charVector);
    }
}


FloatArray *
FiberedCrossSection :: imposeStressConstrainsOnGradient(GaussPoint *gp,
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
    case _1dFiber:
        for ( i = 2; i <= 4; i++ ) {
            gradientStressVector3d->at(i) = 0.;
        }

        break;
    default:
        StructuralCrossSection :: imposeStressConstrainsOnGradient(gp, gradientStressVector3d);
    }

    return gradientStressVector3d;
}


FloatArray *
FiberedCrossSection :: imposeStrainConstrainsOnGradient(GaussPoint *gp,
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
    case _1dFiber:
        for ( i = 2; i <= 4; i++ ) {
            gradientStrainVector3d->at(i) = 0.;
        }

        break;
    default:
        StructuralCrossSection :: imposeStrainConstrainsOnGradient(gp, gradientStrainVector3d);
    }

    return gradientStrainVector3d;
}


void
FiberedCrossSection :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
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
            case _1dFiber:
                answer.resize(3);
                answer.at(1) = 1;
                answer.at(2) = 5;
                answer.at(3) = 6;
                break;
            case _3dBeam:
                answer.resize(6);
                for ( i = 1; i <= 6; i++ ) {
                    answer.at(i) = i;
                }

                break;
            default:
                _error2( "giveStressStrainMask : unknown mode (%s)", __MaterialModeToString(mmode) );
            }
        } else if ( form == FullForm ) {
            switch ( mmode ) {
            case _1dFiber:
                answer.resize(6);
                answer.zero();
                answer.at(1) = 1;
                answer.at(5) = 2;
                answer.at(6) = 3;
                break;
            case _3dBeam:
                answer.resize(6);
                answer.zero();
                for ( i = 1; i <= 6; i++ ) {
                    answer.at(i) = i;
                }

                break;
            default:
                _error2( "giveStressStrainMask : unknown mode (%s)", __MaterialModeToString(mmode) );
            }
        } else {
            _error("giveStressStrainMask : unknown form mode");
        }

        return;
    }
}


int
FiberedCrossSection :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode, Material *mat)
{
    int i;

    if ( mmode == _3dBeam ) {
        if ( type == IST_BeamForceMomentumTensor ) {
            answer.resize(6);
            answer.zero();
            for ( i = 1; i <= 6; i++ ) {
                answer.at(i) = i;
            }

            return 1;
        }
    }

    return 0;
}


int
FiberedCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    StructuralMaterialStatus *status = ( StructuralMaterialStatus * ) ( aGaussPoint->giveMaterial()->giveStatus(aGaussPoint) );

    if ( type == IST_BeamForceMomentumTensor ) {
        answer = status->giveStressVector();
        return 1;
    } else if ( type == IST_BeamStrainCurvatureTensor ) {
        answer = status->giveStrainVector();
        return 1;
    }

    return 0;
}


IRResultType
FiberedCrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

#  ifdef VERBOSE
    // VERBOSE_PRINT1 ("Instanciating cross section ",this->giveNumber())
#  endif

    IR_GIVE_FIELD(ir, numberOfFibers, IFT_FiberedCrossSection_nfibers, "nfibers"); // Macro
    IR_GIVE_FIELD(ir, fiberMaterials, IFT_FiberedCrossSection_fibermaterials, "fibermaterials"); // Macro
    IR_GIVE_FIELD(ir, fiberThicks, IFT_FiberedCrossSection_thicks, "thicks"); // Macro
    IR_GIVE_FIELD(ir, fiberWidths, IFT_FiberedCrossSection_widths, "widths"); // Macro

    // read coordinates of fiber centers from (main pprincipal axes) mid-section
    IR_GIVE_FIELD(ir, fiberYcoords, IFT_FiberedCrossSection_fiberycentrecoords, "fiberycentrecoords"); // Macro
    IR_GIVE_FIELD(ir, fiberZcoords, IFT_FiberedCrossSection_fiberzcentrecoords, "fiberzcentrecoords"); // Macro

    IR_GIVE_FIELD(ir, thick, IFT_FiberedCrossSection_thick, "thick"); // Macro
    IR_GIVE_FIELD(ir, width, IFT_FiberedCrossSection_width, "width"); // Macro

    if ( ( numberOfFibers != fiberMaterials.giveSize() ) ||
        ( numberOfFibers != fiberThicks.giveSize() )    ||
        ( numberOfFibers != fiberWidths.giveSize() )     ||
        ( numberOfFibers != fiberYcoords.giveSize() )    ||
        ( numberOfFibers != fiberZcoords.giveSize() ) ) {
        _error("instanciateFrom : Array size mismatch ");
    }

    if ( numberOfFibers <= 0 ) {
        _error("instanciateFrom : numberOfFibers <= 0 is not alloved");
    }

    return IRRT_OK;
}

GaussPoint *
FiberedCrossSection :: giveSlaveGaussPoint(GaussPoint *masterGp, int i)
//
// return the i-th slave gauss point of master gp
// if slave gp don't exists - create them
//
{
    GaussPoint *slave = masterGp->giveSlaveGaussPoint(i);
    if ( slave == NULL ) {
        // check for proper dimensions - slave can be NULL if index too high or if not
        // slaves previously defined
        if ( i > this->numberOfFibers ) {
            _error("giveSlaveGaussPoint: no such fiber defined");
        }

        // create new slave record in masterGp
        // (requires that this is friend of gp)
        FloatArray *coords;
        // resolve slave material mode
        MaterialMode slaveMode, masterMode = masterGp->giveMaterialMode();
        slaveMode = this->giveCorrespondingSlaveMaterialMode(masterMode);

        masterGp->numberOfGp = this->numberOfFibers;
        masterGp->gaussPointArray = new GaussPoint * [ numberOfFibers ];

        for ( int j = 0; j < numberOfFibers; j++ ) {
            coords = new FloatArray(2);
            coords->at(1) = fiberYcoords.at(j + 1);
            coords->at(2) = fiberZcoords.at(j + 1);
            // in gp - is stored isoparametric coordinate (-1,1) of z-coordinate
            masterGp->gaussPointArray [ j ] = new GaussPoint(masterGp->giveIntegrationRule(), j + 1, coords, 0., slaveMode);
        }

        slave = masterGp->gaussPointArray [ i ];
    }

    return slave;
}


void
FiberedCrossSection :: printYourself()
// Prints the receiver on screen.
{
    printf("Cross Section with properties : \n");
    propertyDictionary->printYourself();
    printf("Fiber Materials: \n");
    fiberMaterials.printYourself();
    printf("Fiber Thicks   : \n");
    fiberThicks.printYourself();
    printf("Fiber Widths   : \n");
    fiberWidths.printYourself();
    printf("Fiber y coordinates: \n");
    fiberYcoords.printYourself();
    printf("Fiber y coordinates: \n");
    fiberZcoords.printYourself();
}


contextIOResultType
FiberedCrossSection :: saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *masterGp)
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

    StructuralMaterial *mat;
    GaussPoint *slaveGP;
    // and now save slave gp of master:
    for ( int i = 1; i <= numberOfFibers; i++ ) {
        slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1);
        mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at(i) ) );
        if ( ( iores = mat->saveIPContext(stream, mode, slaveGP) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType
FiberedCrossSection :: restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *masterGp)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
// restores also slaves of master gp
//
{
    contextIOResultType iores;

    if ( ( iores = CrossSection :: restoreIPContext(stream, mode, masterGp) ) != CIO_OK ) {
        THROW_CIOERR(iores);                                                                   // saved masterGp
    }

    // and now save slave gp of master:
    StructuralMaterial *mat;
    GaussPoint *slaveGP;
    for ( int i = 1; i <= numberOfFibers; i++ ) {
        // creates also slaves if they don't exists
        slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1);
        mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at(i) ) );
        if ( ( iores = mat->restoreIPContext(stream, mode, slaveGP) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


MaterialMode
FiberedCrossSection :: giveCorrespondingSlaveMaterialMode(MaterialMode masterMode)
//
// returns corresponding slave material mode to master mode
//
{
    if ( masterMode == _3dBeam ) {
        return _1dFiber;
    } else {
        _error("giveCorrespondingSlaveMaterialMode : unsupported mode");
    }

    return _Unknown;
}


FloatArray *
FiberedCrossSection :: GiveIntegrated3dBeamStress(GaussPoint *masterGp)
//
// computes integral internal forces for current mode
// this functions assumes that slave gp's have updatet stresses.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x, gamma_xz, gamma_xy, \der{phi_x}{x}, kappa_y, kappa_z}
//
{
    Material *fiberMat;
    StructuralMaterialStatus *fiberStatus;
    FloatArray *answer, fiberStress, reducedFiberStress;
    GaussPoint *fiberGp;
    double fiberThick, fiberWidth, fiberZCoord, fiberYCoord;
    int i;

    answer = new FloatArray(6);
    // perform integration over fibers

    for ( i = 1; i <= numberOfFibers; i++ ) {
        fiberGp = giveSlaveGaussPoint(masterGp, i - 1);
        fiberMat = domain->giveMaterial( fiberMaterials.at(i) );
        fiberStatus = ( ( StructuralMaterialStatus * ) fiberMat->giveStatus(fiberGp) );


        if ( fiberStatus->giveTempStressVector().giveSize() ) { // there exist total sress in gp
            reducedFiberStress = fiberStatus->giveTempStressVector();
            giveFullCharacteristicVector(fiberStress, fiberGp, reducedFiberStress);
        } else { // no total stress
            continue; // skip gp without stress
        }


        //
        // resolve current layer z-coordinate
        //
        fiberThick  = this->fiberThicks.at(i);
        fiberWidth  = this->fiberWidths.at(i);
        fiberYCoord = fiberGp->giveCoordinate(1);
        fiberZCoord = fiberGp->giveCoordinate(2);
        //
        // perform integration
        //
        // 1) membrane terms N, Qz, Qy
        answer->at(1) += fiberStress.at(1) * fiberWidth * fiberThick;
        answer->at(2) += fiberStress.at(5) * fiberWidth * fiberThick;
        answer->at(3) += fiberStress.at(6) * fiberWidth * fiberThick;
        // 2) bending terms mx, my, mxy
        answer->at(4) += ( fiberStress.at(5) * fiberWidth * fiberThick * fiberYCoord -
                           fiberStress.at(6) * fiberWidth * fiberThick * fiberZCoord );
        answer->at(5) += fiberStress.at(1) * fiberWidth * fiberThick * fiberZCoord;
        answer->at(6) -= fiberStress.at(1) * fiberWidth * fiberThick * fiberYCoord;
    }

    return answer;
}


double
FiberedCrossSection :: give(CrossSectionProperty aProperty)
{
    if ( aProperty == CS_Thickness ) {
        return this->thick;
    } else if ( aProperty == CS_Width )   {
        return this->width;
    } else if ( aProperty == CS_Area )   {
        return this->giveArea();
    }

    return CrossSection :: give(aProperty);
}


double FiberedCrossSection :: giveArea()
{
    int i;
    if ( this->area <= 0.0 ) {
        this->area = 0.0;
        for ( i = 1; i <= numberOfFibers; i++ ) {
            this->area += this->fiberThicks.at(i) * this->fiberWidths.at(i);
        }
    }

    return area;
}


void
FiberedCrossSection :: giveFiberMaterialStiffnessMatrix(FloatMatrix &fiberMatrix, MatResponseForm form,
                                                        MatResponseMode rMode, GaussPoint *layerGp,
                                                        TimeStep *tStep)
{
    this->giveMaterialStiffnessMatrixOf(fiberMatrix, form, rMode, layerGp,
                                        dynamic_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at( layerGp->giveNumber() ) ) ),
                                        tStep);
}


void
FiberedCrossSection :: computeStressIndependentStrainVector(FloatArray &answer,
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
