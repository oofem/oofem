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

#include "../sm/CrossSections/fiberedcs.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "material.h"
#include "floatarray.h"
#include "verbose.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_CrossSection(FiberedCrossSection);

void
FiberedCrossSection :: giveRealStress_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
FiberedCrossSection :: giveRealStress_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
FiberedCrossSection :: giveRealStress_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
FiberedCrossSection :: giveRealStress_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
FiberedCrossSection :: giveRealStress_Warping(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported\n");
}


void
FiberedCrossSection :: giveStiffnessMatrix_3d(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
FiberedCrossSection :: giveStiffnessMatrix_PlaneStress(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
FiberedCrossSection :: giveStiffnessMatrix_PlaneStrain(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
FiberedCrossSection :: giveStiffnessMatrix_1d(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}




void
FiberedCrossSection :: giveGeneralizedStress_Beam2d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
FiberedCrossSection :: giveGeneralizedStress_Beam3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    double fiberThick, fiberWidth, fiberZCoord, fiberYCoord;
    FloatArray fiberStrain, reducedFiberStress;
    StructuralElement *element = static_cast< StructuralElement * >( gp->giveElement() );
    FiberedCrossSectionInterface *interface;

    if ( ( interface = static_cast< FiberedCrossSectionInterface * >( element->giveInterface(FiberedCrossSectionInterfaceType) ) ) == NULL ) {
        OOFEM_ERROR("element with no fiber support encountered");
    }

    answer.resize(6);
    answer.zero();

    for ( int i = 1; i <= numberOfFibers; i++ ) {
        GaussPoint *fiberGp = this->giveSlaveGaussPoint(gp, i - 1);
        StructuralMaterial *fiberMat = static_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at(i) ) );
        // the question is whether this function should exist ?
        // if yes the element details will be hidden.
        // good idea also should be existence of element::GiveBmatrixOfLayer
        // and computing strains here - but first idea looks better
        // but treating of geometric non-linearities may become more complicated
        // another approach - use several functions with assumed kinematic constraints

        // resolve current layer z-coordinate
        fiberThick  = this->fiberThicks.at(i);
        fiberWidth  = this->fiberWidths.at(i);
        fiberYCoord = fiberGp->giveNaturalCoordinate(1);
        fiberZCoord = fiberGp->giveNaturalCoordinate(2);

        interface->FiberedCrossSectionInterface_computeStrainVectorInFiber(fiberStrain, strain, fiberGp, tStep);

        fiberMat->giveRealStressVector_Fiber(reducedFiberStress, fiberGp, fiberStrain, tStep);

        // perform integration
        // 1) membrane terms N, Qz, Qy
        answer.at(1) += reducedFiberStress.at(1) * fiberWidth * fiberThick;
        answer.at(2) += reducedFiberStress.at(2) * fiberWidth * fiberThick;
        answer.at(3) += reducedFiberStress.at(3) * fiberWidth * fiberThick;
        // 2) bending terms mx, my, mxy
        answer.at(4) += ( reducedFiberStress.at(2) * fiberWidth * fiberThick * fiberYCoord -
                          reducedFiberStress.at(3) * fiberWidth * fiberThick * fiberZCoord );
        answer.at(5) += reducedFiberStress.at(1) * fiberWidth * fiberThick * fiberZCoord;
        answer.at(6) -= reducedFiberStress.at(1) * fiberWidth * fiberThick * fiberYCoord;
    }

    // now we must update master gp ///@ todo simply chosen the first fiber material as master material /JB
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >
                                       ( domain->giveMaterial( fiberMaterials.at(1) )->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


void
FiberedCrossSection :: giveGeneralizedStress_Plate(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
FiberedCrossSection :: giveGeneralizedStress_Shell(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported");
}


void
FiberedCrossSection :: giveGeneralizedStress_MembraneRot(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
    OOFEM_ERROR("Not supported in given cross-section (yet).");
}

void
FiberedCrossSection :: giveGeneralizedStress_PlateSubSoil(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep)
{
  OOFEM_ERROR("Not supported in given cross-section.");
}

void
FiberedCrossSection :: giveCharMaterialStiffnessMatrix(FloatMatrix &answer,
                                                       MatResponseMode rMode,
                                                       GaussPoint *gp,
                                                       TimeStep *tStep)
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
        OOFEM_ERROR("Not implemented for bulk materials.");
    }
}


void
FiberedCrossSection :: give2dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
}


void
FiberedCrossSection :: give3dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// General strain fiber vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x, gamma_xz, gamma_xy, \der{phi_x}{x}, kappa_y, kappa_z}
//
{
    FloatMatrix fiberMatrix;
    GaussPoint *fiberGp;
    double fiberThick, fiberWidth, fiberZCoord, fiberYCoord;
    double fiberZCoord2, fiberYCoord2, Ip = 0.0, A = 0.0, Ik, G = 0.0;

    // if (form != ReducedForm) error ("give3dShellMaterialStiffness : full form unsupported");

    answer.resize(6, 6);
    answer.zero();
    // perform integration over layers

    for ( int i = 1; i <= numberOfFibers; i++ ) {
        fiberGp = giveSlaveGaussPoint(gp, i - 1);
        this->giveFiberMaterialStiffnessMatrix(fiberMatrix, rMode, fiberGp, tStep);
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

        answer.at(2, 2) += fiberMatrix.at(2, 2) * fiberWidth * fiberThick;

        answer.at(3, 3) += fiberMatrix.at(3, 3) * fiberWidth * fiberThick;

        // 2) bending terms mx, my, mz

        Ip             += fiberWidth * fiberThick * fiberZCoord2 + fiberWidth * fiberThick * fiberYCoord2;
        A              += fiberWidth * fiberThick;
        G               = fiberMatrix.at(2, 2) * fiberWidth * fiberThick;

        answer.at(5, 5) += fiberMatrix.at(1, 1) * fiberWidth * fiberThick * fiberZCoord2;
        answer.at(6, 6) += fiberMatrix.at(1, 1) * fiberWidth * fiberThick * fiberYCoord2;
    }

    ///@todo This must be wrong, it will use the last evaluated G (from the last fiber), outside the loop. FIXME!
    G /= A;
    Ik = A * A * A * A / ( 40.0 * Ip );
    answer.at(4, 4) = G * Ik;
}


void
FiberedCrossSection :: give2dPlateStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
}


void
FiberedCrossSection :: give3dShellStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
}

void
FiberedCrossSection :: giveMembraneRotStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
}

void
FiberedCrossSection :: give2dPlateSubSoilStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
     OOFEM_ERROR("Not supported");
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
    int size = gradientStressVector3d->giveSize();
    if ( size != 6 ) {
        OOFEM_ERROR("gradientStressVector3d size mismatch");
    }

    switch ( mode ) {
    case _Fiber:
        for ( int i = 2; i <= 4; i++ ) {
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
    int size = gradientStrainVector3d->giveSize();
    if ( size != 6 ) {
        OOFEM_ERROR("gradientStrainVector3d size mismatch");
    }

    switch ( mode ) {
    case _Fiber:
        for ( int i = 2; i <= 4; i++ ) {
            gradientStrainVector3d->at(i) = 0.;
        }

        break;
    default:
        StructuralCrossSection :: imposeStrainConstrainsOnGradient(gp, gradientStrainVector3d);
    }

    return gradientStrainVector3d;
}


int
FiberedCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    Material *mat = this->giveDomain()->giveMaterial( fiberMaterials.at(1) ); ///@todo For now, create material status according to the first fiber material
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) );

    if ( type == IST_BeamForceMomentumTensor ) {
        answer = status->giveStressVector();
        return 1;
    } else if ( type == IST_BeamStrainCurvatureTensor ) {
        answer = status->giveStrainVector();
        return 1;
    }
    return CrossSection :: giveIPValue(answer, gp, type, tStep);
}


IRResultType
FiberedCrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

#  ifdef VERBOSE
    // VERBOSE_PRINT1 ("Instanciating cross section ",this->giveNumber())
#  endif

    IR_GIVE_FIELD(ir, numberOfFibers, _IFT_FiberedCrossSection_nfibers);
    IR_GIVE_FIELD(ir, fiberMaterials, _IFT_FiberedCrossSection_fibermaterials);
    IR_GIVE_FIELD(ir, fiberThicks, _IFT_FiberedCrossSection_thicks);
    IR_GIVE_FIELD(ir, fiberWidths, _IFT_FiberedCrossSection_widths);

    // read coordinates of fiber centers from (main pprincipal axes) mid-section
    IR_GIVE_FIELD(ir, fiberYcoords, _IFT_FiberedCrossSection_fiberycentrecoords);
    IR_GIVE_FIELD(ir, fiberZcoords, _IFT_FiberedCrossSection_fiberzcentrecoords);

    IR_GIVE_FIELD(ir, thick, _IFT_FiberedCrossSection_thick);
    IR_GIVE_FIELD(ir, width, _IFT_FiberedCrossSection_width);

    if ( numberOfFibers != fiberMaterials.giveSize() ||
         numberOfFibers != fiberThicks.giveSize()    ||
         numberOfFibers != fiberWidths.giveSize()    ||
         numberOfFibers != fiberYcoords.giveSize()   ||
         numberOfFibers != fiberZcoords.giveSize() ) {
        OOFEM_WARNING("Array size mismatch ");
        return IRRT_BAD_FORMAT;
    }

    if ( numberOfFibers <= 0 ) {
        OOFEM_WARNING("numberOfFibers <= 0 is not allowed");
        return IRRT_BAD_FORMAT;
    }

    return IRRT_OK;
}

void FiberedCrossSection :: createMaterialStatus(GaussPoint &iGP)
{
    for ( int i = 1; i <= numberOfFibers; i++ ) {
        GaussPoint *fiberGp = this->giveSlaveGaussPoint(& iGP, i - 1);
        StructuralMaterial *mat = static_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at(i) ) );
        MaterialStatus *matStat = mat->CreateStatus(fiberGp);
        iGP.setMaterialStatus(matStat);
    }
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
            OOFEM_ERROR("no such fiber defined");
        }

        // create new slave record in masterGp
        // (requires that this is friend of gp)
        // resolve slave material mode
        MaterialMode slaveMode, masterMode = masterGp->giveMaterialMode();
        slaveMode = this->giveCorrespondingSlaveMaterialMode(masterMode);

        masterGp->gaussPoints.resize( this->numberOfFibers );

        for ( int j = 0; j < numberOfFibers; j++ ) {
            // in gp - is stored isoparametric coordinate (-1,1) of z-coordinate
            masterGp->gaussPoints [ j ] = new GaussPoint(masterGp->giveIntegrationRule(), j + 1,
                                                {fiberYcoords.at(j + 1), fiberZcoords.at(j + 1)}, 0., slaveMode);
        }

        slave = masterGp->gaussPoints [ i ];
    }

    return slave;
}


void
FiberedCrossSection :: printYourself()
// Prints the receiver on screen.
{
    printf("Cross Section with properties : \n");
    propertyDictionary.printYourself();
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
FiberedCrossSection :: saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *masterGp)
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
    for ( int i = 1; i <= numberOfFibers; i++ ) {
        GaussPoint *slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1);
        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at(i) ) );
        if ( ( iores = mat->saveIPContext(stream, mode, slaveGP) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType
FiberedCrossSection :: restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *masterGp)
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

    for ( int i = 1; i <= numberOfFibers; i++ ) {
        // creates also slaves if they don't exists
        GaussPoint *slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1);
        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at(i) ) );
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
        return _Fiber;
    } else {
        OOFEM_ERROR("unsupported mode");
    }

    return _Unknown;
}


double
FiberedCrossSection :: give(CrossSectionProperty aProperty, GaussPoint *gp)
{
    if ( aProperty == CS_Thickness ) {
        return this->thick;
    } else if ( aProperty == CS_Width ) {
        return this->width;
    } else if ( aProperty == CS_Area ) { // not given in input
        return this->giveArea();
    }

    return CrossSection :: give(aProperty, gp);
}


double FiberedCrossSection :: giveArea()
{
    if ( this->area <= 0.0 ) {
        this->area = 0.0;
        for ( int i = 1; i <= numberOfFibers; i++ ) {
            this->area += this->fiberThicks.at(i) * this->fiberWidths.at(i);
        }
    }

    return area;
}


void
FiberedCrossSection :: giveFiberMaterialStiffnessMatrix(FloatMatrix &fiberMatrix,
                                                        MatResponseMode rMode, GaussPoint *layerGp,
                                                        TimeStep *tStep)
{
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at( layerGp->giveNumber() ) ) );
    mat->giveStiffnessMatrix(fiberMatrix, rMode, layerGp, tStep);
}


bool FiberedCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode)
{
    ///@todo As far as I can see, it only uses diagonal components for the 3dbeam, but there is no way to check here.

    for ( int i = 1; i <= this->numberOfFibers; i++ ) {
        if ( !this->domain->giveMaterial( this->fiberMaterials.at(i) )->isCharacteristicMtrxSymmetric(rMode) ) {
            return false;
        }
    }
    return true;
}


int
FiberedCrossSection :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    for ( int i = 1; this->fiberMaterials.giveSize(); i++ ) {
        Material *mat = this->giveDomain()->giveMaterial( this->fiberMaterials.at(i) );
        if ( !dynamic_cast< StructuralMaterial * >(mat) ) {
            OOFEM_WARNING("material %s without structural support", mat->giveClassName() );
            result = 0;
            continue;
        }
    }
    return result;
}
} // end namespace oofem
