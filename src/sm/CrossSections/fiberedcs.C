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

#include "sm/CrossSections/fiberedcs.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "material.h"
#include "floatarray.h"
#include "verbose.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_CrossSection(FiberedCrossSection);

FloatArrayF<6>
FiberedCrossSection :: giveRealStress_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<6>();
}


FloatArrayF<4>
FiberedCrossSection :: giveRealStress_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<4>();
}


FloatArrayF<3>
FiberedCrossSection :: giveRealStress_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<3>();
}


FloatArrayF<1>
FiberedCrossSection :: giveRealStress_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<1>();
}


FloatArrayF<2>
FiberedCrossSection :: giveRealStress_Warping(const FloatArrayF<2> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported\n");
    return zeros<2>();
}


FloatMatrixF<6,6>
FiberedCrossSection :: giveStiffnessMatrix_3d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return FloatMatrixF<6,6>();
}


FloatMatrixF<3,3>
FiberedCrossSection :: giveStiffnessMatrix_PlaneStress(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return FloatMatrixF<3,3>();
}


FloatMatrixF<4,4>
FiberedCrossSection :: giveStiffnessMatrix_PlaneStrain(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return FloatMatrixF<4,4>();
}


FloatMatrixF<1,1>
FiberedCrossSection :: giveStiffnessMatrix_1d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return FloatMatrixF<1,1>();
}


FloatArrayF<3>
FiberedCrossSection :: giveGeneralizedStress_Beam2d(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<3>();
}


FloatArrayF<6>
FiberedCrossSection :: giveGeneralizedStress_Beam3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray fiberStrain;
    auto element = static_cast< StructuralElement * >( gp->giveElement() );
    auto interface = static_cast< FiberedCrossSectionInterface * >( element->giveInterface(FiberedCrossSectionInterfaceType) );

    if ( interface == nullptr ) {
        OOFEM_ERROR("element with no fiber support encountered");
    }

    FloatArrayF<6> answer;

    for ( int i = 1; i <= this->fiberMaterials.giveSize(); i++ ) {
        auto fiberGp = this->giveSlaveGaussPoint(gp, i - 1);
        auto fiberMat = static_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at(i) ) );
        // the question is whether this function should exist ?
        // if yes the element details will be hidden.
        // good idea also should be existence of element::GiveBmatrixOfLayer
        // and computing strains here - but first idea looks better
        // but treating of geometric non-linearities may become more complicated
        // another approach - use several functions with assumed kinematic constraints

        // resolve current layer z-coordinate
        double fiberThick  = this->fiberThicks.at(i);
        double fiberWidth  = this->fiberWidths.at(i);
        double fiberYCoord = fiberGp->giveNaturalCoordinate(1);
        double fiberZCoord = fiberGp->giveNaturalCoordinate(2);

        interface->FiberedCrossSectionInterface_computeStrainVectorInFiber(fiberStrain, strain, fiberGp, tStep);

        auto reducedFiberStress = fiberMat->giveRealStressVector_Fiber(fiberStrain, fiberGp, tStep);

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
    auto status = static_cast< StructuralMaterialStatus * >( domain->giveMaterial( fiberMaterials.at(1) )->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);

    return answer;
}


FloatArrayF<5>
FiberedCrossSection :: giveGeneralizedStress_Plate(const FloatArrayF<5> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<5>();
}


FloatArrayF<8>
FiberedCrossSection :: giveGeneralizedStress_Shell(const FloatArrayF<8> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<8>();
}

FloatArrayF<9>
FiberedCrossSection :: giveGeneralizedStress_ShellRot(const FloatArrayF<9> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<9>();
}


FloatArrayF<4>
FiberedCrossSection :: giveGeneralizedStress_MembraneRot(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported in given cross-section (yet).");
    return zeros<4>();
}

FloatArrayF<3>
FiberedCrossSection :: giveGeneralizedStress_PlateSubSoil(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported in given cross-section.");
    return zeros<3>();
}

void
FiberedCrossSection :: giveCharMaterialStiffnessMatrix(FloatMatrix &answer,
                                                       MatResponseMode rMode,
                                                       GaussPoint *gp,
                                                       TimeStep *tStep)
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
    } else {
        OOFEM_ERROR("Not implemented for bulk materials.");
    }
}


FloatMatrixF<3,3>
FiberedCrossSection :: give2dBeamStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not implemented");
    return FloatMatrixF<3,3>();
}


FloatMatrixF<6,6>
FiberedCrossSection :: give3dBeamStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
//
// General strain fiber vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x, gamma_xz, gamma_xy, \der{phi_x}{x}, kappa_y, kappa_z}
//
{
    double Ip = 0.0, A = 0.0, Ik, G = 0.0;

    FloatMatrixF<6,6> beamStiffness;

    // perform integration over layers
    for ( int i = 1; i <= this->fiberMaterials.giveSize(); i++ ) {
        auto fiberGp = giveSlaveGaussPoint(gp, i - 1);
        auto mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at( fiberGp->giveNumber() ) ) );
        auto fiberMatrix = mat->giveFiberStiffMtrx(rMode, fiberGp, tStep);
        //
        // resolve current layer z-coordinate
        //
        double fiberThick  = this->fiberThicks.at(i);
        double fiberWidth  = this->fiberWidths.at(i);
        double fiberZCoord = fiberZcoords.at(i);
        double fiberYCoord = fiberYcoords.at(i);
        double fiberYCoord2 = fiberYCoord * fiberYCoord;
        double fiberZCoord2 = fiberZCoord * fiberZCoord;
        //
        // perform integration
        //
        // 1) membrane terms N, Qz, Qy
        beamStiffness.at(1, 1) += fiberMatrix.at(1, 1) * fiberWidth * fiberThick;

        beamStiffness.at(2, 2) += fiberMatrix.at(2, 2) * fiberWidth * fiberThick;

        beamStiffness.at(3, 3) += fiberMatrix.at(3, 3) * fiberWidth * fiberThick;

        // 2) bending terms mx, my, mz

        Ip += fiberWidth * fiberThick * fiberZCoord2 + fiberWidth * fiberThick * fiberYCoord2;
        A  += fiberWidth * fiberThick;
        G  = fiberMatrix.at(2, 2) * fiberWidth * fiberThick;

        beamStiffness.at(5, 5) += fiberMatrix.at(1, 1) * fiberWidth * fiberThick * fiberZCoord2;
        beamStiffness.at(6, 6) += fiberMatrix.at(1, 1) * fiberWidth * fiberThick * fiberYCoord2;
    }

    ///@todo This must be wrong, it will use the last evaluated G (from the last fiber), outside the loop. FIXME!
    G /= A;
    Ik = A * A * A * A / ( 40.0 * Ip );
    beamStiffness.at(4, 4) = G * Ik;
    return beamStiffness;
}


FloatMatrixF<5,5>
FiberedCrossSection :: give2dPlateStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not implemented");
    return FloatMatrixF<5,5>();
}


FloatMatrixF<8,8>
FiberedCrossSection :: give3dShellStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not implemented");
    return FloatMatrixF<8,8>();
}

FloatMatrixF<9,9>
FiberedCrossSection :: give3dShellRotStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not implemented");
    return FloatMatrixF<9,9>();
}

FloatMatrixF<4,4>
FiberedCrossSection :: giveMembraneRotStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not implemented");
    return FloatMatrixF<4,4>();
}

FloatMatrixF<3,3>
FiberedCrossSection :: give2dPlateSubSoilStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
     OOFEM_ERROR("Not supported");
     return FloatMatrixF<3,3>();
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

    if ( type == IST_BeamForceMomentTensor ) {
        answer = status->giveStressVector();
        return 1;
    } else if ( type == IST_BeamStrainCurvatureTensor ) {
        answer = status->giveStrainVector();
        return 1;
    }
//     return CrossSection :: giveIPValue(answer, gp, type, tStep);
    ///@todo so far this only works for elements where each layer has its own integration rule
    int layer = gp->giveIntegrationRule()->giveNumber();
    return this->giveDomain()->giveMaterial( fiberMaterials.at(layer) )->giveIPValue(answer, gp, type, tStep);
}


void
FiberedCrossSection :: initializeFrom(InputRecord &ir)
{
#  ifdef VERBOSE
    // VERBOSE_PRINT1 ("Instanciating cross section ",this->giveNumber())
#  endif

     CrossSection :: initializeFrom(ir);
  
    IR_GIVE_FIELD(ir, fiberMaterials, _IFT_FiberedCrossSection_fibermaterials);
    IR_GIVE_FIELD(ir, fiberThicks, _IFT_FiberedCrossSection_thicks);
    IR_GIVE_FIELD(ir, fiberWidths, _IFT_FiberedCrossSection_widths);

    // read coordinates of fiber centers from (main pprincipal axes) mid-section
    IR_GIVE_FIELD(ir, fiberYcoords, _IFT_FiberedCrossSection_fiberycentrecoords);
    IR_GIVE_FIELD(ir, fiberZcoords, _IFT_FiberedCrossSection_fiberzcentrecoords);

    IR_GIVE_FIELD(ir, thick, _IFT_FiberedCrossSection_thick);
    IR_GIVE_FIELD(ir, width, _IFT_FiberedCrossSection_width);

    int num = fiberMaterials.giveSize();
    if ( num != fiberThicks.giveSize()    ||
         num != fiberWidths.giveSize()    ||
         num != fiberYcoords.giveSize()   ||
         num != fiberZcoords.giveSize() ) {
        throw ValueInputException(ir, _IFT_FiberedCrossSection_fibermaterials, "Array size mismatch ");
    }

    if ( num <= 0 ) {
        throw ValueInputException(ir, _IFT_FiberedCrossSection_fibermaterials, "number of fibers == 0 is not allowed");
    }

    area = fiberThicks.dotProduct(fiberWidths);
}

void FiberedCrossSection :: createMaterialStatus(GaussPoint &iGP)
{
    for ( int i = 1; i <= fiberMaterials.giveSize(); i++ ) {
        GaussPoint *fiberGp = this->giveSlaveGaussPoint(& iGP, i - 1);
        Material *mat = domain->giveMaterial( fiberMaterials.at(i) );
        MaterialStatus *matStat = mat->CreateStatus(fiberGp);
        iGP.setMaterialStatus(matStat);
    }
}

GaussPoint *
FiberedCrossSection :: giveSlaveGaussPoint(GaussPoint *masterGp, int i) const
//
// return the i-th slave gauss point of master gp
// if slave gp don't exists - create them
//
{
    auto slave = masterGp->giveSlaveGaussPoint(i);
    if ( !slave ) {
        // check for proper dimensions - slave can be NULL if index too high or if not
        // slaves previously defined
        if ( i > this->fiberMaterials.giveSize() ) {
            OOFEM_ERROR("no such fiber defined");
        }

        // create new slave record in masterGp
        // (requires that this is friend of gp)
        // resolve slave material mode
        auto masterMode = masterGp->giveMaterialMode();
        auto slaveMode = this->giveCorrespondingSlaveMaterialMode(masterMode);

        masterGp->gaussPoints.resize( fiberMaterials.giveSize() );

        for ( int j = 0; j < fiberMaterials.giveSize(); j++ ) {
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


void
FiberedCrossSection :: saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *masterGp)
{
    CrossSection :: saveIPContext(stream, mode, masterGp);

    // saved master gp record;
    // and now save slave gp of master:
    for ( int i = 1; i <= fiberMaterials.giveSize(); i++ ) {
        GaussPoint *slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1);
        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at(i) ) );
        mat->saveIPContext(stream, mode, slaveGP);
    }
}


void
FiberedCrossSection :: restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *masterGp)
{
    CrossSection :: restoreIPContext(stream, mode, masterGp);

    for ( int i = 1; i <= fiberMaterials.giveSize(); i++ ) {
        // creates also slaves if they don't exists
        GaussPoint *slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1);
        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at(i) ) );
        mat->restoreIPContext(stream, mode, slaveGP);
    }
}


MaterialMode
FiberedCrossSection :: giveCorrespondingSlaveMaterialMode(MaterialMode masterMode)
{
    if ( masterMode == _3dBeam ) {
        return _Fiber;
    } else {
        throw std::runtime_error("Unsupported mode");
    }

    return _Unknown;
}


double
FiberedCrossSection :: give(CrossSectionProperty aProperty, GaussPoint *gp) const
{
    if ( aProperty == CS_Thickness ) {
        return this->thick;
    } else if ( aProperty == CS_Width ) {
        return this->width;
    } else if ( aProperty == CS_Area ) { // not given in input
        return this->area;
    }

    return CrossSection :: give(aProperty, gp);
}


bool
FiberedCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode) const
{
    ///@todo As far as I can see, it only uses diagonal components for the 3dbeam, but there is no way to check here.

    for ( int imat : this->fiberMaterials ) {
        if ( !this->domain->giveMaterial( imat )->isCharacteristicMtrxSymmetric(rMode) ) {
            return false;
        }
    }
    return true;
}


int
FiberedCrossSection :: checkConsistency()
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

Material *
FiberedCrossSection :: giveMaterial(IntegrationPoint *ip) const
{
    ///@todo We should keep track in integration point (integration rule) what material from layer is assigned. Otherwise difficulties due to different elements and IP numbering.
    if ( ip->giveIntegrationRule()->giveIntegrationDomain() == _Cube ||
        ip->giveIntegrationRule()->giveIntegrationDomain() == _Wedge
    ) {
        return domain->giveMaterial( fiberMaterials.at(1) );
        //return this->domain->giveMaterial( this->giveLayerMaterial(ip->giveNumber()) );
    }
    
    if (ip->hasSlaveGaussPoint()) {
        return domain->giveMaterial( fiberMaterials.at(1) );//virtual master, has no material assigned in input file
    } else {
        return domain->giveMaterial( fiberMaterials.at(1) );//virtual master, has no material assigned in input file
        //OOFEM_ERROR("Not implemented.")
    }
    return nullptr;
}
} // end namespace oofem
