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

#include "intmatcoulombcontact.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(IntMatCoulombContact);

IntMatCoulombContact :: IntMatCoulombContact(int n, Domain *d) : StructuralInterfaceMaterial(n, d) { }

IntMatCoulombContact :: ~IntMatCoulombContact() { }


void
IntMatCoulombContact :: giveEngTraction_3d( FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
{
    // Returns the (engineering) traction vector in 3d based on the spatial jump.

    IntMatCoulombContactStatus *status = static_cast< IntMatCoulombContactStatus * >( this->giveStatus( gp ) );
    this->initGpForNewStep( gp );

    double normalJump = jump.at( 3 );
    FloatArray shearJump;
    shearJump.setValues( 2, jump.at( 1 ), jump.at( 2 ) );

    double normalStress = 0.0;
    FloatArray shearStress, tempShearStressShift;
    tempShearStressShift = status->giveTempShearStressShift();
    this->computeEngTraction( normalStress, shearStress, tempShearStressShift,
                              normalJump, shearJump );
    
    // Set stress components in the traction vector
    answer.resize( 3 );
    answer.at( 1 ) = shearStress.at( 1 );
    answer.at( 2 ) = shearStress.at( 2 );
    answer.at( 3 ) = normalStress;

    // Update gp
    status->setTempShearStressShift( tempShearStressShift );
    status->letTempJumpBe( jump );
    status->letTempTractionBe( answer );

#if 0
    IntMatCoulombContactStatus *status = static_cast< IntMatCoulombContactStatus * >( this->giveStatus(gp) );
    this->initGpForNewStep(gp);
    FloatArray shearJump(2), shearStress; 

    FloatArray tempShearStressShift = status->giveTempShearStressShift();
    double normalJump = jump.at( jump.giveSize() );
    
    double normalStress, maxShearStress, dp;
    double shift = -this->kn * this->stiffCoeff * normalClearance;

    answer.zero();
    if ( normalJump + normalClearance <= 0. ) {
        normalStress = this->kn * ( normalJump + normalClearance ) + shift; //in compression and after the clearance gap closed
        maxShearStress = fabs(normalStress) * this->frictCoeff;
    } else {
        normalStress = this->kn * this->stiffCoeff * ( normalJump + normalClearance ) + shift;
        maxShearStress = 0.;
    }

    //switch ( mMode ) {
    switch( jump.giveSize() ) {
    //case _1dInterface:
    case 1:
        answer.resize(1);
        break;
    //case _2dInterface:
    case 2:
        answer.resize(2);
        shearJump.at(1) = jump.at(2);
        shearStress.beScaled(this->kn, shearJump);
        shearStress.subtract(tempShearStressShift);
        dp = shearStress.dotProduct(shearStress, 1);
        if ( dp > maxShearStress * maxShearStress ) {
            shearStress.times( maxShearStress / sqrt(dp) );
        }

        tempShearStressShift.beScaled(this->kn, shearJump);
        tempShearStressShift.subtract(shearStress);
        //answer.at(2) = shearStress.at(1);
        answer.at( 1 ) = shearStress.at( 1 );
        break;
    //case _3dInterface:
    case 3:
        answer.resize(3);
        //shearJump.at(1) = jump.at(2);
        //shearJump.at(2) = jump.at(3);
        shearJump.at( 1 ) = jump.at( 1 );
        shearJump.at( 2 ) = jump.at( 2 );
        shearStress.beScaled(this->kn, shearJump);
        shearStress.subtract(tempShearStressShift);
        dp = shearStress.dotProduct(shearStress, 2);
        if ( dp > maxShearStress * maxShearStress ) {
            shearStress.times( maxShearStress / sqrt(dp) );
        }

        tempShearStressShift.beScaled(this->kn, shearJump);
        tempShearStressShift.subtract(shearStress);
        //answer.at(2) = shearStress.at(1);
        //answer.at(3) = shearStress.at(2);
        answer.at( 1 ) = shearStress.at( 1 );
        answer.at( 2 ) = shearStress.at( 2 );
        break;
    default:
        _error("giveMaterialMode: Unsupported interface mode");
    }

    double lim = 1.e+50;
    //answer.at(1) = min(normalStress, lim);  //threshold on maximum
    //answer.at(1) = max(answer.at(1), -lim);  //threshold on minimum
    answer.at( jump.giveSize() ) = min( normalStress, lim );  //threshold on maximum
    answer.at( jump.giveSize() ) = max( answer.at( 1 ), -lim );  //threshold on minimum

    //answer.at(1) = normalStress > lim ? lim : normalStress < -lim ? -lim : normalStress;

    // update gp
    status->setTempShearStressShift(tempShearStressShift);
    //status->letTempStrainVectorBe(strainVector);
    //status->letTempStressVectorBe(answer);

    status->letTempJumpBe( jump );
    status->letTempTractionBe( answer );
#endif

}




void
IntMatCoulombContact::giveEngTraction_2d( FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep )
{
    // Returns the (engineering) traction vector in 2d based on the spatial jump.

    IntMatCoulombContactStatus *status = static_cast< IntMatCoulombContactStatus * >( this->giveStatus( gp ) );
    this->initGpForNewStep( gp );
    double normalJump = jump.at( 2 );
    FloatArray shearJump;
    shearJump.setValues( 1, jump.at( 1 ) );

    double normalStress = 0.0;
    FloatArray shearStress, tempShearStressShift;
    tempShearStressShift = status->giveTempShearStressShift();
    this->computeEngTraction( normalStress, shearStress, tempShearStressShift,
                              normalJump, shearJump );

    // Set stress components in the traction vector
    answer.resize( 2 );
    answer.at( 1 ) = shearStress.at( 1 );
    answer.at( 2 ) = normalStress;

    // Update gp
    status->setTempShearStressShift( tempShearStressShift );
    status->letTempJumpBe( jump );
    status->letTempTractionBe( answer );

}


void
IntMatCoulombContact::giveEngTraction_1d( FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep )
{
    // Returns the (engineering) traction vector (normal stress only) in 1d based on the 
    // spatial jump. The shear stress is not relevant in this case.

    IntMatCoulombContactStatus *status = static_cast< IntMatCoulombContactStatus * >( this->giveStatus( gp ) );
    this->initGpForNewStep( gp );

    double normalJump = jump.at( 1 );
    FloatArray shearJump(0);

    double normalStress = 0.0;
    FloatArray shearStress, tempShearStressShift(0);
    this->computeEngTraction( normalStress, shearStress, tempShearStressShift,
                              normalJump, shearJump );

    // Set stress components in the traction vector
    answer.resize( 1 );
    answer.at( 1 ) = normalStress;

    // Update gp
    status->letTempJumpBe( jump );
    status->letTempTractionBe( answer );
}


void 
IntMatCoulombContact::computeEngTraction( double &normalStress, FloatArray &shearStress,
FloatArray &tempShearStressShift, const double normalJump, const FloatArray &shearJump )
{

    double maxShearStress = 0.0;
    double shift = -this->kn * this->stiffCoeff * normalClearance;

    if(normalJump + normalClearance <= 0.) {
        normalStress = this->kn * ( normalJump + normalClearance ) + shift; //in compression and after the clearance gap closed
        maxShearStress = fabs( normalStress ) * this->frictCoeff;
    }
    else {
        normalStress = this->kn * this->stiffCoeff * ( normalJump + normalClearance ) + shift;
        maxShearStress = 0.;
    }

    if(shearJump.giveSize() != 0) {

        shearStress = this->kn * shearJump - tempShearStressShift;
        double dp = shearStress.dotProduct( shearStress );
        const double eps = 1.0e-15; // small number
        if(dp > maxShearStress * maxShearStress) {
            shearStress.times( maxShearStress / ( sqrt( dp ) + eps ) );
        }
        tempShearStressShift = this->kn * shearJump - shearStress;

    } else { // 1d -> no shear stresses
        return;
    }
}


void
IntMatCoulombContact :: giveGeneralStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode,
                                               GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{

    IntMatCoulombContactStatus *status = static_cast< IntMatCoulombContactStatus * > ( this->giveStatus(gp) );
    FloatArray jump;
    jump = status->giveTempJump();

    const int numSpaceDim = jump.giveSize();
    const double normalJump = jump.at( numSpaceDim );

    answer.resize( numSpaceDim, numSpaceDim );
    answer.beUnitMatrix();
    if(rMode == SecantStiffness || rMode == TangentStiffness) {
        if( normalJump + normalClearance <= 0) {
            answer.times( this->kn ); //in compression and after the clearance gap closed
        } else {
            answer.times( this->kn * this->stiffCoeff );
        }
    } else {
        answer.times( this->kn );
    }

#if 0
    switch ( jump.giveSize() ) {
    //case _1dInterface:
    case 1:
        answer.resize(1, 1);
        if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
            //if ( normalJump + normalClearance <= 0 ) {
            if( jump.at(1) + normalClearance <= 0) {
                answer.at(1, 1) = this->kn; //in compression and after the clearance gap closed
            } else {
                answer.at(1, 1) = this->kn * this->stiffCoeff;
            }
        } else {
            if ( rMode == ElasticStiffness ) {
                answer.at(1, 1) = this->kn;
            } else {
                _error2( "give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode (%s)", __MatResponseModeToString(rMode) );
            }
        }

        return;

    //case _2dInterface:
    case 2:
        answer.resize(2, 2);
        if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
            //if ( normalJump + normalClearance <= 0. ) {
            if( jump.at(2) + normalClearance <= 0.) {
                answer.at(1, 1) = answer.at(2, 2) = this->kn; //in compression and after the clearance gap closed
            } else {
                answer.at(1, 1) = answer.at(2, 2) = this->kn * this->stiffCoeff;
            }
        } else {
            if ( rMode == ElasticStiffness ) {
                answer.at(1, 1) = answer.at(2, 2) = this->kn;
            } else {
                _error2( "give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode (%s)", __MatResponseModeToString(rMode) );
            }
        }

        return;

    //case _3dInterface:
    case 3:
        answer.resize(3, 3);
        if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
            //if ( normalJump + normalClearance <= 0. ) {
            if( jump.at(3) + normalClearance <= 0.) {
                answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = this->kn; //in compression and after the clearance gap closed
            } else {
                answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = this->kn * this->stiffCoeff;
            }
        } else {
            if ( rMode == ElasticStiffness ) {
                answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = this->kn;
            } else {
                _error2( "give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode (%s)", __MatResponseModeToString(rMode) );
            }
        }

        return;

    default:
        //StructuralMaterial :: giveStiffnessMatrix(answer, rMode, gp, tStep);
        OOFEM_ERROR1( "IntMatCoulombContact :: giveStiffnessMatrix - size of jump is not 1, 2 or 3");
    }
#endif
}




IRResultType
IntMatCoulombContact :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    frictCoeff = 0.;
    stiffCoeff = 0.;
    normalClearance = 0.;
    IR_GIVE_FIELD(ir, kn, _IFT_IntMatCoulombContact_kn);
    IR_GIVE_OPTIONAL_FIELD(ir, frictCoeff, _IFT_IntMatCoulombContact_frictCoeff);
    IR_GIVE_OPTIONAL_FIELD(ir, stiffCoeff, _IFT_IntMatCoulombContact_stiffCoeff);
    IR_GIVE_OPTIONAL_FIELD(ir, normalClearance, _IFT_IntMatCoulombContact_normalClearance);

    return StructuralInterfaceMaterial :: initializeFrom( ir );
}


void
IntMatCoulombContact :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);
    input.setField(this->kn, _IFT_IntMatCoulombContact_kn);
    input.setField(this->frictCoeff, _IFT_IntMatCoulombContact_frictCoeff);
    input.setField(this->stiffCoeff, _IFT_IntMatCoulombContact_stiffCoeff);
    input.setField(this->normalClearance, _IFT_IntMatCoulombContact_normalClearance);
}


IntMatCoulombContactStatus :: IntMatCoulombContactStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
{
    int size = d->giveNumberOfSpatialDimensions() - 1;
    shearStressShift.resize(size);
    tempShearStressShift.resize(size);
    shearStressShift.zero();
    tempShearStressShift.zero();
    //this->initTempStatus();
}


IntMatCoulombContactStatus :: ~IntMatCoulombContactStatus()
{ }


void
IntMatCoulombContactStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus::printOutputAt( file, tStep );
    fprintf(file, "status { ");
    if(this->shearStressShift.giveSize() == 2) {
        fprintf( file, "shearStressShift (%f, %f)", this->shearStressShift.at( 1 ), this->shearStressShift.at( 2 ) );
    } else if(this->shearStressShift.giveSize() == 1) {
        fprintf( file, "shearStressShift (%f)", this->shearStressShift.at( 1 ) );
    }
    fprintf(file, "}\n");
}


void
IntMatCoulombContactStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus::initTempStatus( );
    tempShearStressShift = shearStressShift;
}


void
IntMatCoulombContactStatus :: updateYourself(TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus::updateYourself( tStep );
    shearStressShift = tempShearStressShift;
}


FloatArray
IntMatCoulombContactStatus :: giveTempShearStressShift()
{
    FloatArray answer = tempShearStressShift;
    return answer;
}


contextIOResultType
IntMatCoulombContactStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if(( iores = StructuralInterfaceMaterialStatus::saveContext( stream, mode, obj ) ) != CIO_OK) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    //if ( !stream->write(& kappa, 1) ) {
    //THROW_CIOERR(CIO_IOERR);
    //}

    return CIO_OK;
}


contextIOResultType
IntMatCoulombContactStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if(( iores = StructuralInterfaceMaterialStatus::restoreContext( stream, mode, obj ) ) != CIO_OK) {
        THROW_CIOERR(iores);
    }

    // read raw data
    //if ( !stream->read(& kappa, 1) ) {
    //THROW_CIOERR(CIO_IOERR);
    //}

    return CIO_OK;
}
} // end namespace oofem
